#!/usr/bin/python

""" Implement Distance Method for HGT detection based on algorithm described
    in:
        Wei. X et al., "A Distance-Based Method for Detecting HGT in Whole Genomes",
        International Symposium on Bioinformatics Research and Applications (ISBRA),
        2008, pages 26-37

    The workflow follows the algorithm:
    1. For each gene in target genome,
        i.    BLAST sequence against all other genes in the reference genomes;
        ii.   Go to step 3 if gene has more than threshold number of homologs
              (min-num-homologs), otherwise go to next gene in target genome;
        iii.  Compute multiple sequence alignment on homolog genes using CLUSTAL;
        iv.   Compute pairwise distance matrix using PHYLIP's protdist function
              and Z-score normalize the set of pairwise distances for each gene
              family and species;
        v.    Add distance matrix for all pairwise distances into global distance
              matrix storing results for all genes

    2. Cluster gene families by species,
        vi.   Compute all species sets (sets of genes whose orthologs are detectable
              in exactly the same subset of the considered species) using the
              Hamming distance clustering algorithm;
        vii.  Cluster genes to each core species set cluster;
        viii. Run outlier detection algorithm on each cluster (paragraph 2
              of section 'Detecting Outlier Genes' in original paper)
"""

import sys
import click
import glob
import numpy
import operator
import threading
import subprocess
import traceback
import shlex
from os.path import join, splitext, basename
from string import maketrans
from itertools import imap, chain


from skbio.parse.sequences import parse_fasta


class Command(object):
    """
    Enables to run subprocess commands in a different thread
    with TIMEOUT option.

    Based on jcollado's solution:
    http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
    https://gist.github.com/kirpit/1306188
    """
    command = None
    process = None
    status = None
    output, error = '', ''
 
    def __init__(self, command):
        if isinstance(command, basestring):
            command = shlex.split(command)
        self.command = command
 
    def run(self, timeout=None, **kwargs):
        """ Run a command then return: (status, output, error). """
        def target(**kwargs):
            try:
                self.process = subprocess.Popen(self.command, **kwargs)
                self.output, self.error = self.process.communicate()
                self.status = self.process.returncode
            except:
                self.error = traceback.format_exc()
                self.status = -1
        # default stdout and stderr
        if 'stdout' not in kwargs:
            kwargs['stdout'] = subprocess.PIPE
        if 'stderr' not in kwargs:
            kwargs['stderr'] = subprocess.PIPE
        # thread
        thread = threading.Thread(target=target, kwargs=kwargs)
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
        return self.status, self.output, self.error


def hamming(str1, str2):
    "Compute the Hamming distance between two strings"
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(imap(ne, str1, str2))


def preprocess_data(working_dir,
                    target_proteomes_dir,
                    verbose=False):
    """ Map each gene to sudo name (ex. S1_G1) for easier
        output comparison and the 10 character name limitation
        in PHYLIP output. This format is limited up to 999
        species and 9999 genes per species.
    """
    gene_map = {}
    ref_db = {}
    if verbose:
        sys.stdout.write("Target organism\tNumber of genes\n")
    # each file contains genes for species
    species = 0
    for _file in glob.glob(join(target_proteomes_dir, "*.fa")):
        if verbose:
            sys.stdout.write("%s. %s\t" % (species+1, basename(_file)))
        with open(_file, 'rb') as readfile:
            gene = 0
            for label, seq in parse_fasta(readfile):
                label = label.split()[0]
                # Blast parses reference ids to remove the |cl| from <|cl|ref_id
                # (update this for other parsing requirements)
                #label = label.split('|')[1]
                ref_db[label] = seq
                sudo_label = "%s_%s" % (species, gene)
                gene_map[label] = sudo_label
                gene_map[sudo_label] = label
                gene+=1
        if verbose:
            sys.stdout.write("%s\n" % (gene+1))
        species+=1
    return gene_map, ref_db, species
            

def launch_blast(query_proteome_fp,
                 ref_fp,
                 working_dir,
                 e_value=10e-20,
                 threads=1):
    """ Launch BLASTp given a query proteome and a reference database of proteomes
    """
    # build blast database
    makeblastdb_command = ["makeblastdb",
                           "-in", ref_fp,
                           "-dbtype", "prot",
                           "-parse_seqids"]
    proc = subprocess.Popen(makeblastdb_command,
                 stdout=subprocess.PIPE,
                 stderr=subprocess.PIPE,
                 close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr

    # launch blast
    out_file_fp = join(working_dir, "%s.blast" % basename(splitext(query_proteome_fp)[0]))
    blastp_command = ["blastp",
                      "-db", ref_fp,
                      "-query", query_proteome_fp,
                      "-evalue", str(e_value),
                      "-num_threads", str(threads),
                      "-outfmt", "6 std qcovs",
                      "-task", "blastp",
                      "-out", out_file_fp]
    proc = subprocess.Popen(blastp_command,
                 stdout=subprocess.PIPE,
                 stderr=subprocess.PIPE,
                 close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr

    return out_file_fp


def parse_blast(alignments_fp,
                hits,
                gene_map,
                evalue_cutoff):
    """ Parse BLASTp alignment file into a dictionary where
        the keys are the queries and the values are all the
        reference sequences to which the query mapped with
        E-value cutoff score
    """
    # read blastp results
    with open(alignments_fp, 'U') as alignments_f:
        for line in alignments_f:
            line = line.split('\t')
            query = line[0]
            ref = line[1]
            e_value = line[10]
            # do not store alignments where the query
            # gene matched itself in the database
            if query == ref:
                continue
            if query not in hits:
                hits[query] = [ref]
            else:
                # check that the query mapped to a different species
                # since we only want the best homolog per species
                current_species = gene_map[ref].split('_')[0]
                add_alignment = True
                for gene in hits[query]:
                    species = gene_map[gene].split('_')[0]
                    if species == current_species:
                        add_alignment = False
                        break
                if add_alignment:
                  if float(e_value) <= float(evalue_cutoff):
                    hits[query].append(ref)


def launch_msa(fasta_in_fp,
               clustal_command_fp,
               gene_map,
               ref_db,
               hits,
               query,
               timeout,
               verbose=False,
               warnings=False):
    """ Create multiple sequence alignments for all gene othologs
    """
    with open(fasta_in_fp, 'w') as in_f:
        sorted(hits[query])
        for ref in hits[query]:
          in_f.write(">%s\n%s\n" % (gene_map[ref], ref_db[ref]))

    with open(clustal_command_fp, 'U') as clustal_command_f:
      clustalw_command = Command("clustalw")
      status, output, error = clustalw_command.run(timeout=timeout,
        stdin=clustal_command_f,
        close_fds=True)
      if status < 0:
        sys.stdout.write("status: %s\noutput: %s\terror: %s\t" % (status, output, error))


def compute_distances(phylip_command_fp,
                      phylip_fp,
                      verbose=False,
                      warnings=False):
    """ Compute distances between each pair of sequences in MSA
    """
    with open(phylip_command_fp, 'U') as phylip_command_f:
        proc = subprocess.Popen("protdist",
                     stdin=phylip_command_f,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     close_fds=True)
        proc.wait()
        stdout, stderr = proc.communicate()
        if stderr and warnings:
            print stderr


def normalize_distances(phylip_fp,
                        full_distance_matrix,
                        num_species,
                        full_distance_matrix_offset,
                        species_set_dict,
                        gene_bitvector_map):
    """ - Parse PHYLIP alignments and Z-score normalize each
        alignment entry
        - Compute bitvectors for all species sets
    """
    # assume a pairwise alignment exists for all species
    missing_species = [str(x) for x in range(0,num_species)]
    # scan through file and remove species that exist
    # from missing_species
    with open(phylip_fp, 'U') as phylip_f:
      for line in phylip_f:
        alignment_dist = line.strip().split()
        if len(alignment_dist) == 1:
          continue
        if not line.startswith(' '):
          species = line.split()[0].split('_')[0]
          missing_species.remove(species)

    # scan through file again, collecting alignment
    # distances
    orig_order_labels = []
    p = numpy.empty(shape=(num_species,num_species))
    p.fill(numpy.nan)
    idx_p = 0
    ind_a = 0
    with open(phylip_fp, 'U') as phylip_f:
        alignment_list = []
        # skip first line containing number of lines in
        # the file
        next(phylip_f)
        for line in phylip_f:
            alignment_dist = line.strip().split()
            if line.startswith(' '):
                alignment_list.extend(alignment_dist)
            else:
                # new species alignment pairs
                if alignment_list:
                    for i in range(0, len(missing_species)):
                      alignment_list.append(None)                
                    a = numpy.asarray(alignment_list[1:], dtype=float)
                    a[ind_a] = numpy.nan
                    ind_a += 1
                    mean = numpy.nanmean(a)
                    stdev = numpy.nanstd(a)
                    for x in numpy.nditer(a, op_flags=['readwrite']):
                      x[...] = (x-mean)/stdev
                    p[idx_p] = a
                    idx_p += 1
                    orig_order_labels.append(alignment_list[0])
                alignment_list = alignment_dist

    # add distance on final line
    for i in range(0, len(missing_species)):
      alignment_list.append(None)
    a = numpy.asarray(alignment_list[1:], dtype=float)
    a[ind_a] = numpy.nan
    mean = numpy.nanmean(a)
    stdev = numpy.nanstd(a)
    for x in numpy.nditer(a, op_flags=['readwrite']):
      x[...] = (x-mean)/stdev
    p[idx_p] = a
    orig_order_labels.append(alignment_list[0])

    # add the missing species names to the labels array
    bitvector_gene = 'I' * num_species
    for species in missing_species:
      orig_order_labels.append("%s_X" % species)
      # indicate missing gene for current species
      l = list(bitvector_gene)
      l[int(species)] = 'O'
      bitvector_gene = ''.join(l)

    # update species set counts
    if bitvector_gene not in species_set_dict:
      species_set_dict[bitvector_gene] = 1
    else:
      species_set_dict[bitvector_gene] += 1

    gene_bitvector_map[full_distance_matrix_offset] = bitvector_gene

    # sort the distance matrix based on species names (S1, S2, S3 ..)
    # in order to be consistent across all gene families
    sorted_order_labels = sorted(orig_order_labels)
    map_orig_sorted = {}
    for idx_orig, label in enumerate(orig_order_labels):
      idx_sorted = sorted_order_labels.index(label)
      map_orig_sorted[idx_orig] = idx_sorted

    # re-order rows and columns by ordered species name (0,1,2 ..)
    p2 = numpy.zeros(shape=(num_species,num_species))
    for idx_a, arr in enumerate(p):
      t = numpy.zeros(shape=num_species)
      for idx_b, el in enumerate(arr):
        t[map_orig_sorted[idx_b]] = el
      p2[map_orig_sorted[idx_a]] = t
    del p
 
    # add normalized distance matrix for current gene
    # to full distance matrix
    full_distance_matrix[full_distance_matrix_offset] = p2
                

def cluster_distances(species_set_dict,
                      species_set_size,
                      hamming_distance):
    """ Cluster gene families by species with detectable
        orthologs in exactly the same subset of the
        considered species
    """
    sorted_species_set = sorted(species_set_dict.items(),
      key=operator.itemgetter(1), reverse=True)

    # determine core clusters (initial species sets
    # with more than species_set_size genes)
    gene_clusters_dict = {}
    # if the largest species set contains less than
    # threshold (species_set_size) elements, set the
    # only core cluster to the largest species set
    if sorted_species_set[0][1] < species_set_size:
      gene_clusters_dict[sorted_species_set[0][0]] = []
    else:
      for bitvector in sorted_species_set:
        if bitvector[1] >= species_set_size:
          gene_clusters_dict[bitvector[0]] = []      

    # assign species sets with fewer than species_set_size
    # species to core clusters if the Hamming distance
    # between the two bitvectors is less than
    # hamming_distance
    species_set_assigned = []
    for cluster_core in gene_clusters_dict:
      for bitvector in sorted_species_set:
        if hamming(cluster_core, bitvector[0]) <= hamming_distance:
          gene_clusters_dict[cluster_core].append(bitvector[0])
          species_set_assigned.append(bitvector[0])

    # assign the remaining species sets to the
    # cluster with the closest core Hamming distance
    for bitvector in sorted_species_set:
      if bitvector[0] not in species_set_assigned:
        min_hamming_cluster = ""
        min_hamming_distance = sys.maxint
        # find cluster core with smallest Hamming distance
        # to species set
        for cluster_core in gene_clusters_dict:
          dist = hamming(cluster_core, bitvector[0])
          if dist < min_hamming_distance:
            min_hamming_distance = dist
            min_hamming_cluster = cluster_core
        gene_clusters_dict[min_hamming_cluster].append(bitvector[0])

    return gene_clusters_dict


def detect_outlier_genes(species_set,
        gene_bitvector_map,
        full_distance_matrix,
        stdev_offset,
        outlier_hgt,
        num_species,
        total_genes):
    """ Detect outlier genes
    """
    outlier_flag_matrix = numpy.zeros(shape=(total_genes,num_species,num_species), dtype=bool)
    distance_vector = numpy.zeros(total_genes)
    for i in range(num_species):
      for j in range(num_species):
        if i != j:
          for k in range(total_genes):
            distance_vector[k] = full_distance_matrix[k][i][j]
          mean = numpy.nanmean(distance_vector)
          stdev = numpy.nanstd(distance_vector)
          low_bound = mean - stdev_offset*stdev
          up_bound = mean + stdev_offset*stdev
          for k, distance in enumerate(distance_vector):
            if (distance != numpy.nan and ((distance < low_bound) or (distance > up_bound))):
              outlier_flag_matrix[k][i][j] = 1

    # traverse outlier_matrix by gene and count the number of
    # outlier distances by species
    outlier_count_matrix = numpy.zeros(shape=(total_genes,num_species), dtype=int)
    for i in range(total_genes):
      for j in range(num_species):
        for k in range(num_species):
          if outlier_flag_matrix[i][j][k]:
            outlier_count_matrix[i][k] += 1

    # if number of outlier distances exceeds threshold, label
    # gene as outlier
    outlier_genes = set()
    for i in range(total_genes):
      for j in range(num_species):
        if outlier_count_matrix[i][j] > num_species*outlier_hgt:
          outlier_genes.add(i)

    return outlier_genes


def output_full_matrix(matrix, num_species):
  """ Output distance matrix to stdout
  """
  for i in range(num_species):
    for j in range(num_species):
      # for gene number
      for k in range(len(matrix)):
        sys.stdout.write("%s\t" % matrix[k][i][j])
      sys.stdout.write("\n")



@click.command()
@click.argument('query-proteome-fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('target-proteomes-dir', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('working-dir', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=False,
                                file_okay=True))
@click.option('--min-num-homologs', type=int, required=False, default=3, show_default=True,
              help=('The mininum number of homologs (determined by BLAST search) '
                    'for each gene to test'))
@click.option('--e-value', type=float, required=False, default=10e-20, show_default=True,
              help=("The E-value cutoff to identify orthologous genes using BLASTP"))
@click.option('--threads', type=int, required=False, default=1, show_default=True,
              help=("Number of threads to use"))
@click.option('--stdev-offset', type=float, required=False, default=2.326, show_default=True,
              help=("The number of standard deviations a gene's normalized distance "
                    "is from the mean to identify it as an outlier for a species pair"))
@click.option('--outlier-hgt', type=float, default=0.5, show_default=True,
              required=False,
              help=('The fraction (value between (0,1]) of normalized pairwise '
                    'distances over all species-pair vectors belonging to the same '
                    'gene that are z-score standard deviations from the mean'))
@click.option('--species-set-size', type=int, required=False, default=30, show_default=True,
              help=('Threshold number of genes to consider a species set large'
                    '(a species set is a set of genes whose orthologs are detectable '
                    'in exactly the same subset of the considered species)'))
@click.option('--hamming-distance', type=int, required=False, default=2, show_default=True,
              help=('Distance between two binary vectors indicating the '
                    'species in which the corresponding ortholog gene appears'))
@click.option('--verbose', type=bool, required=False, default=False, show_default=True,
              help=("Run in verbose mode"))
@click.option('--warnings', type=bool, required=False, default=False, show_default=True,
              help=("Print program warnings"))
@click.option('--timeout', type=int, required=False, default=120, show_default=True,
              help=("Number of seconds to allow Clustalw to run per call"))
def _main(query_proteome_fp,
          target_proteomes_dir,
          working_dir,
          min_num_homologs,
          e_value,
          threads,
          stdev_offset,
          outlier_hgt,
          species_set_size,
          hamming_distance,
          verbose,
          warnings,
          timeout):
    """
    """
    if verbose:
        sys.stdout.write("Begin whole-genome HGT detection using the Distance method.\n\n")
        sys.stdout.write("Query genome: %s\n" % query_proteome_fp)

    gene_map, ref_db, num_species = preprocess_data(working_dir=working_dir,
        target_proteomes_dir=target_proteomes_dir,
        verbose=verbose)

    if verbose:
        sys.stdout.write("\nRunning BLASTp ..\n")
    hits = {}
    for _file in glob.glob(join(target_proteomes_dir, "*.fa")):
        # launch BLASTp
        alignments_fp = launch_blast(query_proteome_fp=query_proteome_fp,
            ref_fp=_file,
            working_dir=working_dir,
            e_value=e_value,
            threads=threads)

        # generate a dictionary of orthologous genes
        parse_blast(alignments_fp=alignments_fp,
                    hits=hits,
                    gene_map=gene_map,
                    evalue_cutoff=e_value)

    # keep only genes with >= min_num_homologs
    hits_min_num_homologs = {}
    max_homologs = 0
    for query in hits:
        len_hits = len(hits[query])
        if len_hits >= min_num_homologs:
            if query in hits_min_num_homologs:
                raise ValueError("Duplicate gene names found: %s" % query)
            hits_min_num_homologs[query] = hits[query]
            if len_hits > max_homologs:
                max_homologs = len_hits
    hits.clear()

    if verbose:
        sys.stdout.write("Total number of orthologous gene families with at least %s genes: %s\n" % (min_num_homologs, len(hits_min_num_homologs)))
    # generate command for CLUSTALW
    phy_msa_fp = join(working_dir, "msa.phy")
    dnd_msa_fp = join(working_dir, "msa.dnd")
    phylip_fp = join(working_dir, "msa.dis")
    # create fasta file for each gene family and run CLUSTALW
    fasta_in_fp = join(working_dir, "input.faa")
    clustal_command_fp = join(working_dir, "clustal_command.txt")
    with open(clustal_command_fp, 'w') as clustal_command_f:
        clustal_command_f.write('1\n%s\n2\n9\n1\n4\n\n1\n%s\n%s\nX\n\nX\n' % (fasta_in_fp, phy_msa_fp, dnd_msa_fp))
    phylip_command_fp = join(working_dir, "phylip_command.txt")
    with open(phylip_command_fp, 'w') as phylip_command_f:
        phylip_command_f.write('%s\nF\n%s\nR\nY\n' % (phy_msa_fp, phylip_fp))

    total_genes = len(hits_min_num_homologs)
    if verbose:
        sys.stdout.write("\nRunning CLUSTALW and PROTDIST ..\n")
    if max_homologs > num_species:
        raise ValueError("max_homologs > num_species: %s > %s " % (max_homologs, num_species))
    # distance matrix containing distances between all ortholog genes
    full_distance_matrix = numpy.zeros(shape=(total_genes,num_species,num_species), dtype=float)
    # dictionary to store all subsets of orthologs (keys) and
    # their number of occurrences (values) (maximum occurrences
    # is equal to the number of genes)
    species_set_dict = {}
    gene_bitvector_map = {}
    gene_id = {}
    for i, query in enumerate(hits_min_num_homologs):
      if verbose:
          print "Computing MSA and distances for gene %s .. (%s/%s)" % (query, i+1, total_genes)
      gene_id[i] = query
      # generate a multiple sequence alignment
      # for each orthologous gene family
      launch_msa(fasta_in_fp=fasta_in_fp,
          clustal_command_fp=clustal_command_fp,
          ref_db=ref_db,
          gene_map=gene_map,
          hits=hits_min_num_homologs,
          query=query,
          timeout=timeout,
          verbose=verbose,
          warnings=warnings)

      # compute distances between each pair of sequences in MSA
      compute_distances(phylip_command_fp=phylip_command_fp,
          phylip_fp=phylip_fp,
          verbose=verbose,
          warnings=warnings)

      # Z-score normalize distance matrix and add results
      # to full distance matrix (for all genes)
      normalize_distances(phylip_fp=phylip_fp,
          full_distance_matrix=full_distance_matrix,
          num_species=num_species,
          full_distance_matrix_offset=i,
          species_set_dict=species_set_dict,
          gene_bitvector_map=gene_bitvector_map)

    #output_full_matrix(full_distance_matrix, num_species)

    # cluster gene families by species
    gene_clusters_dict = cluster_distances(species_set_dict=species_set_dict,
      species_set_size=species_set_size,
      hamming_distance=hamming_distance)

    # detect outlier genes per core cluster of genes
    sys.stdout.write("Candidate HGT genes: \n")
    for core_cluster in gene_clusters_dict:
      outlier_genes = detect_outlier_genes(species_set=gene_clusters_dict[core_cluster],
        gene_bitvector_map=gene_bitvector_map,
        full_distance_matrix=full_distance_matrix,
        stdev_offset=stdev_offset,
        outlier_hgt=outlier_hgt,
        num_species=num_species,
        total_genes=total_genes)

      for gene in outlier_genes:
        sys.stdout.write("%s\n" % gene_id[gene])          

    #output_full_matrix(outlier_genes, num_species)





if __name__ == "__main__":
    _main()