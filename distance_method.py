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
        iv.   Compute pairwise distance matrix using PHYLIP's protdist function;
        v.    Add distance matrix for all pairwise distances into global distance
              matrix storing results for all genes

    2. Cluster gene families by species,
        vi.   Cluster gene families according to Hamming distance;
        vii.  For each gene family and species, z-score normalize the set of
              pairwise distances between the gene in that species and all other
              species;
        viii. Cluster pairwise distances to each gene cluster;
        ix.   Run outlier detection algorithm on each cluster (paragraph 2
              of section 'Detecting Outlier Genes' in original paper)
"""

import sys
import click
import glob
from os.path import join, splitext, basename
from subprocess import Popen, PIPE
from string import maketrans
import numpy

from skbio.parse.sequences import parse_fasta


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
    for _file in glob.glob(join(target_proteomes_dir, "*.faa")):
        if verbose:
            sys.stdout.write("%s. %s\t" % (species+1, basename(_file)))
        with open(_file, 'rb') as readfile:
            gene = 0
            for label, seq in parse_fasta(readfile):
                label = label.split()[0]
                # Blast parses reference ids to remove the |cl| from <|cl|ref_id
                # (update this for other parsing requirements)
                label = label.split('|')[1]
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
    proc = Popen(makeblastdb_command,
                 stdout=PIPE,
                 stderr=PIPE,
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
    proc = Popen(blastp_command,
                 stdout=PIPE,
                 stderr=PIPE,
                 close_fds=True)
    proc.wait()
    stdout, stderr = proc.communicate()
    if stderr:
        print stderr

    return out_file_fp


def parse_blast(alignments_fp,
                hits,
                gene_map):
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
                    hits[query].append(ref)


def launch_msa(fasta_in_fp,
               clustal_command_fp,
               gene_map,
               ref_db,
               hits,
               query,
               verbose=False,
               warnings=False):
    """ Create multiple sequence alignments for all gene othologs
    """
    with open(fasta_in_fp, 'w') as in_f:
        sorted(hits[query])
        for ref in hits[query]:
            in_f.write(">%s\n%s\n" % (gene_map[ref], ref_db[ref]))
    # run CLUSTALW
    with open(clustal_command_fp, 'U') as clustal_command_f:
        proc = Popen("clustalw",
                     stdin=clustal_command_f,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)
        proc.wait()
        stdout, stderr = proc.communicate()
        if stderr and warnings:
            print stderr


def compute_distances(phylip_command_fp,
                      phylip_fp,
                      verbose=False,
                      warnings=False):
    """ Compute distances between each pair of sequences in MSA
    """
    with open(phylip_command_fp, 'U') as phylip_command_f:
        proc = Popen("protdist",
                     stdin=phylip_command_f,
                     stdout=PIPE,
                     stderr=PIPE,
                     close_fds=True)
        proc.wait()
        stdout, stderr = proc.communicate()
        if stderr and warnings:
            print stderr


def normalize_distances(phylip_fp,
                        full_distance_matrix,
                        num_species,
                        full_distance_matrix_offset):
    """ Parse PHYLIP alignments and Z-score normalize each
        alignment entry.
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
        for line in phylip_f:
            alignment_dist = line.strip().split()
            if len(alignment_dist) == 1:
                continue
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
    for species in missing_species:
      orig_order_labels.append("%s_X" % species)

    # sort the distance matrix based on species names (S1, S2, S3 ..)
    # in order to be consistent across all gene families
    sorted_order_labels = sorted(orig_order_labels)
    map_orig_sorted = {}
    for idx_orig, label in enumerate(orig_order_labels):
      idx_sorted = sorted_order_labels.index(label)
      map_orig_sorted[idx_orig] = idx_sorted

    # re-order rows by ordered species name (0,1,2 ..)
    p2 = numpy.zeros(shape=(num_species,num_species))
    for idx, arr in enumerate(p):
      p2[map_orig_sorted[idx]] = arr
    del p
 
    # re-order columns by ordered species name (0,1,2 ..)
    p3 = numpy.zeros(shape=(num_species,num_species))
    for idx_a, arr in enumerate(p2):
      t = numpy.zeros(shape=num_species)
      for idx_b, el in enumerate(arr):
        t[map_orig_sorted[idx_b]] = el
      p3[idx_a] = t
    del p2

    # add normalized distance matrix for current gene
    # to full distance matrix
    full_distance_matrix[full_distance_matrix_offset] = p3
                


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
@click.option('--z-score', type=float, required=False, default=1.5, show_default=True,
              help=("The number of standard deviations a gene's distance is from the mean "
                    "to identify it as an outlier for a species pair"))
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
def _main(query_proteome_fp,
          target_proteomes_dir,
          working_dir,
          min_num_homologs,
          e_value,
          threads,
          z_score,
          outlier_hgt,
          species_set_size,
          hamming_distance,
          verbose,
          warnings):
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
    for _file in glob.glob(join(target_proteomes_dir, "*.faa")):
        # launch BLASTp
        alignments_fp = launch_blast(query_proteome_fp=query_proteome_fp,
            ref_fp=_file,
            working_dir=working_dir,
            e_value=e_value,
            threads=threads)

        # generate a dictionary of orthologous genes
        parse_blast(alignments_fp=alignments_fp,
                    hits=hits,
                    gene_map=gene_map)

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
    i = 0
    if verbose:
        sys.stdout.write("\nRunning CLUSTALW and PROTDIST ..\n")
    if max_homologs > num_species:
        raise ValueError("max_homologs > num_species: %s > %s " % (max_homologs, num_species))
    # distance matrix containing distances between all ortholog genes
    full_distance_matrix = numpy.zeros(shape=(total_genes,num_species,num_species), dtype=float)
    for query in hits_min_num_homologs:
        if verbose:
            print "Computing MSA and distances for gene %s .. (%s/%s)" % (query, i, total_genes)
        # generate a multiple sequence alignment
        # for each orthologous gene family
        launch_msa(fasta_in_fp=fasta_in_fp,
            clustal_command_fp=clustal_command_fp,
            ref_db=ref_db,
            gene_map=gene_map,
            hits=hits_min_num_homologs,
            query=query,
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
            full_distance_matrix_offset=i)

        i += 1




            



if __name__ == "__main__":
    _main()