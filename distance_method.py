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
                    reference_proteomes_dir,
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
    for _file in glob.glob(join(reference_proteomes_dir, "*.faa")):
        if verbose:
            sys.stdout.write("%s\t" % basename(_file))
        with open(_file, 'rb') as readfile:
            gene = 0
            for label, seq in parse_fasta(readfile):
                label = label.split()[0]
                # Blast parses reference ids to remove the |cl| from <|cl|ref_id
                # (update this for other parsing requirements)
                label = label.split('|')[1]
                ref_db[label] = seq
                sudo_label = "S%s_G%s" % (species, gene)
                gene_map[label] = sudo_label
                gene_map[sudo_label] = label
                gene+=1
        if verbose:
            sys.stdout.write("%s\n" % (gene+1))
        species+=1
    return gene_map, ref_db, species
            

def launch_blast(target_proteome_fp,
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
    out_file_fp = join(working_dir, "%s.blast" % basename(splitext(target_proteome_fp)[0]))
    blastp_command = ["blastp",
                      "-db", ref_fp,
                      "-query", target_proteome_fp,
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
                        full_distance_matrix):
    """ Parse PHYLIP alignments and Z-score normalize each
        alignment entry
    """
    with open(phylip_fp, 'U') as phylip_f:
        alignment_list = []
        alignment_dict = {}
        alignment_count = 0
        for line in phylip_f:
            sys.stdout.write(line)
            alignment_dist = line.strip().split()
            if len(alignment_dist) == 1:
                continue
            if line.startswith(' '):
                alignment_list.extend(alignment_dist)
            else:
                # list is empty (first alignment in file)
                if alignment_list:
                    # remove the 0.0 distance for the alignment of a gene vs. itself
                    del alignment_list[alignment_count]
                    a = numpy.asarray(alignment_list[1:], dtype=float)
                    mean = numpy.average(a)
                    stdev = numpy.std(a)
                    for x in numpy.nditer(a, op_flags=['readwrite']):
                        x[...] = (x-mean)/stdev
                    alignment_dict[alignment_list[0]] = a
                alignment_list = alignment_dist
                alignment_count+=1

        # sort by species name (ex. S1, S2, S3 .. Sn)
        print "alignment_dict = ", alignment_dict

                


@click.command()
@click.argument('target-proteome-fp', required=True,
                type=click.Path(resolve_path=True, readable=True, exists=True,
                                file_okay=True))
@click.argument('reference-proteomes-dir', required=True,
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
def _main(target_proteome_fp,
          reference_proteomes_dir,
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
        sys.stdout.write("Query genome: %s\n" % target_proteome_fp)

    gene_map, ref_db, num_species = preprocess_data(working_dir=working_dir,
        reference_proteomes_dir=reference_proteomes_dir,
        verbose=verbose)

    if verbose:
        sys.stdout.write("\nRunning BLASTp ..\n")
    hits = {}
    for _file in glob.glob(join(reference_proteomes_dir, "*.faa")):
        # launch BLASTp
        alignments_fp = launch_blast(target_proteome_fp=target_proteome_fp,
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
    full_distance_matrix = numpy.zeros(total_genes*num_species*(num_species-1))
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

        normalize_distances(phylip_fp=phylip_fp,
            full_distance_matrix=full_distance_matrix)
        i += 1




            



if __name__ == "__main__":
    _main()