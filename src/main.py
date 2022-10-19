#!/usr/bin/env python
import argparse
import pandas as pd
from os.path import exists
from Bio import SeqIO
import os
import subprocess
from tqdm import tqdm

def read_organism_table(org_table_filename):
    """
    Takes as argument the organism table filename with full path.
    Returns a dictionary indexed by 3-character organism code, mapping to a list
    of contig ids.
    :param str org_table_filename: full path to organism filename
    :return: dictionary indexed by 3-character organism code, mapping to a list
    of contig ids.
    :rtype: dic
    """
    # read dataframe
    df = pd.read_csv(org_table_filename, delimiter='\t')

    # keep only bacteria
    df['lineage'] = df['lineage'].str.lower()
    df = df[ df['lineage'].str.contains('bacteria') ]

    # take interesting things
    gb_ids = df['gb_ncbi_seq_id'].tolist()
    org_codes = df['org_code'].tolist()
    names = df['name'].tolist()

    # keep only non-empty ncbi things here
    ret_dic = {}
    for gb_id, org_code, name in list( zip(gb_ids, org_codes, names) ):
        if gb_id == '[]':
            continue
        else:
            gb_id = gb_id[1:-1].replace('\'', '').split(',')
            ret_dic[org_code] = gb_id

    return ret_dic

def read_organism_table_for_single_chr_bacteria(org_table_filename):
    """
    Takes as argument the organism table filename with full path.
    Returns a dictionary indexed by 3-character organism code, mapping to a list
    of contig ids. Only those bacteria with single chromosome are kept.
    :param str org_table_filename: full path to organism filename
    :return: dictionary indexed by 3-character organism code, mapping to a list
    of contig ids.
    :rtype: dic
    """
    # read dataframe
    df = pd.read_csv(org_table_filename, delimiter='\t')

    # keep only bacteria
    df['lineage'] = df['lineage'].str.lower()
    df = df[ df['lineage'].str.contains('bacteria') ]

    # take interesting things
    gb_ids = df['gb_ncbi_seq_id'].tolist()
    org_codes = df['org_code'].tolist()
    names = df['name'].tolist()

    # keep only non-empty ncbi things here
    ret_dic = {}
    for gb_id, org_code, name in list( zip(gb_ids, org_codes, names) ):
        if gb_id == '[]':
            continue
        else:
            gb_id = gb_id[1:-1].replace('\'', '').split(',')
            if len(gb_id) > 1:
                continue
            ret_dic[org_code] = gb_id

    return ret_dic

def read_gene_ids(kegg_gene_filename):
    """
    Given a kegg gene filename, returns list of all genes
    :param str kegg_gene_filename: filename with full path
    :return: all gene ids in that file
    :rtype: list
    """
    f = open(kegg_gene_filename, 'r')
    lines = f.readlines()[1:]
    f.close()
    return [ line.split(' ')[0].split('\t')[0] for line in lines ]

def read_gene_id_with_ko(kegg_gene_filename):
    """
    Given a kegg gene filename, returns list of all genes with KO ids available
    :param str kegg_gene_filename: filename with full path
    :return: all gene ids in that file
    :rtype: list
    """
    df = pd.read_csv(kegg_gene_filename, delimiter='\t')
    df = df[ ~pd.isnull(df['koid']) ]
    return(df['kegg_gene_id'].tolist())

def read_gene_id_with_kos_labeled(kegg_gene_filename):
    """
    Given a kegg gene filename, returns list of all genes with KO ids available
    :param str kegg_gene_filename: filename with full path
    :return: all gene ids in that file, along with their assigned KO ids
    :rtype: list of 4-tuples, tuple[0] is the gene id, tuple[1] is the KOid,  [2] is ntseq, [3] is aaseq
    """
    df = pd.read_csv(kegg_gene_filename, delimiter='\t')
    df = df[ ~pd.isnull(df['koid']) ]
    return(  list( zip( df['kegg_gene_id'].tolist(), df['koid'].tolist(), df['ntseq'].tolist(), df['aaseq'] ) )  )

def make_kegg_gene_file_name(org_code):
    """
    Using org_code (3-character code for organisms), returns kegg gene filename
    :param str org_code: 3-character code for organisms, internal to kegg database
    :return: kegg gene filename
    """
    return org_code + '_kegg_genes.txt'

def parse_args(): # pragma: no cover
    parser = argparse.ArgumentParser(description="This script works on a directory named 'extracted_genomes'. The script generates a list of all the genes in the mapping files, and lists the KOIDs of these genes.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = parser.parse_args()
    return args

def main(): # pragma: no cover
    args = parse_args()
    print('Test!')

    # get a list of all genomes in our directory
    genome_dir_names = [x[0] for x in os.walk(genome_path)][1:]
    for genome_dir in genome_dir_names[:10]:
        genome_name = genome_dir.split('/')[-1]
        print(genome_name)

    # for all of these genomes:
        # get list of all the gene ids
        # get gene_and_ko_list
        # index by gene_id
        # for only the genes present in ours:
            # get the ko_id
            # add <gene_id, ko_id> in a dataframe
    # write dataframe


if __name__ == '__main__': # pragma: no cover
    main()
