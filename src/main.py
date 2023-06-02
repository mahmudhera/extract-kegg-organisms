#!/usr/bin/env python
import argparse
import pandas as pd
from os.path import exists
from Bio import SeqIO
import os
import subprocess
from tqdm import tqdm
import sys


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

    genome_path = '/scratch/mbr5797/genomes_extracted_from_kegg'
    directory_with_kegg_gene_files = '/scratch/shared_data/archive/KEGG_data/organisms/kegg_gene_info'

    genome_names = []
    present_genes_and_kos = []
    genome_dir_names = [x[0] for x in os.walk(genome_path)][1:]
    for genome_dir in tqdm(genome_dir_names):
        genome_name = genome_dir.split('/')[-1]
        mapping_filename = genome_name + '_mapping.csv'

        genome_names.append(genome_name)

        df = pd.read_csv(genome_dir+'/'+mapping_filename, delimiter=',')
        genes_present_in_mapping_file = df['gene_name'].tolist()

        org_code = genome_name
        gene_filename = make_kegg_gene_file_name(org_code)
        gene_file_with_path = directory_with_kegg_gene_files + '/' + gene_filename
        gene_and_ko_list = read_gene_id_with_kos_labeled(gene_file_with_path)

        gene_id_to_ko_id = {}
        for (gene_id, ko_id, nt_seq, aa_seq) in gene_and_ko_list:
            gene_id_to_ko_id [gene_id] = ko_id

        for present_gene in genes_present_in_mapping_file:
            present_genes_and_kos.append( (present_gene, gene_id_to_ko_id [present_gene]) )

    df = pd.DataFrame(present_genes_and_kos, columns=['gene_id', 'ko_id'])
    df.to_csv('present_genes_and_koids.csv')

    with open('list_of_genomes', 'w') as sys.stdout:
        for genome_name in genome_names:
            print(genome_name)


if __name__ == '__main__': # pragma: no cover
    main()
