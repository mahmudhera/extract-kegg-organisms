import argparse
import pandas as pd

def parse_args(): # pragma: no cover
    parser = argparse.ArgumentParser(description="This script will take a list of <gene_id, ko_id> mapping, and also take the gene abundance information in a simulation. Then, this script will summarize the abundane information from gene level to KO level.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--present_genes', type=str, help='file containing list of present genes and KO ids. Generated by main.py of the repo: https://github.com/mahmudhera/extract-kegg-organisms')
    parser.add_argument('--gene_abundance', type=str, help='file containing gene abundance information of a simulation. Generated by the script find_genes_in_sim.py of the repo: https://github.com/mahmudhera/KEGG_sketching_annotation')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    present_genes_filename = args.present_genes
    gene_abundance_filename = args.gene_abundance

    gene_abundance_df = pd.read_csv(gene_abundance_filename)
    gene_names = gene_abundance_df['gene_name'].tolist()
    gene_lengths = gene_abundance_df['gene_length'].tolist()
    nt_ovelaps = gene_abundance_df['nucleotide_overlap'].tolist()
    list_reads_mapped = gene_abundance_df['reads_mapped'].tolist()

    total_num_reads = sum(list_reads_mapped)
    total_nucleotides_covered = sum(nt_ovelaps)

    test = 0
    test2 = 0

    for gene_name, gene_length, nt_overlap, num_reads_mapped in list( zip(gene_names, gene_lengths, nt_ovelaps, list_reads_mapped) ):
        abundance_estimate_1 = 1.0 * num_reads_mapped / total_num_reads
        abundance_estimate_2 = 1.0 * nt_overlap / total_nucleotides_covered

        test += abundance_estimate_1
        test2 += abundance_estimate_2

    print(test)
    print(test2)
