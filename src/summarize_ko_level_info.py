import argparse

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

    print(present_genes_filename, gene_abundance_filename)
