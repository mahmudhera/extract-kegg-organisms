There are 4496 genomes in total. These genomes are:

1. All bacterial genomes
1. There is only a single chromosome
1. Genome dir name contains 3/4 digit org code
1. Have fasta files in genome dir named as: org_code.fasta
1. All genome dirs have a mapping file
1. The mapping file contains gene names, start and end positions, strand info etc.

As of June 2, 2023: these genomes are here in the GPU machine:
/scratch/mbr5797/genomes_extracted_from_kegg


Nore: Some gene info is not in the mapping files.
These are genes that originate from two different locations in the genome.
Out of 7964903 genes, count of missing spliced genes: only 5772

## Install
```
conda create -y --name extract_kegg
conda install -y --name extract_kegg -c conda-forge -c bioconda --file requirements.txt
conda activate extract_kegg
```

## Run
`python src/main.py`

## Present genomes
Listed in list_of_genomes
