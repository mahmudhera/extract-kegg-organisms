There are 114 genomes in total.
In every genome, there is the fasta sequence, and the mapping file.
The mapping file contains information of a gene -- genes that are present within that organism.
The start location, end location, nucleotide sequence, and amino acid sequence -- all are these.

All genomes have a single chromosome.

Some gene info is not in the mapping files.
These are genes that originate from two different locations in the genome.
Out of these 114 genomes, there are only 137 such genes.

## Install
```
conda create -y --name extract_kegg
conda install -y --name extract_kegg -c conda-forge -c bioconda --file requirements.txt
conda activate extract_kegg
```
