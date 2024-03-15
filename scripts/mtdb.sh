#!/bin/bash

cd MTDB/ || exit

# make sure to stash your credentials first...
mtdb manage -p

# read -s PASSWORD
# <type password>
# printf '{"username":"myname","password":"%s"}' $PASSWORD | mtdb update -u

# create the predb file
Rscript prep-mtdb.R

# ! must contain absolute file paths

# make sure you are not linked to a mtdb already
mtdb -u

# Create a mycotoolsdb from a predb file:
mtdb predb2mtdb ships_predb.tsv -s

# Create a mycotoolsdb referencing an alternative master database:
# predb2mtdb.py <PREDBFILE> <REFERENCEDB>

# add local genomes
mtdb update -a ./predb2mtdb_*/predb2mtdb.mtdb

# search
# obtain gene homologs using db2search
db2search -o cargo-search -c 8 -a hmmsearch -d ./ships/mycotoolsdb/mtdb/*.mtdb -qd ../SQL/data/hmm/cargo

# If there are more than 1000 homologs, truncate the resulting homologs around a gene of interest between 50 and 500 genes:
# fa2clus -f <FASTA>.fa --min_seq 50 --max_seq 500 -i <FOCAL_GENE>

# If there are less than 1000 homologs, construct a fastree, -f, if there are many samples, or remove that parameter for a robust IQ-TREE.
# fa2tree -i <FASTA>.fa 

# If the phylogeny can be improved, extract a highly supported node from 4 and reconstruct, or systematically truncate the dataset with 3 and reconstruct.
