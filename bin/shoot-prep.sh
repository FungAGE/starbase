#!/bin/bash

# TODO: add input flags
# TODO: make options for prot/nucl

BASEDIR=/home/adrian/Systematics
threads=8
DIAMOND=~/diamond

QUERY_FASTA=$1
query_name=$(basename QUERY_FASTA .fa)

DB_NAME=$2

SHOOT=$BASEDIR/bin/SHOOT
SHOOT_RESULTS=$BASEDIR/Starship_Database/shoot/$DB_NAME

FASTA_DIR="$SHOOT_RESULTS/OrthoFinder/input/OrthoFinder"

cd "$SHOOT_RESULTS" || exit

#################
# run shoot
# find the folder of the most recent run
# ORTHO_OUT=$(find $FASTA_DIR/ -name "Results*" -type d -printf "%t - %p\n" | sort -n | tail -1 | awk '{ print NF }')
ORTHO_OUT=$(find "$FASTA_DIR"/ -name "Results*" -type d)

if [[ $type == "protein" ]];
  then
    # ! can't find where this file is created
    # * possible that this file is created by cat'ing orthogroup sequences and renaming headers? 
    DIAMOND_INPUT=$ORTHO_OUT/profile_sequences.all.fa
    elif [[ $type == "nucl" ]];
      then
        # DIAMOND_INPUT=...
        exit
      else
        exit
  fi

# create shoot database from orthofinder db
# options are: full, profiles, mmseqs
python $SHOOT/shoot/create_shoot_db.py "$ORTHO_OUT" full

# Resolve polytomies (only necessary if using EPA-ng)
# python $SHOOT/shoot/bifurcating_trees.py $ORTHO_OUT

# 2. get protein sequence of fasta query and make a tree with all seqs in database
# cannot find a way to grab top query seq without having to run diamond again later

# filter the database by list of accessions `--seqidlist`
# this is only supported with blast databases
# generate blast database
mkdir "$SHOOT_RESULTS"/blastdb

cp "$DIAMOND_INPUT" "$SHOOT_RESULTS"/blastdb/shoot_database
if [[ $type == "protein" ]];
  then
    makeblastdb -in "$SHOOT_RESULTS"/blastdb/shoot_database -input_type fasta -dbtype prot -parse_seqids
    elif [[ $type == "nucl" ]];
      then
        # makeblastdb -in "$SHOOT_RESULTS"/blastdb/shoot_database -input_type fasta -dbtype nucl -parse_seqids
        exit
      else
        exit
  fi

# so have to make a separate db and prep it with `diamond prepdb`
$DIAMOND prepdb --db "$SHOOT_RESULTS"/blastdb/shoot_database --threads "$threads"