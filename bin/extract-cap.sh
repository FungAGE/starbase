#!/bin/bash

BASEDIR=/home/adrian/Systematics
# ! conda version was not compiled with BLAST support
DIAMOND=$BASEDIR/bin/diamond-linux64/diamond

# problem with headers being too long...
# seqkit replace -p "^(\w+.*)" -r {nr} $BASEDIR/Starship_RNAseq/expressionMetadata/Afumigatus/aspfum524/from-emile/aspfum524.extend-final.starships.clean.proteins.faa > $BASEDIR/Starship_Database/blastdb/starship_aa.faa

# makeblastdb -in blastdb/YRsuperfamRefs.faa -input_type fasta -dbtype prot -parse_seqids

$DIAMOND makedb --in blastdb/YRsuperfamRefs.faa --db blastdb/YRsuperfamRefs.faa --threads 4

# blast and fix fasta headers
# TODO: add tyr family to header?
$DIAMOND blastx --db blastdb/YRsuperfamRefs.faa \
  -q blastdb/concatenated.fa \
  -f 6 qseqid qseq \
  -e 0.001 \
  -p 4 \
  -k 1 | sed 's/^/>/' | tr "\t" "\n" > blastdb/YRsuperfamRefs.fa

makeblastdb -in blastdb/YRsuperfamRefs.fa -input_type fasta -dbtype nucl -parse_seqids

# TODO: add hmmer option?