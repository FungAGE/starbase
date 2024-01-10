#!/bin/bash

# ! conda version was not compiled with BLAST support

diamond makedb --in Starships/blastdbYRsuperfamRefs.faa --db Starships/blastdbYRsuperfamRefs.faa --threads 4

# blast and fix fasta headers
# TODO: add tyr family to header?
diamond blastx --db Starships/blastdbYRsuperfamRefs.faa \
  -q Starships/blastdbconcatenated.fa \
  -f 6 qseqid qseq \
  -e 0.001 \
  -p 4 \
  -k 1 | sed 's/^/>/' | tr "\t" "\n" > Starships/blastdbYRsuperfamRefs.fa

makeblastdb -in Starships/blastdbYRsuperfamRefs.fa -input_type fasta -dbtype nucl -parse_seqids

# TODO: add hmmer option?