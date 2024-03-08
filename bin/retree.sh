#!/bin/bash
# phylogenetic placement

threads=8

QUERY=$1
type=$2

if [[ $type == "prot" ]]; then
  REF="Starships/captain/tyr/alignments/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit"
  TREE="Starships/captain/tyr/tree/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile"
  HMM="Starships/SQL/data/hmm/captain/YRsuperfams.p1-512.hmm"
else
  REF="Starships/captain/tyr/YRsuperfamRefs.fa"
  TREE="Starships/captain/tyr/tree/YRsuperfamRefs.trimal.fa.iqtree"
fi

NAME=$(basename $REF)
QRY_HEADER=$(head -n 1 "$QUERY" | sed 's/>//g')

# rename and align
seqkit replace -p "$QRY_HEADER" -r "QUERY" "$QUERY"  > tmp/"$QRY_HEADER".fa
mafft --adjustdirection --auto --thread "$threads" --add tmp/"$QRY_HEADER".fa $REF > tmp/"$QRY_HEADER".mafft

# trim
# trimal -in tmp/$NAME.mafft -out tmp/$NAME.trimal

# make a tree
# iqtree2 -s tmp/$QRY_HEADER.mafft -T AUTO

seqkit grep -n -p "QUERY" tmp/"$QRY_HEADER".mafft > tmp/"$QRY_HEADER".mafft.query
seqkit grep -v -n -p "QUERY" tmp/"$QRY_HEADER".mafft > tmp/"$QRY_HEADER".mafft.ref

raxml-ng --redo --evaluate --msa $REF --tree $TREE --prefix tmp/info --model LG+G8+F --threads $threads
epa-ng --redo --ref-msa tmp/"$QRY_HEADER".mafft.ref --tree $TREE --query tmp/"$QRY_HEADER".mafft.query -T $threads --model tmp/info.raxml.bestModel

# pplacer using HMMER
# ? reference alignment in Stockholm format?
# Then we can use it to make a combined alignment with the reference sequences and the reads:

hmmalign -o tmp/combo.sto --mapali $REF "$HMM" tmp/"$QRY_HEADER".fa

# Now we can run pplacer:
pplacer -t tmp/combo.tre -s tmp/info.raxml.bestModel
