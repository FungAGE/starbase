#!/bin/bash

BASEDIR=/home/adrian/Systematics
threads=8

QUERY_FASTA=$1
query_name=$3
DB_NAME=$3

if [[ $DB_NAME == "tyr" ]]; then
  SHOOT_RESULTS=$BASEDIR/Starship_Database/shoot/YRsuperfamRefs/blastdb/
  cd "$SHOOT_RESULTS" || exit
  ORTHO_OUT=$BASEDIR/Starship_Database/shoot/YRsuperfamRefs/OrthoFinder/input/OrthoFinder/Results_Aug03
else
  exit
fi

SHOOT=$BASEDIR/bin/SHOOT

$DIAMOND blastx --db "$DB" \
  -q "$QUERY_FASTA" \
  -f 6 qseq_translated \
  -e 0.001 \
  -p "$threads" \
  --seqidlist cap_tyr.seq \
  --skip-missing-seqids \
  -k 1 | sed "1i >$query_name" > "${QUERY_FASTA%.fa}".diamond.fa

# error in running shoot on custom database
# try running diamond manually:
# https://github.com/davidemms/$SHOOT/issues/3#issuecomment-1210659139

$DIAMOND blastp --db "$DIAMOND_INPUT" \
  -q "${QUERY_FASTA%.fa}".diamond.fa \
  -o "${QUERY_FASTA%.fa}".diamond.sh.ogs.txt.gz \
  --quiet \
  -e 0.001 \
  --compress 1 \
  -p "$threads"

# 3. use shoot.bio to place user sequences in a phylogeny
# now run shoot to create a phylogeny of just the orthologous sequences
# note that the query protein fasta file cannot contain non-alpha characters, i.e. `*`
sed -i '1!s/[^[:alpha:]]//g' "${QUERY_FASTA%.fa}".diamond.fa

# must specify -t and -n otherwise it will fail
python "$SHOOT"/shoot --profiles_all -f -t iqtree -n "$threads" "${QUERY_FASTA%.fa}".diamond.fa "$ORTHO_OUT"
