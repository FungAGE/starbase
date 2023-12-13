#!/bin/bash
# use shoot.bio to place user sequences in a phylogeny

threads=8

QUERY_FASTA=$1
DB_NAME=$2

if [[ $DB_NAME == "tyr" ]]; then
  DB=shoot/YRsuperfamRefs/blastdb/YRsuperfamRefs
  ORTHO_OUT=shoot/YRsuperfamRefs/OrthoFinder/input/OrthoFinder/Results_Aug03
else
  exit
fi

# ? .sh.ogs.txt.gz is required for orthofinder?
diamond blastp --db "$DB" \
  -q "$QUERY_FASTA" \
  -o "$QUERY_FASTA".sh.ogs.txt.gz \
  --quiet \
  -e 0.001 \
  --compress 1 \
  -p "$threads"

# note that the query protein fasta file cannot contain non-alpha characters, i.e. `*`
sed -i '1!s/[^[:alpha:]]//g' "${QUERY_FASTA}"

# seem to need absolute path to folder
SHOOT=/home/adrian/Systematics/Starship_Database/starbase/bin/SHOOT/shoot/
# now run shoot to create a phylogeny of just the orthologous sequences
# must specify -t and -n otherwise it will fail
python "$SHOOT" --profiles_all -f -t iqtree -n "$threads" "${QUERY_FASTA}" "$ORTHO_OUT"
