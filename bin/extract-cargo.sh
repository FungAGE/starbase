#!/bin/bash

# ! conda version was not compiled with BLAST support
BASEDIR=/home/adrian/Systematics
DIAMOND=$BASEDIR/bin/diamond-linux64/diamond

cargo_genes="duf3723 fre nlr plp"

# * query file cannot contain certain characters in the header
query=Starships/blastdb/concatenated.fa
query_file=$(basename -- "$query")
query_name=$(echo "$query_file" | cut -d . -f 1)
query_ext="${query_file##*.}"

for cargo in $cargo_genes
  do
    ref=Starships/cargo/$cargo/$cargo.mycoDB.faa
    ref_file=$(basename -- "$ref")
    ref_name=$(echo "$ref_file" | cut -d . -f 1)
    ref_ext="${ref_file##*.}"

    # TODO: add a check to prevent overwriting?
    # FIXME: fa != fna
    output=Starships/blastdb/$ref_name"."$query_ext

    if [[ $ref_ext == "fa" ]]; then
      exit
    elif [[ $ref_ext == "faa" ]]; then
      db_type="prot"
      if [[ $query_ext == "fa" ]]; then
        DIAMOND_CMD="$DIAMOND blastx "
      fi
      if [[ $query_ext == "faa" ]]; then
        DIAMOND_CMD="$DIAMOND blastp "
      fi
    fi
    
    ref_out=Starships/blastdb/$ref_file
    seqkit rmdup -n "$ref" > "$ref_out"

    makeblastdb -in "$ref_out" -input_type fasta -dbtype "$db_type" -parse_seqids
    $DIAMOND makedb --in "$ref_out" --db "$ref_out" --threads 4

    # blastx/blastp and fix fasta headers
    # only keep the top hits for each query
    
    # TODO: add gene family info to header?

    $DIAMOND_CMD --db "$ref_out" \
      -q $query \
      -f 6 qseqid qseq \
      -e 0.001 \
      --max-target-seqs 1 \
      -p 4 \
      -k 1 | sed 's/^/>/g' | tr "\t" "\n" > "$output"

    # makeblastdb -in "$output" -input_type fasta -dbtype "$db_type" -parse_seqids

    # TODO: add hmmer option?
done