#!/bin/bash

threads=8

# ! conda version was not compiled with BLAST support
BASEDIR=/home/adrian/Systematics
PROJECT_DIR=$BASEDIR/Starship_Database/starbase/Starships
DIAMOND=$BASEDIR/bin/diamond-linux64/diamond

genes="duf3723 fre nlr plp tyr"
ships=$(find $PROJECT_DIR/ships/fna/manual $PROJECT_DIR/ships/fna/starfish -type f \( -name "*.fa" -o -name "*.fna" \))

# make sure outdirs exist
# find $BASEDIR/Starship_Database/Starships/$gene_dir/*/ -name "faa" -type d -exec mkdir -p {}/blast \;

for gene in $genes; do
  if [[ $gene == "tyr" ]]; then
    ref=YRsuperfamRefs.faa
    gene_dir="captain"
  else
    ref=$gene.mycoDB.faa
    gene_dir="cargo"
  fi
  BLAST_DIR=$PROJECT_DIR/$gene_dir/$gene/faa/blastdb
  ref_file="$BLAST_DIR"/"$ref"
  ref_name=$(basename -- "$ref_file")
  ref_ext="${ref_name##*.}"

  if [[ $ref_ext == "fa" ]]; then
    exit
  elif [[ $ref_ext == "faa" ]]; then
    db_type="prot"
  fi
  ref_out=$BLAST_DIR/${ref_name%."$ref_ext"}.dd.$ref_ext
  # problems with duplicates and headers being too long...
  seqkit replace -p "^(\w+.*)" -r {nr} "$ref_file" | seqkit rmdup -n - > "$ref_out"

  # make BLAST database and format it for diamond
  makeblastdb -in "$ref_out" -input_type fasta -dbtype "$db_type" -parse_seqids
  $DIAMOND makedb --in "$ref_out" --db "$ref_out" --threads $threads

  for query in $ships; do
    # * query file cannot contain certain characters in the header
    query_file=$(basename -- "$query")
    query_name=$(echo "$query_file" | cut -d . -f 1)
    query_ext="${query_file##*.}"

    # remove gap characters
    sed -i '2,$s/-//g' "$query"

    if [[ $query_ext == "fa" ]]; then
        DIAMOND_CMD="$DIAMOND blastx "
    elif [[ $query_ext == "faa" ]]; then
        DIAMOND_CMD="$DIAMOND blastp "
    fi

    # TODO: add a check to prevent overwriting?
    # FIXME: fa != fna
    fas_output=$PROJECT_DIR/$gene_dir/$gene/fna/blast/$query_name".$gene."$query_ext
    faa_output=${fas_output%."$query_ext"}.faa
    tab_output=${fas_output%."$query_ext"}.txt

    # TODO: if there is only one sequence in DB with a matching header, pull that sequence instead of running diamond
    query_header=$(grep "^>" "$query" | sed 's/>//g' | cut -d "_" -f 1)
    seqkit grep -r -p "$query_header.*" "$ref_file" > "$fas_output"
    n_seq=$(grep -c "^>" "$fas_output")

    if [[ $n_seq != 1 ]]; then
      rm -f "$fas_output"

      # blastx/blastp and fix fasta headers
      # only keep the top hits for each query
      $DIAMOND_CMD --db "$ref_out" \
        -q "$query" \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qseq_translated \
        -e 0.001 \
        --max-target-seqs 1 \
        -p $threads > "$tab_output"
      
      # TODO: add gene family info to header?
      # TODO: add in a check for low quality hits
      # nucl fasta
      cut -f1,13 "$tab_output" | sed 's/^/>/g' | tr "\t" "\n" > "$fas_output"
      # protein fasta
      cut -f1,14 "$tab_output" | sed 's/^/>/g' | tr "\t" "\n" > "$faa_output"

    else
      echo "$gene gene already found in $query_name"
    fi
  done
  # TODO: add hmmer option?
done