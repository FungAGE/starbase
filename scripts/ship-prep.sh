#!/bin/bash

# This script takes the main ship fastas creates blastdbs and split the fasta into individual files

STARSHIPS=$(find ships/* -type d -maxdepth 1 -mindepth 1)
blastdb=blastdb/concatenated.fa

# TODO: once tax info is available for mycodb set...
for SHIPDIR in $STARSHIPS; do
  ANNO_VERSION=$(basename "$SHIPDIR")
  
  if [[ $SHIPDIR == *"manual"* ]]; then
    SHIPFASTA=$SHIPDIR/Starships.fa
    SHIPSQL=SQL/data/fna/ships/"$ANNO_VERSION"
  
    # ? should all existing fastas be removed? to avoid conflicts
    rm -f SQL/data/fna/ships/"$ANNO_VERSION"/*.{fa,fna}
    
    # split fasta up, output is in a new directory for use in SQL
    seqkit split -f -i "$SHIPFASTA" -O "$SHIPSQL"

    # fill in taxonomic information
    tail -n +2 "$ANNO_VERSION"/Starships.csv | tr "," "\t" | taxonkit name2taxid -i 2 -r > "$SHIPDIR"/Starships.withtaxaids.tsv
    
    # FIXME: manually added missing taxids here

    taxonkit lineage -i 12 -R "$SHIPDIR"/Starships.withtaxaids.tsv > "$SHIPDIR"/Starships.withtaxa.tsv

    # TODO: replace spaces and other characters (i.e. parentheses) in species names
    # creates $SHIPDIR/Starships.fulltaxa.csv
    Rscript scripts/format_taxonomy.R

    # problem with headers being too long...
    seqkit replace -p "^(\w+.*)" -r {nr} "$SHIPFASTA" > "$SHIPFASTA".new

    # remove characters that mess up fasta headers
    sed -i -E 's/[^\x00-\x7F]+/ /g' "$SHIPFASTA".new

  else
    SHIPDIR=$SHIPDIR/output
    SHIPFASTA=$SHIPDIR/mycodb.final.starships.fna
    # fix spaces and "|" in fasta headers
    sed 's/ /_/g' "$SHIPFASTA" | seqkit replace -p "^(\w+.*)" -r {nr} "$SHIPFASTA" | sed -E 's/[^\x00-\x7F]+/ /g' > "$SHIPFASTA".new
  fi
done

find ships/ -type f -name "*.new" -exec cat {} + > "$blastdb"
makeblastdb -in "$blastdb" -dbtype nucl -input_type fasta -parse_seqids

# clean up
find ships/ -type f -name "*.new" -delete
