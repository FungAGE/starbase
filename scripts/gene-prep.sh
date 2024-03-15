#!/bin/bash

# script for pre-processing of genes before storing in a database

# TODO: bed -> gff, and splitting
# TODO: ? fill out databases for nucl sequences for all genes ?
# TODO: ? prot for ships is not needed, right ?

files="fna
faa"

genes="tyr
duf3723
fre
nlr
plp"

sources="mycodb
manual"

types="cargo"

# create gff files here

# saves them in SQL dir...
python scripts/bed2gff.py MTDB/mycodb.final.starships.bed

# ! this is restricted to mycodb set, until we get gffs/beds for manual set
fnas=$(find SQL/data/fna/ships/mycodb -name "*.fna")
for fna in $fnas
  do
    # make sure fasta headers match gff
    ship_code=$(basename "$fna" | cut -d _ -f 3)
    gff=$(find metadata/ships/starfish/gff/starfish/ -name "*$ship_name.gff")
    header=$(grep "$ship_code" "$gff" | head -n 1 | cut -f 1)
    sed -i "1s/^>.*/>$header/" "$fna"
  done

# fna/faa files
for file in $files
  do
  for type in $types
    do
    if [[ $file == "fna" ]]; then
      sql_dir="SQL/data/$file/$type"
      for source in $sources
        do
          repo_dir=$type/$source
          mkdir -p "$sql_dir"/"$source"
          seqkit split -j 4 -f -i "$repo_dir"/*/*."$file" -O "$sql_dir"/"$source"
        done
    elif [[ $file == "faa" ]]; then
      for gene in $genes
        do
          repo_dir=$type/$gene
          sql_dir="SQL/data/$file/$type/$gene"
          mkdir -p "$sql_dir/mycodb"
          
          # skip because we've already done this
          if [[ $gene != "tyr" ]]; then
            # automate splitting of multifastas into separate fastas
            seqkit split -j 4 -f -i "$repo_dir"/*."$file" -O "$sql_dir/mycodb"
          fi

        done
    fi
  done
done