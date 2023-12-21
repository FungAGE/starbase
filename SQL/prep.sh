#!/bin/bash

# script for pre-processing of sequences before storing in a database

# TODOs:
# bed -> gff, and splitting
# ? fill out databases for nucl sequences for all genes ?
# ? prot for ships is not needed, right ?

# generate checksum based on sequence
# ? strandedness ?
# ? adapting/adding information to headers ?

files="fna
faa"

genes="tyr
duf3723
fre
nlr
plp"

sources="mycodb
manual"

types="cargo
ships"

# gff files

# saves them in SQL dir...
python bin/bed2gff.py MTDB/mycodb.final.starships.bed

# ! this is restricted to mycodb set, until we get gffs/beds for manual set
fnas=$(find SQL/data/fna/ships/mycodb -name "*.fna")
for fna in $fnas
  do
    # make sure fasta headers match gff
    ship_code=$(basename "$fna" | cut -d _ -f 3)
    gff=$(find SQL/data/gff/mycodb/ -name "*$ship_name.gff")
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
          repo_dir=Starships/$type/$source
          mkdir -p "$sql_dir"/"$source"
          seqkit split -j 4 -f -i "$repo_dir"/*/*."$file" -O "$sql_dir"/"$source"
        done
      seqkit sum -j 4 "$sql_dir"/*/* | sed 's/seqkit.v0.1_DLS_k0_|seqkit.v0.1_PLS_k0_//g' > "$sql_dir"/"$file".checksums.txt

      # check for duplicates
      # wc -l $sql_dir/$type.checksums.txt
      # cut -f 1 $sql_dir/$type.checksums.txt | sort | uniq | wc -l
      # cut -f 1 $sql_dir/$type.checksums.txt | rev | cut -d _ -f 1 | rev | sort | uniq | wc -l

    elif [[ $file == "faa" ]]; then
      for gene in $genes
        do
          repo_dir=Starships/$type/$gene
          sql_dir="SQL/data/$file/$type/$gene"
          mkdir -p "$sql_dir/mycodb"
          
          # skip because we've already done this
          if [[ $gene != "tyr" ]]; then
            # automate splitting of multifastas into separate fastas
            seqkit split -j 4 -f -i "$repo_dir"/*."$file" -O "$sql_dir/mycodb"
          fi

          # processing of checksums for all genes
          seqkit sum -j 4 "$sql_dir"/*/* | sed 's/seqkit.v0.1_DLS_k0_|seqkit.v0.1_PLS_k0_//g' > "$sql_dir"/"$file".checksums.txt
        done
    fi
  done
done

# remove empty files
find SQL/ -name "*.checksums.txt" -type f -empty -delete