#!/bin/bash

dirs="captain
cargo
ships"

files="fna
faa"

# generate checksum based on sequence
# ? strandedness ?
# ? adapting/adding information to headers ?

for dir in $dirs; do
  if [[ $dir == "ships" ]]; then
    out_dir=$dir/$file
    checksum_out="$out_dir"/"$file".checksums.txt
  else
    if [[ $dir == "cargo" ]]; then
      genes="nlr
      fre
      plp
      duf3723"
    elif [[ $dir == "captain" ]]; then
      gene="tyr"
    fi

    out_dir=$dir/$gene/$file
    checksum_out="$out_dir"/"$file".checksums.txt

  fi
  type=$(cat "$dir" | cut -d "/" -f 2)
  if [[ $type == "ships" ]]; then
    file=$(cat "$dir" | cut -d "/" -f 3)
  else
    gene=$(cat "$dir" | cut -d "/" -f 3)
    file=$(cat "$dir" | cut -d "/" -f 4)
  fi

  # return list of existing directories
  anno_dirs=$(find $out_dir -type d \( -name "manual" -o -name "starfish" \) 2>/dev/null)
  
  for anno_dir in $anno_dirs; do
    blast_dir=$anno_dir/blastdb
    anno_dir_name=$(dirname $anno_dir)
    # automate splitting of multifastas into separate fastas
    seqkit split -j 4 -f -i "$blast_dir"/*."$file" -O "$out_dir/$anno_dir_name"
  done

  # processing of checksums for all genes
  if [[ $file == "faa" ]]; then  
    find "$dir"/starbase -type f -name "*.faa" | seqkit sum -j 4 - | sed 's/seqkit.v0.1_DLS_k0_|seqkit.v0.1_PLS_k0_//g' > $checksum_out
  else
    find "$dir"/starbase -type f \( -name "*.fa" -o -name "*.fna" \) | seqkit sum -j 4 - | sed 's/seqkit.v0.1_DLS_k0_|seqkit.v0.1_PLS_k0_//g' > $checksum_out
  fi
done

# remove empty files
find ./ -name "*.checksums.txt" -type f -mindepth 3 -maxdepth 3 -empty -delete

CHECKSUMS=$(find ./ -name "*.checksums.txt" -type f -mindepth 3 -maxdepth 3)

# check for duplicates
for CHECKSUM in $CHECKSUMS; do
  DIR=$(basedir "$CHECKSUM")
  # number of entries
  wc -l "$DIR"/"$type".checksums.txt
  # number of unique entries
  cut -f 1 "$DIR"/"$type".checksums.txt | sort | uniq | wc -l
done
