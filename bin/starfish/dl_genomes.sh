#!/bin/bash

# Initialize variables
API="12c0c894e9b2d5bb11e5a5236b3b0d2a9d09"

species_flag=0
accession_flag=0
species=""
accession=""

QUERY="$1"
GENOME_PATH=$2
run_name=$3

# Function to display usage
usage() {
  echo "Usage: $0 [--species <species_name>] [--accession <accession_number(s)>] [<query input>] [<output dir>]"
  exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --species)
      if [ $accession_flag -eq 1 ]; then
        echo "Error: Cannot use --species and --accession together."
        usage
      fi
      species_flag=1
      # must be quoted for genus species query
      QUERY='"'$QUERY"'"
      shift 2
      ;;
    --accession)
      if [ $species_flag -eq 1 ]; then
        echo "Error: Cannot use --species and --accession together."
        usage
      fi
      accession_flag=1
      shift 2
      ;;
    *)
      echo "Error: Invalid argument '$1'."
      usage
      ;;
  esac
done

# Check if one of the flags is set and the corresponding input is provided
if [ $species_flag -eq 1 ]; then
  if [ -z "$QUERY" ]; then
    echo "Error: Missing species name after --species flag."
    usage
  fi
elif [ $accession_flag -eq 1 ]; then
  if [ -z "$QUERY" ]; then
    echo "Error: Missing accessions after --accession flag."
    usage
  fi
else
  echo "Error: Either --species or --accession flag must be provided."
  usage
fi

# redundant, but just to be safe
ROOT_DIR=$(pwd)

GENOME_DIR=$ROOT_DIR/Genomes/$GENOME_PATH
mkdir -p "$GENOME_DIR"

cd "$GENOME_DIR" || exit

# download assemblies and annotations from NCBI
# ! currently this is a very lax way to download genomes
# perhaps a more concious way would be to first use `--dehydrated` and then check to see how many genomes exist, what quality, which ones have annotations, then decide which ones and how many to dl

# * for `download genome accession`:
# * --include string(,string)   Specify the data files to include (comma-separated).
# * --inputfile string          Read a list of NCBI Assembly or BioProject accessions from a file to use as input

INCLUDES="genome,gff3,protein,cds"

# TODO: how to handle "invalid" accessions?

if [ $species_flag -eq 1 ]; then
  $ROOT_DIR/bin/datasets download genome taxon --api-key $API "$QUERY" --include "$INCLUDES" --exclude-atypical --tax-exact-match --assembly-version 'latest' --assembly-source 'RefSeq' --filename $run_name.zip
elif [ $accession_flag -eq 1 ]; then
  $ROOT_DIR/bin/datasets download genome accession --api-key $API "$QUERY" --include "$INCLUDES" --exclude-atypical --assembly-version 'latest' --assembly-source 'RefSeq' --filename $run_name.zip
fi

# unzip
unzip $run_name.zip

n_files=$(find "$GENOME_DIR"/$run_name/ncbi_dataset/ -type f \( -name "*.fna" -o -name "*.gff*" -o -name "*.faa" \) -size +0c | wc -l)

if [ "$n_files" -gt 0 ]; then
  # rename files (to match fna file)
  find "$GENOME_DIR"/$run_name/ncbi_dataset/* -type d -exec mv -n -- {}/cds*.fna {}.cds.fna \; -empty -delete
  find "$GENOME_DIR"/$run_name/ncbi_dataset/* -type d -exec mv -n -- {}/*.gff* {}.gff \; -empty -delete
  find "$GENOME_DIR"/$run_name/ncbi_dataset/* -type d -exec mv -n -- {}/*.faa {}.faa \; -empty -delete

  # subroutine for checking if gff's, or protein fastas exist
  # TODO: handle empty results

  # Create a TSV file to store the results
  output_file=$GENOME_DIR/$run_name/"file_list.tsv"

  # Clear the output file if it exists
  > "$output_file"

  # Loop through subdirectories and search for files
  while IFS= read -r -d '' sub_dir; do
    # Get the directory name
    dir_name=$(basename "$sub_dir")

    # Find files with the specified extension and store them in an array
    cds_file=($(find "$search_dir" -type f -name "*.cds.fna"))
    fna_file=($(find "$search_dir" -type f -name "*.fna"))
    gff_file=($(find "$search_dir" -type f -name "*.gff"))
    faa_file=($(find "$search_dir" -type f -name "*.faa"))

    echo -e "$dir_name\tcds_file\t$fna_file\t$gff_file\t$faa_file" >> $output_file
      
  done < <(find "$GENOME_DIR/$run_name/ncbi_dataset" -type d -print0)
  # clean up
  mv $GENOME_DIR/$run_name/ncbi_dataset/* $GENOME_DIR/$run_name/
  rm -fr $GENOME_DIR/$run_name/ncbi_dataset* $GENOME_DIR/$run_name/README.md
fi

