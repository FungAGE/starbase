#!/bin/bash

###################
# variables
###################

# input
while getopts 'gsfanlm:' OPTION; do
  case "$OPTION" in
    g)
      genus="$OPTARG"
      ;;
    s)
      species="$OPTARG"
      ;;
    f)
      FASTA_IN="$OPTARG"
      ;;
    a)
      GFF_IN="$OPTARG"
      ;;
    n)
      ACCESSION_IN="$OPTARG"
      ;;
    l)
      liftover="yes"
      ;;
    m)
      models="$OPTARG"
      ;;
    ?)
      echo "Usage: $0 [-g <genus_name>] [-s <species_name>] [-f <fasta file>] [ -a <gff file>] [-n <genome assembly acccession>] [-l liftover] [-m model]" >&2
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

# TODO: add input error handling

GENOME_PATH=$genus/$species
genus_species=$genus" "$species
run_name=$genus'_'$species'_starfish'

# starfish parameters
MISSING=1
MAXCOPY=5
PID=90
HSP=1000
FLANK=6

# directories
ROOT_DIR=$(pwd)
TEMP=/tmp/starfish/$run_name
INPUT_DIR=$TEMP/input
FASTA_DIR=$INPUT_DIR/fasta
mkdir -p "$FASTA_DIR"
GFF_DIR=$INPUT_DIR/gff
mkdir -p "$GFF_DIR"
OUTPUT_DIR=$TEMP/output
mkdir -p "$OUTPUT_DIR"

# starfish
SCRIPT_DIR=$ROOT_DIR/bin/starfish
STARFISH=/home/adrian/Systematics/bin/starfish

# starfish environment
conda activate "$ROOT_DIR"/bin/starfish-env

# PATH variables
export PATH=$PATH:$ROOT_DIR/bin/starfish/
export PATH=$PATH:$ROOT_DIR/bin/starfish/CNEFinder/
export EGGNOG_DATA_DIR=/mnt/sda/johannesson_lab/eggnog/

# TODO: take from input and create assemblies.txt
# TODO: make temp dirs for starfish run
ASSEMBLIES=$ROOT_DIR/$run_name/assemblies.txt
GFFS=$ROOT_DIR/$run_name/gffs.txt

find "$FASTA_DIR" -type f -exec sh -c 'echo "$(basename "\$0") $(readlink -f "\$0")"' {} \; > "$ASSEMBLIES"
find "$GFF_DIR" -type f -exec sh -c 'echo "$(basename "\$0") $(readlink -f "\$0")"' {} \; > "$GFFS"

##### dl genomes
# download assemblies and annotations from NCBI
# TODO: or handle $ACCESSION_IN
if [[ -n $ACCESSION_IN ]]; then
  bash "$SCRIPT_DIR"/dl_genomes.sh --accession "$ACCESSION_IN" "$GENOME_PATH" "$run_name"
else
  bash "$SCRIPT_DIR"/dl_genomes.sh --species "$genus_species" "$GENOME_PATH" "$run_name"
fi

# ! not necessary to perform liftover annotations if missing gff file. emapper is the more important step here
##### liftover annotations
if [[ $liftover -eq 1 ]]; then
  echo "
  Starting liftover annotations
  "
  bash "$SCRIPT_DIR"/liftover-annotations.sh "$GENOME_PATH" "$REF_CODE"

  echo "
  liftover annotations finished
  "
fi

##### prep
# format fasta/gffs before running emapper on CDS
echo "
Running starfish_prep.sh
"
bash "$SCRIPT_DIR"/starfish_prep.sh "$species" "$run_name" "$GENOME_PATH"
echo "
finished running starfish_prep.sh
"

# ##### emapper
# * huge difference in required resources between different modes and databases!
# - use diamond and the default $EGGNOG_DATA_DIR

echo "
launching eggnog mapper jobs
"
PROTEOMES=$(find "$ROOT_DIR"/Genomes/"$GENOME_PATH"/*/*.faa)
for PROTEOME in $PROTEOMES
  do
    STRAIN=$(dirname "$PROTEOME" | rev | cut -d / -f 1 | rev)
    for file in $ROOT_DIR/$run_name/emapper/$STRAIN/*.emapper.annotations
      do
        if [ ! -e "$file" ]
        then
          bash "$SCRIPT_DIR"/emapper.sh "$run_name" "$PROTEOME" "$STRAIN"
        fi
      done
  done

# run eggnog mapper for assemblies that do no have proteomes
ASSEMBLIES_MISSING_PROTEOMES=$(find "$ROOT_DIR"/Genomes/"$GENOME_PATH" -type d -exec sh -c 'test -z "$(find "{}" -maxdepth 1 -name "*.faa" -print -quit)"' \; -print | tail -n +2)

for MISSING_PROTEOME in $ASSEMBLIES_MISSING_PROTEOMES
  do
    STRAIN_MISSING_PROTEOMES=$( echo "$MISSING_PROTEOME" | rev | cut -d / -f 1 | rev)
    for file in $ROOT_DIR/$run_name/emapper/$STRAIN_MISSING_PROTEOMES/*.emapper.annotations
      do
        if [ ! -e "$file" ]
        then
          # TODO: get cds sequences rather than running on the whole genome
          GFF=$(find "$MISSING_PROTEOME"/*.gff3 | grep "starfish_format")
          FASTA=$(find "$MISSING_PROTEOME"/*.fna | grep "starfish_format")
          BED=${GFF%.gff3}.bed
          gff2bed < "$GFF" > "$BED"
          bedtools sort -i "$BED" > "${BED%.bed}".sorted.bed
          bedtools getfasta -fi "$FASTA" -bed "${BED%.bed}".sorted.bed -name -split -s -fo "${FASTA%.*}".cds.fasta

          bash "$SCRIPT_DIR"/emapper.sh "$run_name" "${FASTA%.*}".cds.fasta "$STRAIN_MISSING_PROTEOMES"
      
          # rm -f "$BED"* "${FASTA%.*}".cds.fasta
        fi
      done
  done

##### prep
echo "
Running starfish_prep_annotations.sh
"
bash "$SCRIPT_DIR"/starfish_prep_annotations.sh "$run_name" "$GENOME_PATH" $MISSING $MAXCOPY
echo "
finished running starfish_prep_annotations.sh
"

##### run
# separated cpu-intesive steps into modules
for model in $models
  do
    mkdir -p "$OUTPUT_DIR"/"$model"
    cd "$OUTPUT_DIR"/"$model" || exit

    echo "
    Starfish annotate....
    "

    bash "$SCRIPT_DIR"/starfish_annotate.sh $STARFISH "$species" "$model"

    echo "
    Starfish Element Finder....
    "

    bash "$SCRIPT_DIR"/starfish_element_finder.sh $STARFISH "$species" "$model"

    echo "
    Starfish Region Finder....
    "

    bash "$SCRIPT_DIR"/starfish_region_finder.sh $STARFISH "$species" "$model" $MAXCOPY $MISSING $FLANK

  done