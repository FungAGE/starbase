#!/bin/bash
# liftover annotations

threads=8

GENOME_PATH=$1
REF_CODE=$2

if [[ $USER == "adrianf" ]]; then
  module load bioinfo-tools minimap2
  ROOT_DIR=/home/adrianf/systbio-project-folder
  STARFISH=$ROOT_DIR/bin/starfish
else
  ROOT_DIR=/mnt/sda/johannesson_lab/adrian
  STARFISH=$ROOT_DIR/starfish
fi

GENOME_DIR=$ROOT_DIR/Genomes/$GENOME_PATH
REFERENCE_GENOME=$(find "$GENOME_DIR"/"$REF_CODE"/ -name "*.fna" -type f)
REFERENCE_GTF=$(find "$GENOME_DIR"/"$REF_CODE"/ -name "*.gff3" -type f)

assemblies=$(find "$GENOME_DIR"/ \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \)  -type f | grep -v "$REF_CODE")

cd "$GENOME_DIR" || exit

for assembly in $assemblies
    do
      ASSEMBLY_CODE=$(dirname $assembly | rev | cut -d / -f 1 | rev)
      BASEDIR=$GENOME_DIR/$ASSEMBLY_CODE
      cd "$BASEDIR" || exit

      LIFTOVER_GTF=$BASEDIR/${ASSEMBLY_CODE}.liftover.gff3

      liftoff \
          -p "$threads" \
          -g "$REFERENCE_GTF" \
          -o "$LIFTOVER_GTF" \
          "$assembly" \
          "$REFERENCE_GENOME"

      rm -fr ./intermediate_files
    done