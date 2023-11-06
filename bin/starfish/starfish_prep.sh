#!/bin/bash
###################
# variables
###################
species=$1
run_name=$2
GENOME_PATH=$3

###################
# environment
###################

if [[ $USER == "adrianf" ]]; then
  source $(which activate)
  conda activate "$STARFISH"-env || exit
  ROOT_DIR=/home/adrianf/systbio-project-folder
  module load bioinfo-tools blast samtools BEDTools
else
  ROOT_DIR=/mnt/sda/johannesson_lab/adrian
  # conda activate "$STARFISH"-env || exit
fi

STARFISH=$ROOT_DIR/bin/starfish
export PATH=$PATH:$ROOT_DIR/bin/starfish/
export PATH=$PATH:$ROOT_DIR/bin/starfish/CNEFinder/

BASEDIR=$ROOT_DIR/$run_name
cd "$BASEDIR" || exit

# TODO: add a filter for certain species assemblies here? or just be organised with folders in GENOME_DIR?
GENOME_DIR=$ROOT_DIR/Genomes/$GENOME_PATH

export PATH=$PATH:$STARFISH/
export PATH=$PATH:$STARFISH/CNEFinder/

###################
# prepare for starfish
###################
# clean up
find "$GENOME_DIR" \( -name "*.header.*" -o -name "*.starfish_format" \) -delete
# rm -f /*/*.header.* $GENOME_DIR/*/*.starfish_format.*

# make sure dir names don't contain '_'
find $GENOME_DIR -type d -name '*_*' -exec sh -c 'mv "$0" "$(echo "$0" | sed "s/_//g")"' {} \;

find "$GENOME_DIR"/*/ \( -name "*.fna" -o -name "*.fasta" -o -name "*.fa" \) | grep -v "cds_regions" | while read -r FASTA
  do
    header=$(dirname "$FASTA" | rev | cut -d "/" -f 1 | rev)
    header_squished=$(echo $header | sed 's/_//g')
    echo "Adding $header to $FASTA..."
    
    # Check if the header is already present in the file
    if grep -E -q "^>${header}|^>$header_squished" "$FASTA"; then
      awk '{print $1}' "$FASTA" | awk '{sub(/scaffold_/, "scaffold", $1)} 1' | sed -e "s/>${header}/>${header_squished}_/g" -e "s/>${header_squished}/>${header_squished}_/g" > "${FASTA%.fna}".header.fna
    else
      # Add the header if it's not present
      awk '{print $1}' "$FASTA" | awk '{sub(/scaffold_/, "scaffold", $1)} 1' | sed "s/>/>${header_squished}_/g" > "${FASTA%.fna}".header.fna
    fi

    GFF=$(find "$GENOME_DIR"/"$header_squished"/*.gff3)
    echo "adding $header_squished to $GFF..."
    awk '{sub(/scaffold_/, "scaffold", $1)} 1' "$GFF" | sed "/^#/! s/^/${header_squished}_/" | sed 's/jgi.p|Paevar1|//g; s/jgi.p|Paevar_HGY_1|//g' > "${GFF%.gff3}".header.gff3
  done

# create .txt files detailing the absolute path to each genome's gff3 and assembly:
# realpath $GENOME_DIR/*/*.header.fna | perl -pe 's/^(.+?([^\/]+?).fasta)$/\2\t\1/' > assemblies.txt
find "$GENOME_DIR"/* -name "*.header.fna" > assemblies.txt
# realpath $GENOME_DIR/*/*.gff3 | perl -pe 's/^(.+?([^\/]+?).final.gff3)$/\2\t\1/' > gffs.txt
find "$GENOME_DIR"/* -name "*.header.gff3" > gffs.txt

# add columns containing genome code
cat assemblies.txt | rev | cut -d "/" -f2 | rev | sed 's/_//g' > genome-codes.txt

files="assemblies.txt
gffs.txt"

for i in $files
do
  echo "formatting $i..."
  paste genome-codes.txt "$i" > "$i".temp
  mv "$i".temp "$i"
done

echo "formatting gffs to starfish format..."
# will crash if there are extra delimeter characters in fasta headers or 1st column of gff
$STARFISH/starfish format -f assemblies.txt -g gffs.txt -s '_'

# concatenate all gff3 files into a single file (a useful shortcut for some analyses):
rm -f "$BASEDIR"/"$species".gff3
find "$GENOME_DIR"/* -name "*starfish_format.gff3" -exec sh -c "cat {} >> $BASEDIR/$species.gff3" \;

# concatenate all assembly files and make a blastn database:
rm -fr blastdb
mkdir blastdb

# make sure that this file now points to formated fasta's and gff's
sed -i 's/.fna/.starfish_format.fna/g' assemblies.txt
sed -i 's/.gff3/.starfish_format.gff3/g' gffs.txt

# concatenate
find "$GENOME_DIR"/* -name "*.starfish_format.fna" -exec sh -c "cat {} >> blastdb/$species.assemblies.fna" \;
makeblastdb -in blastdb/"$species".assemblies.fna -out blastdb/"$species".assemblies -parse_seqids -dbtype nucl

# calculate %GC content across all genomes (useful for visualizing elements later):
$STARFISH/scripts/seq-gc.sh -Nbw 1000 blastdb/"$species".assemblies.fna > "$species".assemblies.gcContent_w1000.bed
rm -f blastdb/"$species".assemblies.fna