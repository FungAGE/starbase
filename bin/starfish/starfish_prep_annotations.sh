#!/bin/bash
###################
# variables
###################
run_name=$1
GENOME_PATH=$2
MISSING=$3
MAXCOPY=$4

###################
# environment
###################

# source $(which activate)
# conda activate "$STARFISH"-env || exit

if [[ $USER == "adrianf" ]]; then
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/CNEFinder/
  ROOT_DIR=/home/adrianf/systbio-project-folder
  STARFISH=$ROOT_DIR/bin/starfish
  module load bioinfo-tools blast samtools BEDTools
else
  export PATH=$PATH:/home/adrian/bin/starfish/
  export PATH=$PATH:/home/adrian/bin/starfish/CNEFinder/
  ROOT_DIR=/mnt/sda/johannesson_lab/adrian
  STARFISH=$ROOT_DIR/bin/starfish
fi

BASEDIR=$ROOT_DIR/$run_name
cd "$BASEDIR" || exit

# TODO: add a filter for certain species assemblies here? or just be organised with folders in GENOME_DIR?
GENOME_DIR=$ROOT_DIR/Genomes/$GENOME_PATH

# ! eggnog mapper annotations need to have unique prefixes for each genome
EGGNOG_DIR="$BASEDIR/emapper/"

# find "$EGGNOG_DIR" -name "*.annotations" | while read -r file
#   do
#     name=$(basename "$file" .emapper.annotations)
#     name+=_
#     sed "/^[^#]/ s/^gene_/$name/" "$file" |\

#     sed 's/_0::[^\t]*//g'
#      > "${file}".header
#   done

mkdir "$BASEDIR"/ann/

# parse the provided eggnog mapper annotations (NB the format of the output file has changed in more recent emapper versions):
find "$EGGNOG_DIR"/*.annotations.header -exec cut -f1,8 {} \; | grep -v  '#' | perl -pe 's/\t/\tEMAP\t/' | grep -vP '\tNA' > "$BASEDIR"/ann/gene2emap.txt

# retrieve the narrowest eggnog ortholog group per sequence and convert to mcl format:
# Important! Note that the above command only works with older emapper annotation output, such as the output provided as part of this tutorial.
# In order to parse the output from more recent versions of emapper, you must use the following command to retrieve the narrowest eggnog ortholog group per sequence:
# cut -f1,5 "$EGGNOG_DIR"/*/*emapper.annotations | grep -v '#' | perl -pe 's/(^.+?)\t.+,([^,]+)$/\1\t\2/' | perl -pe 's/@/\t/' > "$BASEDIR"/ann/gene2og.txt

# EDIT: there was an extra column at the end, containing taxonomic information. so just kept the first 2 here:
find "$EGGNOG_DIR"/*.annotations.header -exec cut -f1,5 {} \; | grep -v '#' | perl -pe 's/(^.+?)\t.+,([^,]+)$/\1\t\2/' | perl -pe 's/@/\t/' | cut -f 1,2 > "$BASEDIR"/ann/gene2og.txt

# BUG: for this command to run, the second column cannot be empty
# convert to .mcl format:
$STARFISH/scripts/geneOG2mclFormat.pl -i "$BASEDIR"/ann/gene2og.txt -o "$BASEDIR"/ann/

# make sure to keep formatted copies of assembly and annotation files for later
# rm -f "$GENOME_DIR"/*/*.header.fna "$GENOME_DIR"/*/*.header.gff3 "$EGGNOG_DIR"/*emapper.annotations.header

# you can increase confidence in region homology by only looking at gene ortholog groups with low copy numbers missing from few genomes
# It is more useful to play around with copy number thresholds than with genome absence thresholds 
# because dereplicate will automatically filter OGs to retain those that have more taxonomic information i.e., are present in a greater number of individuals.

# ? eggmapper results based on the same genome will have higher similarity of OGs?
$STARFISH/scripts/filterOG.pl \
    -O "$BASEDIR"/ann/gene2og.mcl \
    -a "$MISSING" \
    -c "$MAXCOPY" \
    -o "$BASEDIR"/ann/