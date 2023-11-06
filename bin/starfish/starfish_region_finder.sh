#!/bin/bash

STARFISH=$1
species=$2
model=$3
MAXCOPY=$4
MISSING=$5
FLANK=$6

# source $(which activate)
# conda activate "$STARFISH"-env || exit

if [[ $USER == "adrianf" ]]; then
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/CNEFinder/
else
  export PATH=$PATH:/home/adrian/bin/starfish/
  export PATH=$PATH:/home/adrian/bin/starfish/CNEFinder/
fi

threads=8

rm -fr regionFinder locusViz
mkdir regionFinder locusViz

# group all candidates into families using mmseqs2 easy-clust with a very permissive 50% percent ID/ 25% coverage threshold (families with only a single member will automatically be assigned the prefix 'sng'):
mmseqs easy-cluster geneFinder/"$species"_"$model".filt_intersect.fas regionFinder/"$species"_"$model" regionFinder/ \
  --threads "$threads" \
  --min-seq-id 0.5 \
  -c 0.25 \
  --alignment-mode 3 \
  --cov-mode 0 \
  --cluster-reassign

"$STARFISH"/scripts/mmseqs2mclFormat.pl \
  -i regionFinder/"$species"_"$model"_cluster.tsv \
  -g fam \
  -o regionFinder/

# use sourmash and mcl to group all elements into haplotypes based on pairwise k-mer similarities across entire elements:
"$STARFISH"/starfish sim \
  -m element \
  -t nucl \
  -b elementFinder/"$species"_"$model".elements.bed \
  -x "$species" \
  -o regionFinder/ \
  -a ../assemblies.txt

mv regionFinder/"$species".element.nucl.mat regionFinder/"$species"_"$model".element.nucl.mat
mv regionFinder/"$species".element.nucl.sim regionFinder/"$species"_"$model".element.nucl.sim

"$STARFISH"/starfish group \
  -s regionFinder/"$species"_"$model".element.nucl.sim \
  -i hap \
  -o regionFinder/ \
  -t 0.05

# replace captainIDs with elementIDs in the captain groups file:
grep -P '\tcap\t' elementFinder/"$species"_"$model".elements.bed | cut -f4,7 > regionFinder/"$species"_"$model".cap2ship.txt

"$STARFISH"/scripts/searchReplace.pl \
  -i regionFinder/"$species"_"$model"_cluster.mcl \
  -r regionFinder/"$species"_"$model".cap2ship.txt > regionFinder/"$species"_"$model".element_cluster.mcl

# merge captain family with element haplotype info:
"$STARFISH"/scripts/mergeGroupfiles.pl \
  -t regionFinder/"$species"_"$model".element_cluster.mcl \
  -q regionFinder/"$species"_"$model".element.nucl.I1.5.mcl > regionFinder/"$species"_"$model".element.fam-hap.mcl

# create a file with candidates that are not found in any elements (will let us assign them to fragmented haplotypes in the dereplicate analysis):
grep -f <(comm -23 <(cut -f1 geneFinder/"$species"_"$model".filt_intersect.ids | sort) <(grep -P "\tcap\t|\t$model\t" elementFinder/"$species"_"$model".elements.bed | cut -f4| sort)) geneFinder/"$species"_"$model".bed > regionFinder/unaffiliated_"$model".bed

"$STARFISH"/starfish dereplicate \
  --separator "_" \
  -e regionFinder/"$species"_"$model".element.fam-hap.mcl \
  -t regionFinder/unaffiliated_"$model".bed \
  -F elementFinder/"$species"_"$model".elements.feat \
  -S elementFinder/"$species"_"$model".elements.named.stats \
  -O ../ann/gene2og.a"$MISSING".c"$MAXCOPY".txt \
  -g ../gffs.txt \
  -x "$species"_"$model" \
  -o regionFinder/ \
  --flanking "$FLANK" \
  --distance 600000 \
  --mismatching 1 

# I would normally recommend going with the default --flanking 6 but because the Defiant insertion is in a gene-sparse region, it can only be recovered with ---flanking 3

# Element haplotypes consist of predicted mobile element sequences.
# Empty haplotypes consist of a contiguous sequence formed by the flanking regions of an element.
# Fragmented haplotypes consist of a non-empty sequence flanked by the flanking regions of an element but missing a predicted element.
# Coordinates of predicted insertion sites from the elementFinder module are cross referenced with empty and fragmented haplotypes and are considered to be 'verified' if their coordinates overlap.

# Use gggenomes to visualize nucmer alignments (takes ~5min):
"$STARFISH"/starfish locus-viz \
  -T "$threads" \
  -m region-align \
  -a ../assemblies.txt \
  -b elementFinder/"$species"_"$model".elements.bed \
  -x "$species"_"$model" \
  -o locusViz/ \
  -A nucmer \
  -r regionFinder/"$species"_"$model".fog"$FLANK".d600000.m1.regions.txt \
  -d regionFinder/"$species"_"$model".fog"$FLANK".d600000.m1.dereplicated.txt \
  -j regionFinder/"$species"_"$model".fog"$FLANK".d600000.m1.haplotype_jaccard.sim  \
  -g "$species"_"$model".consolidatedGFF.txt \
  --tags geneFinder/"$species"_"$model".filt_intersect.ids \
  --gc ../"$species".assemblies.gcContent_w1000.bed
