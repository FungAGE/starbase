#!/bin/bash

STARFISH=$1
species=$2
model=$3

# source $(which activate)
# conda activate "$STARFISH"-env || exit

if [[ $USER == "adrianf" ]]; then
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/CNEFinder/
else
  export PATH=$PATH:/home/adrian/bin/starfish/
  export PATH=$PATH:/home/adrian/bin/starfish/CNEFinder/
fi

# if [[ $SLURM_NTASKS = "" ]]; then
#   threads=2
# else
  threads=8
# fi

# rm -fr elementFinder pairViz
# mkdir elementFinder pairViz

# search for insertions containing at least one predicted tyr (takes ~1min):
"$STARFISH"/starfish insert \
  --separator "_" \
  -T "$threads" \
  -a ../assemblies.txt \
  -d ../blastdb/"$species".assemblies \
  -b geneFinder/"$species"_"$model".goi.bed \
  -i "$model" \
  -x "$species"_"$model" \
  -o elementFinder

# search for flanking repeats around predicted element boundaries:
"$STARFISH"/starfish flank \
  -a ../assemblies.txt \
  -b elementFinder/"$species"_"$model".insert.bed \
  -x "$species"_"$model" \
  -o elementFinder

# summarize the element metadata, identify overlaps, name sites and identify all captain and cargo genes:
"$STARFISH"/starfish summarize \
  --separator "_" \
  -a ../assemblies.txt \
  -b elementFinder/"$species"_"$model".flank.bed \
  -x "$species"_"$model" \
  -o elementFinder \
  -S elementFinder/"$species"_"$model".insert.stats \
  -f elementFinder/"$species"_"$model".flank.singleDR.stats \
  -g "$species"_"$model".consolidatedGFF.txt \
  -A ../ann/gene2emap.txt \
  -t geneFinder/"$species"_"$model".filt_intersect.ids 

# it is strongly recommended to look at an alignment of each element against its 'best' insertion site to manually filter out false positives. Use circos to visualize nucmer alignments (takes ~2min):
"$STARFISH"/starfish pair-viz \
  -m all \
  -t empty \
  -T "$threads" \
  -A nucmer \
  -a ../assemblies.txt \
  -b elementFinder/"$species"_"$model".elements.bed \
  -f elementFinder/"$species"_"$model".flank.singleDR.stats \
  -S elementFinder/"$species"_"$model".elements.named.stats \
  -o pairViz
