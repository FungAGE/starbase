#!/bin/bash

STARFISH=$1
species=$2
model=$3

# source $(which activate) || exit
# conda activate "$STARFISH"-env || exit

if [[ $USER == "adrianf" ]]; then
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/
  export PATH=$PATH:/crex/proj/naiss2023-23-128/nobackup/ADRIAN/bin/starfish/CNEFinder/
  TEMP=$SNIC_TMP
else
  export PATH=$PATH:/home/adrian/bin/starfish/
  export PATH=$PATH:/home/adrian/bin/starfish/CNEFinder/
  TEMP=/tmp/starfish
  mkdir "$TEMP"
fi

threads=4

# Gene finder module
# we can de novo annotate any gene we want, as long as we have an HMM file of a predicted domain within that gene and a multifasta of amino acid sequences of that gene (the more predicted sequences the better).
# first, clean up from last time and create a dedicated directory for good housekeeping:
# rm -fr geneFinder

STARFISH_DB=$STARFISH/database

# de novo annotate genes
if [[ $model = "tyr" ]]; then
  HMM=$STARFISH_DB/YRsuperfams.p1-512.hmm
  PROTEINS=$STARFISH_DB/YRsuperfamRefs.faa
else
  HMM=$STARFISH_DB/"$model".hmm
  PROTEINS=$STARFISH_DB/"$model".mycoDB.faa
fi

# BUG: sometimes "$species"_"$model".filt_intersect.gff is not created
"$STARFISH"/starfish annotate --noCheck --force -T "$threads" \
  --separator "_" \
  -x "$species"_"$model" \
  -a ../assemblies.txt \
  -g ../"$species".gff3 \
  -p "$HMM" \
  -P "$PROTEINS" \
  -i "$model" \
  --tempdir "$TEMP" \
  -o geneFinder  2>&1 | tee ./log_file.txt

# now consolidate the newly predicted gene coordinates with the existing gff3:
"$STARFISH"/starfish consolidate --noCheck \
  --separator "_" \
  -o ./ \
  -g ../"$species".gff3 \
  -G geneFinder/"$species"_"$model".filt_intersect.gff

if [[ -f geneFinder/"$species"_"$model".filt_intersect.gff && -s geneFinder/"$species"_"$model".filt_intersect.gff ]]; then
# create a .txt file with the path to the new consolidated gff file:
realpath "$species"_"$model".filt_intersect.consolidated.gff | perl -pe "s/^/""$species""_$model\t/" > "$species"_"$model".consolidatedGFF.txt
fi

if [[ -f "$species"_"$model".consolidatedGFF.txt && -s "$species"_"$model".consolidatedGFF.txt ]]; then
# organize hits into mutually exclusive neighbourhoods separated by at least 10kb to avoid adjacent hits messing up subsequent analyses:
"$STARFISH"/starfish sketch \
  --separator "_" \
  --noCheck \
  -m 10000 \
  -q geneFinder/"$species"_"$model".filt_intersect.ids \
  -g "$species"_"$model".consolidatedGFF.txt \
  -i s \
  -x "$species"_"$model" \
  -o geneFinder/
fi

if [[ -f geneFinder/"$species"_"$model".bed && -s geneFinder/"$species"_"$model".bed ]]; then
# neighbourhoods will often contain intervening genes located between genes of interest so pull out the coordinates of candidate captains only:
grep -P "\t$model\t" geneFinder/"$species"_"$model".bed > geneFinder/"$species"_"$model".goi.bed 
fi
