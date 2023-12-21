#!/bin/bash
out_dir=SQL/data/PaperBLAST

# run through the separate gene models for each captain family
cd SQL/data/hmm/captain/ || exit
awk -v RS="//" '{print $0 > ("YRsuperfams." NR ".hmm")}' YRsuperfams.p1-512.hmm
cd ../../../../ || exit

captain_hmm_files=$(find SQL/data/hmm/captain/ -name "YRsuperfams.*.hmm" -type f -size +100b | grep -v "YRsuperfams.p1-512.hmm")

for captain_hmm_file in $captain_hmm_files; do
  out_file=$(basename $captain_hmm_file .hmm)
  # make sure the file ends with "//"
  sed -i '1{/^$/d}' $captain_hmm_file 
  echo "//" >> $captain_hmm_file
  python bin/HMMpaperblast.py $captain_hmm_file $out_dir/$out_file.tsv
done

# now combine results and find the best hit for each hit?
Rscript bin/tophits-HMMpaperblast.R

# run gene models for cargo genes
cargo_hmm_files=$(find SQL/data/hmm/cargo -name "*.hmm")
for cargo_hmm_file in $cargo_hmm_files; do
  out_file=$(basename $cargo_hmm_file .hmm)
  python bin/HMMpaperblast.py $cargo_hmm_file $out_dir/$out_file.tsv
done