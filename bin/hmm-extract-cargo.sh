#!/bin/bash

output_path=test.txt
parsed_output=parsed.txt
cpu=4
hmm_path=
fa=blastdb/YRsuperfamRefs.fa

hmmsearch -o $output_path --cpu $cpu --domE 0.001 "$hmm_path" $fa
python bin/hmm.py -i $output_path -q $fa