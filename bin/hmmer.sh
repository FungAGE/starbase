#!/bin/bash
HMM=/home/adrian/Systematics/bin/starfish/database/YRsuperfams.p1-512.hmm

# prepare
hmmpress $HMM

# search protein database with protein query
jackhmmer MADE1.hmm protein_target.fa

# search nucl database with nucl query
nhmmer MADE1.hmm dna_target.fa