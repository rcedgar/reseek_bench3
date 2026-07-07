#!/bin/bash -e

# algo=$1
# ref=$2
# score_or_evalue=$3

./hits_to_edf_and_topcat.bash foldseek scop40 evalue
./hits_to_edf_and_topcat.bash foldseek scop40x evalue
./hits_to_edf_and_topcat.bash foldseek cath40 evalue
