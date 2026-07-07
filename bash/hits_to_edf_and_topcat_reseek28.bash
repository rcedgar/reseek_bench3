#!/bin/bash -e

# algo=$1
# ref=$2
# score_or_evalue=$3

./hits_to_edf_and_topcat.bash reseek28_fast scop40 evalue
./hits_to_edf_and_topcat.bash reseek28_fast scop40x evalue
./hits_to_edf_and_topcat.bash reseek28_sensitive scop40 evalue
./hits_to_edf_and_topcat.bash reseek28_sensitive scop40x evalue
