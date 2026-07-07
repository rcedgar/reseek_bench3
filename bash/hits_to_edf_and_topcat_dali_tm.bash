#!/bin/bash -e

# algo=$1
# ref=$2
# score_or_evalue=$3

./hits_to_edf_and_topcat.bash dali scop40 score
./hits_to_edf_and_topcat.bash dali scop40x score
./hits_to_edf_and_topcat.bash tm scop40 score
./hits_to_edf_and_topcat.bash tm scop40x score
