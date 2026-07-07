#!/bin/bash -e

# algo=$1
# searchref=$2
# ref=$3

./cache_hits.bash dali scop40 scop40
./cache_hits.bash dali scop40 scop40x
./cache_hits.bash tm scop40 scop40
./cache_hits.bash tm scop40 scop40x
