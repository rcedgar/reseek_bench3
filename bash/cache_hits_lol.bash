#!/bin/bash -e

# algo=$1
# searchref=$2
# ref=$3

for x in prefilter exhaustive
do
	./cache_hits.bash lol_$x scop40 scop40
	./cache_hits.bash lol_$x scop40 scop40x
	# ./cache_hits.bash lol_$x cath40 cath40
done