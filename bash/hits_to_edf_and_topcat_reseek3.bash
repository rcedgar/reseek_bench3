#!/bin/bash -e

# mode=$1
# truth=$2
# ref=$3

for ref in scop40 scop40x cath40
do
	for mode in fast sensitive
	do
		for truth in family superfamily superfamilyx fold
		do
			./hits_to_edf_and_topcat_reseek3_one.bash $mode $truth $ref
		done
	done
done