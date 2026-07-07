#!/bin/bash -e

# algo=$1
# searchref=$2
# ref=$3
for ref in scop40 cath40
do
	for mode in fast sensitive
	do
		for truth in family superfamily fold
		do
			./cache_hits.bash reseek3_${mode}_${truth} $ref $ref
			if [ $ref == scop40 ] ; then
				./cache_hits.bash reseek3_${mode}_${truth} scop40 scop40x
			fi
		done
	done
done
