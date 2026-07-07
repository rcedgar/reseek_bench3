#!/bin/bash -e

# algo=$1
# searchref=$2
# ref=$3
for ref in scop40 cath40
do
	for mode in fast sensitive
	do
		./cache_hits.bash reseek28_${mode} $ref $ref
		if [ $ref == scop40 ] ; then
			./cache_hits.bash reseek28_${mode} scop40 scop40x
		fi
	done
done
