#!/bin/bash -e

cd ../scaling

for x in foldseek reseek_qidx reseek_dbidx reseek_dbidx_shardonly
do
	python ../py/linefit.py \
	--input watson_pdb_subsets_afdb50_${x}_64_threads.tsv \
	--output watson_pdb_subsets_afdb50_${x}_64_threads.svg \
	--xlabel Queries \
	--ylabel Seconds \
	--title $x \
	| tee watson_pdb_subsets_afdb50_${x}_64_threads.fit.txt
done
