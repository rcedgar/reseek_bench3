#!/bin/bash -e

outdir=../raw_precision
mkdir -p $outdir
cd $outdir

for truth in fam sf fold
do
	python ../py/hits_precision.py \
		--input ../big_reseek3_search_scop40x_raw/$truth.hits \
		--truth $truth \
		--output $truth.tsv
done
