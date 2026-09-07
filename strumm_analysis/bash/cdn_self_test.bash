#!/bin/bash -e

outdir=../cdn_self_test
mkdir -p $outdir
cd $outdir

muscle \
	-strumm_search ../data/cdn.bcb \
	-strumm ../strumm/cdn.strumm \
	-pvalue 0.001 \
	-a3m search.a3m \
	-output search.hits

python ../py/cdn_classify.py \
	search.a3m \
	--motif_cols ../strumm/cdn_strumm_seed_motif_cols.tsv \
	--selfreport self_test.txt \
	-o search_classify.tsv

cat self_test.txt
