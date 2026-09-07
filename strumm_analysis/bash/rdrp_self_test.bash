#!/bin/bash -e

outdir=../rdrp_self_test
mkdir -p $outdir
cd $outdir

muscle \
	-strumm_search ../data/rdrp.bcb \
	-strumm ../strumm/rdrp_abc.strumm \
	-a3m search.a3m \
	-output search.hits

python ../py/rdrp_classify.py \
	search.a3m \
	--selfreport self_test.txt \
	-o search_classify.tsv

cat self_test.txt
