#!/bin/bash -e

outdir=../big_reseek3_search_scop40x_raw
mkdir -p $outdir
cd $outdir

db=../big_reseek_dbs/scop40x.bcb

for stats in fam sf fold
do
	reseek \
		-search $db \
		-db $db \
		-sensitive \
		-pvalue 1 \
		-stats $stats \
		-columns query+target+raw+p_tp+p_fp \
		-output $stats.hits
done
