#!/bin/bash -e

outdir=../search_hits
mkdir -p $outdir
cd $outdir

refdb=../data/cdn.bcb

for bigdb in afdb50 bfvd esm30
do
	../binaries/reseek \
		-search $refdb \
		-db ../big_search_dbs/$bigdb.bcb \
		-sensitive \
		-stats sf \
		-pvalue 1e-6 \
		-columns query+target+pvalue \
		-output cdn_$bigdb.hits \
		-log cdn_search.$bigdb.log
done
