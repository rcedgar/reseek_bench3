#!/bin/bash -e

outdir=../search_scop_and_cath
mkdir -p $outdir
cd $outdir

for scop in scop cath
do
	reseek \
		-search ../data/rdrp_abc_palmcore.bcb \
		-db $c/data/${scop}40/${scop}40_annot.bcb \
		-fast \
		-stats sf \
		-pvalue 1e-6 \
		-columns query+target+qlo+qhi+tlo+thi+qcovpct+tcovpct+pctid+pvalue \
		-output rdrp_${scop}.hits \
		-log rdrp_${scop}.log

	reseek \
		-search ../data/cdn.bcb \
		-db $c/data/${scop}40/${scop}40_annot.bcb \
		-fast \
		-stats sf \
		-pvalue 1e-6 \
		-columns query+target+qlo+qhi+tlo+thi+qcovpct+tcovpct+pctid+pvalue \
		-output cdn_${scop}.hits \
		-log cdn_${scop}.log
done
