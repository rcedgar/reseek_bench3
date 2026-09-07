#!/bin/bash -e

outdir=../search_afdb_cdn
mkdir -p $outdir
cd $outdir

db=$c/data/afdb/afdb50.bcb

reseek \
	-search ../data/cdn.bcb \
	-db $db \
	-stats sf \
	-sensitive \
	-columns query+target+pvalue \
	-output search_afdb_cdn.tsv

grep -i Uncharacterized search_afdb_cdn.tsv \
	| cut -f2 \
	| sort \
	| uniq \
	> unch.labels

reseek \
	-getchains c:/data/afdb/afdb50.bcb \
	-labels \
	unch.labels \
	-bcb unch.bcb \
	-cal unch.cal

muscle \
	-strumm_search unch.bcb \
	-strumm ../strumm/cdn.strumm \
	-a3m unch.a3m \
	-pvalue 0.001

python \
	../py/a3m_delete_inserts.py \
	unch.a3m \
	> unch.afa
