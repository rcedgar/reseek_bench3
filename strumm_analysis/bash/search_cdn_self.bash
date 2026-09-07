#!/bin/bash -e

outdir=../search_cdn_self
mkdir -p $outdir
cd $outdir

muscle \
	-strumm_search ../data/cdn.bcb \
	-strumm ../strumm/cdn.strumm \
	-a3m self.a3m \
	-output self.hits \
	-pvalue 0.01

python \
	../py/a3m_delete_inserts.py \
	self.a3m \
	> self.afa
