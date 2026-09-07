#!/bin/bash -e

dbdir=../big_search_dbs
hitsdir=../search_hits
outdir=../rdrp_analyze_hits
mkdir -p $outdir
cd $outdir

for db in afdb50 bfvd esm30
do
	python $src/palm_annot/py/palm_annot.py \
		--input $db.rdrp.palmdb_notmatchedid90.fasta \
		--seqtype aa \
		--rdrp $db.rdrp.palmdb_notmatchedid90.palm_annot.fasta \
		--fev $db.rdrp.palmdb_notmatchedid90.palm_annot.fev
done
