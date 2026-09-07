#!/bin/bash -e

dbdir=../big_search_dbs
hitsdir=../search_hits
outdir=../cdn_analyze_hits
mkdir -p $outdir
cd $outdir

for db in afdb50 bfvd esm30
do
	cut -f2 $hitsdir/cdn_$db.hits \
		> $db.labels
	
	reseek \
		-getchains $dbdir/$db.bcb \
		-labels $db.labels \
		-bcb $db.hits.bcb \
		-fasta $db.hits.fasta
	
	muscle \
		-strumm_search $db.hits.bcb \
		-strumm ../strumm/cdn.strumm \
		-pvalue 1e-3 \
		-output $db.cdn_strumm.hits \
		-a3m $db.cdn.a3m
done

for db in afdb50 bfvd esm30
do
	python ../py/cdn_classify.py \
		$db.cdn.a3m \
		--motif_cols ../strumm/cdn_strumm_seed_motif_cols.tsv \
		--cdnfasta $db.cdn.fasta \
		-o $db.cdn_classify.tsv
done

cut -f1,2 afdb50.cdn_classify.tsv | grep "cdn$" | cut -f1 > afdb50_cdn.accs
fgrep -Ff afdb50_cdn.accs $c/data/afdb/acc_desc.txt > afdb50_cdn_acc_desc.tsv
