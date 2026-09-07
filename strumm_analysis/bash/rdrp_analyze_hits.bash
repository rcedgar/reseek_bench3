#!/bin/bash -e

dbdir=../big_search_dbs
hitsdir=../search_hits
outdir=../rdrp_analyze_hits
mkdir -p $outdir
cd $outdir

for db in afdb50 bfvd esm30
do
	cut -f2 $hitsdir/rdrp_$db.hits \
		> $db.labels
	
	reseek \
		-getchains $dbdir/$db.bcb \
		-labels $db.labels \
		-bcb $db.hits.bcb \
		-fasta $db.hits.fasta
	
	muscle \
		-strumm_search rdrp_$db.hits.bcb \
		-strumm ../strumm/rdrp_abc.strumm \
		-pvalue 1e-3 \
		-output $db.rdrp_abc_strumm.hits \
		-a3m $db.rdrp_abc.a3m

	python ../py/rdrp_classify.py \
		$db.rdrp_abc.a3m \
		--motif_cols ../strumm/rdrp_abc_strumm_seed_motif_cols.tsv \
		--rdrpfasta $db.rdrp.fasta \
		-o $db.rdrp_classify.tsv

	usearch \
		-usearch_global $db.rdrp.fasta \
		-db $dbdir/palmdb-2023-04-23_sotus.fa \
		-id 0.9 \
		-maxrejects 128 \
		-userout $db.rdrp.palmdb_id90.hits \
		-userfields query+target+id \
		-notmatched $db.rdrp.palmdb_notmatchedid90.fasta
done

rm -f rdrp_collected_palmdb_notmatchedid90.fasta

for db in afdb50 bfvd esm30
do
	cat $db.rdrp.palmdb_notmatchedid90.fasta \
		| sed "-es/^>/>$db./" \
		>> rdrp_collected_palmdb_notmatchedid90.fasta
done

usearch \
	-cluster_fast rdrp_collected_palmdb_notmatchedid90.fasta \
	-id 0.9 \
	-maxrejects 128 \
	-uc rdrp_cluster_collected_palmdb_notmatchedid90.uc \
	-centroids rdrp_novel.id90.fasta
