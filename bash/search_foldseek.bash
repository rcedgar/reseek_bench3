#!/bin/bash -e

foldseek=foldseek-v10

mkdir -p ../big_hits
mkdir -p ../time

for ref in scop40 cath40
do
	tmpdir=../tmp_foldseek
	rm -rf $tmpdir
	mkdir -p $tmpdir
	db=../big_foldseek_dbs/$ref/$ref

	/bin/time -vo ../time/foldseek.$ref \
	$foldseek \
		  easy-search \
		  $db \
		  $db \
		  ../big_hits/foldseek.$ref.hits \
		  $tmpdir \
		  --format-output query,target,evalue

	rm -rf $tmpdir
done
