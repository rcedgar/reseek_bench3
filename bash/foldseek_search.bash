#!/bin/bash -e

foldseek=foldseek-v10

mkdir -p ../hits
mkdir -p ../time

for ref in scop40 cath40
do
	tmpdir=../tmp_foldseek
	rm -rf $tmpdir
	mkdir -p $tmpdir
	db=../foldseek_dbs/$ref

	/bin/time -vo ../time/foldseek.$ref \
	$foldseek \
		  easy-search \
		  $db \
		  $db \
		  ../hits/foldseek.$ref \
		  $tmpdir \
		  --format-output query,target,evalue

	rm -rf $tmpdir
done
