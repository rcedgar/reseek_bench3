#!/bin/bash -e

foldseek=foldseek-v10

mkdir -p ../big_hits
mkdir -p ../time

tmpdir=../tmp_foldseek
db=$c/data/pdb/foldseek_db/pdb

for n in 100 1000
do
	rm -rf $tmpdir
	mkdir -p $tmpdir
	
	q=../big_foldseek_dbs/$n/$n

	/bin/time -vo ../time/foldseek_pdb_subset$n \
	$foldseek \
		  easy-search \
		  $db \
		  $db \
		  ../big_hits/foldseek_pdb_subset$n.hits \
		  $tmpdir \
		  --format-output query,target,evalue
done