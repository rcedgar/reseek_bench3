#!/bin/bash -e

db=../big_reseek_dbs/pdb.bcb

mkdir -p ../big_hits
mkdir -p ../time
mkdir -p ../logs

for N in 100 1000 # 10000
do
	q=../big_reseek_dbs/pdb_subset$N.bcb
	name=reseek_search_pdbsubset${N}_pdb
	reseek \
		-search $q \
		-db $db \
		-fast \
		-stats sf \
		-output ../big_hits/$name.hits \
		-log ../logs/$name.log
done
