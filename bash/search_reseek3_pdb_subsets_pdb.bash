#!/bin/bash -e

db=../big_reseek_dbs/pdb.bcb

threads=64

mkdir -p ../big_hits
mkdir -p ../time
mkdir -p ../logs

for N in 100 1000 # 10000
do
	q=../big_reseek_dbs/pdb_subset$N.bcb
	name=reseek_search_pdbsubset${N}_afdb50_threads$threads
	/bin/time -v -o ../time/$name.time \
		reseek \
			-flat_search_kappa $q \
			-db $db \
			-threads $threads \
			-stats sf \
			-output ../big_hits/$name.hits \
			-log ../logs/$name.log
done
