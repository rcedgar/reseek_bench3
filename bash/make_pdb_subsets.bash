#!/bin/bash -e

db=../big_reseek_dbs/pdb.bcb
mkdir -p ../logs

for N in 100 1000 4000 10000
do
	pdbcaoutdir=../big_pdb_subsets/$N
	bcb=../big_reseek_dbs/pdb_subset$N.bcb
	mkdir -p $pdbcaoutdir

	reseek \
		-bcx_subsample $db \
		-n $N \
		-randseed 1 \
		-bcb $bcb \
		-log ../logs/pdb_subset_$N.log

	reseek \
		-flat_convert $bcb \
		-pdbcaoutdir $pdbcaoutdir \
		-log ../logs/pdb_subset_${N}_caonly.log
done
