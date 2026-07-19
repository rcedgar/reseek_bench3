#!/bin/bash -e

for N in 80 100 150 250 1000 10000
do
	pdbdir=../big_pdb_subsets/$N
	outdir=../big_foldseek_dbs/$N
	rm -rf $outdir
	mkdir -p $outdir
	foldseek createdb $pdbdir $outdir/$N
done
