#!/bin/bash -e

cd z:/int/muscle_benchmark/

outdir=c:/src/muscle3d_balibase/bcb/
####################################
rm -rf $outdir
##############
mkdir $outdir

for x in `cat c:/src/muscle3d_balibase/info/sets.txt`
do
	reseek \
		-convert z:/int/muscle_benchmark/list_files/balibase/$x.files \
		-bcb $outdir/$x.bcb
	ls -lh $outdir/$x.bcb
done
