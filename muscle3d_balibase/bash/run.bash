#!/bin/bash

outdir=../afa
##############
rm -rf $outdir
##############
mkdir -p $outdir
mkdir -p ../log_align

mkdir -p ../results

for x in `cat ../info/sets.txt`
do
	echo $x
	muscle \
		-quiet \
		-align ../bcb/$x.bcb \
		-output $outdir/$x \
		-log ../log_align/$x.log
done

muscle \
	-qscoredir ../info/sets.txt \
	-testdir ../afa \
	-refdir ../ref \
	-bysequence \
	-log ../results/qscoredir_defaults.log \
	-output ../results/qscoredir_defaults.tsv

grep avgtc= ../results/qscoredir_defaults.tsv
