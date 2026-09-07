#!/bin/bash

name=$1
opts=$2

outdir=../afa/
##############
rm -rf $outdir
##############
mkdir -p $outdir
mkdir -p ../qscoredir

for x in `cat ../info/sets.txt`
do
	echo $name $opts $x
	muscle \
		-quiet \
		-align ../bcb/$x.bcb \
		-output $outdir/$x \
		$opts
done

muscle \
	-qscoredir ../info/sets.txt \
	-testdir ../afa \
	-refdir ../ref \
	-bysequence \
	-log $name \
	-output ../qscoredir/$name
