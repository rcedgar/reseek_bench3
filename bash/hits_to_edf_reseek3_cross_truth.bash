#!/bin/bash -e

truth_algo=$1
truth_eval=$2

if [ -z "$truth_eval" ] ; then
	echo Missing arg
	exit 1
fi

algo=reseek3_sensitive_$truth_algo

ref=scop40x
cache=../cache/${algo}.${ref}.cache
edf=../cross_edf/${algo}.${ref}.$truth_eval

mkdir -p ../cross_edf

if [ ! -s $cache ] ; then
	./cache_hits.bash $algo $ref $ref
fi

../cpp/bin/hits_to_edf \
	--bin $cache \
	--evalue \
	--truth $truth_eval \
	--derived-info ../derived_info/$ref.json \
	--algo $algo \
	--output $edf
