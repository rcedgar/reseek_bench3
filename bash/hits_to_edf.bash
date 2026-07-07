#!/bin/bash -e

algo=$1
ref=$2
score_or_evalue=$3

if [ "$score_or_evalue" != score -a "$score_or_evalue" != evalue ] ; then
	echo "Must set score or evalue"
	exit 1
fi

derived_info=../derived_info/$ref.json

if [ ! -s $derived_info ] ; then
	echo Not found derived_info=$derived_info
	exit 1
fi

CPP_BIN=../cpp/bin
CACHE_DIR=../cache/
HITS_TO_EDF=../cpp/bin/hits_to_edf

if [ ! -x $HITS_TO_EDF ]; then
	echo "Build C++ tools first: cd cpp && make"
	exit 1
fi

edf=../edf/${algo}.${ref}
cache=../cache/${algo}.${ref}.cache

if [ $edf -nt $cache ] ; then
	echo === up to date edf
	exit 0
fi

mkdir -p ../edf

$HITS_TO_EDF \
	--bin "$cache" \
	--$score_or_evalue \
	--all-truths \
	--derived-info "../derived_info/$ref.json" \
	--algo "$algo" \
	--output $edf
