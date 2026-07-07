#!/bin/bash -e

algo=$1
searchref=$2
ref=$3

# searchref is database used for search (e.g. scop40)
# ref may be same as searchref (scop40) or subset (scop40x)

lookup=../info/$ref.lookup

CPP_BIN=../cpp/bin
CACHE_DIR=../cache/
HITS_TO_BIN=../cpp/bin/hits_to_bin

if [ ! -x $HITS_TO_BIN ]; then
	echo "Build C++ tools first: cd cpp && make"
	exit 1
fi

if [ ! -s $lookup ] ; then
	echo Not found lookup=$lookup >> /dev/stderr
	exit 1
fi

hits=../big_hits/${algo}.${searchref}.hits
cache=../cache/${algo}.${ref}.cache

if [ $cache -nt $hits ] ; then
	echo
	ls -lh $cache $hits
	echo === up to date $cache
	exit 0
fi

mkdir -p $CACHE_DIR

$HITS_TO_BIN \
	--hits $hits \
	--lookup $lookup \
	--fields 1,2,3 \
	--output $cache
