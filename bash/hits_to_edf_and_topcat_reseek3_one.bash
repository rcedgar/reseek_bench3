#!/bin/bash -e

mode=$1
truth=$2
ref=$3
score_or_evalue=evalue

if [ "$score_or_evalue" != score -a "$score_or_evalue" != evalue ] ; then
	echo "Must set score or evalue"
	exit 1
fi

if [ "$mode" != fast -a "$mode" != sensitive ] ; then
	echo "mode must be fast or sensitive"
	exit 1
fi

algo=reseek3_${mode}
derived_info=../derived_info/$ref.json

if [ ! -s $derived_info ] ; then
	echo Not found derived_info=$derived_info
	exit 1
fi

CPP_BIN=../cpp/bin
CACHE_DIR=../cache/
HITS_TO_EDF=../cpp/bin/hits_to_edf
HITS_TO_TOPCAT=../cpp/bin/hits_to_topcat

if [ ! -x $HITS_TO_EDF -o ! -x $HITS_TO_TOPCAT ]; then
	echo "Build C++ tools first: cd cpp && make"
	exit 1
fi

truthx=$truth
if [ $truth == "superfamilyx" ] ; then
	truthx=superfamily
fi

cache=../cache/${algo}_${truthx}.${ref}.cache
edf=../edf/${algo}.${ref}.$truth

if [ ! -s "$cache" ] ; then
	echo Not found cache=$cache
	exit 1
fi

mkdir -p ../edf
mkdir -p ../topcat

$HITS_TO_EDF \
	--bin "$cache" \
	--$score_or_evalue \
	--truth $truth \
	--derived-info "../derived_info/$ref.json" \
	--algo "$algo" \
	--output $edf

if [ $truth == superfamily ] ; then
	topcat=../topcat/${algo}.${ref}.topsf
	$HITS_TO_TOPCAT \
		--bin "$cache" \
		--$score_or_evalue \
		--truth topsf \
		--derived-info "../derived_info/$ref.json" \
		--algo "$algo" \
		--output $topcat
elif [ $truth == fold ] ; then
	topcat=../topcat/${algo}.${ref}.topfold
	$HITS_TO_TOPCAT \
		--bin "$cache" \
		--$score_or_evalue \
		--truth topfold \
		--derived-info "../derived_info/$ref.json" \
		--algo "$algo" \
		--output $topcat
fi
