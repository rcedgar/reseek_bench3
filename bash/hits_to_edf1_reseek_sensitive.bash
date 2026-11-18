#!/bin/bash -e

algo=reseek3_sensitive
score_or_evalue=evalue
ref=scop40x

for truth in fold superfamily family superfamilyx
do
	if [ $truth == fold ] ; then
		hits=../big_hits/reseek3_sensitive_fold.scop40.hits
	elif [ $truth == superfamily -o $truth == superfamilyx ] ; then
		hits=../big_hits/reseek3_sensitive_superfamily.scop40.hits
	elif [ $truth == family ] ; then
		hits=../big_hits/reseek3_sensitive_family.scop40.hits
	else
		echo "Invalid truth=$truth"
		exit 1
	fi

	if [ ! -s $hits ] ; then
		echo "Not found hits=$hits"
		exit 1
	fi

	if [ "$score_or_evalue" != score -a "$score_or_evalue" != evalue ] ; then
		echo "Must set score or evalue"
		exit 1
	fi

	derived_info=../derived_info/$ref.json

	if [ ! -s $derived_info ] ; then
		echo Not found derived_info=$derived_info
		exit 1
	fi

	mkdir -p ../edf1

	edf=../edf1/${algo}.${ref}.$truth
	python ../py/hits_to_edf1.py \
		--hits $hits \
		--derived-info "../derived_info/$ref.json" \
		--lookup ../data/scop40x.lookup \
		--algo "$algo" \
		--truth $truth \
		--$score_or_evalue \
		--output $edf
	ls -lh $edf
done