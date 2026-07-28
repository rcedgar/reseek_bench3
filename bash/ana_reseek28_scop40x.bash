#!/bin/bash -e

uno=`uname -o`
if [ $uno == Cygwin ] ; then
	os=cygwin
else
	os=linux
fi

mkdir -p ../ana_log_scop40x

rm -f ../ana_log_scop40x/reseek28_${os}*

for mode in fast sensitive
do
	for truth in family superfamily fold
	do
		hitsname=reseek28_${os}_${mode}
		ananame=reseek28_${os}_${mode}_${truth}
		hits=../big_hits_scop40x/$hitsname.hits
		if [ ! -s $hits ] ; then
			echo Not found hits=$hits
			exit 1
		fi

		reseek \
			-fast_bench_hits $hits \
			-truth $truth \
			-log ../ana_log_scop40x/$ananame
	done
done
