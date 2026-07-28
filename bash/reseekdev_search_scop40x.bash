#!/bin/bash -e

uno=`uname -o`
if [ $uno == Cygwin ] ; then
	os=cygwin
else
	os=linux
fi

if [ ! -x $reseek ] ; then
	echo Not found reseek=$reseek
	exit 1
fi

mkdir -p ../big_hits_scop40x
mkdir -p ../search_log_scop40x

ref=scop40x
db=../big_reseek_dbs/$ref.bcb
if [ ! -s $db ] ; then
	echo Not found db=$db
	exit 1
fi

################################################
rm -rf ../big_hits_scop40x/reseekdev_${os}*
rm -rf ../search_log_scop40x/reseekdev_${os}*
###################################################

for mode in fast sensitive
do
	for truth in family superfamily fold
	do
		name=reseekdev_${os}_${mode}_${truth}
		hits=../big_hits_scop40x/$name

		$reseek \
				-flat_search_kappa $db \
				-$mode \
				-stats $truth \
				-db $db \
				-output $hits \
				-columns query+target+pvalue \
				-log ../search_log_scop40x/$name
	done
done

echo SECONDS=$SECONDS MINUTES=$(($SECONDS/60))
