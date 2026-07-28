#!/bin/bash -e

uno=`uname -o`
if [ $uno == Cygwin ] ; then
	os=cygwin
	reseek=$src/reseek/github_releases/reseek-v2.8-win64.exe
else
	os=linux
	reseek=$src/reseek/github_releases/reseek-v2.8-linux-x86
fi

if [ ! -x $reseek ] ; then
	echo Not found reseek=$reseek
	exit 1
fi

mkdir -p ../big_hits_scop40x
mkdir -p ../search_log_scop40x

ref=scop40x
db=../big_reseek_dbs/$ref.bca
if [ ! -s $db ] ; then
	echo Not found db=$db
	exit 1
fi

################################################
rm -rf ../big_hits_scop40x/reseek28_${os}*
rm -rf ../search_log_scop40x/reseek28_${os}*
###################################################

for mode in fast sensitive
do
	name=reseek28_${os}_${mode}
	hits=../big_hits_scop40x/$name

	$reseek \
			-search $db \
			-$mode \
			-db $db \
			-output $hits \
			-columns query+target+evalue \
			-log ../logs/$name
done

echo SECONDS=$SECONDS MINUTES=$(($SECONDS/60))
