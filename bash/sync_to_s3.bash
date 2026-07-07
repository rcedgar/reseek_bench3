#!/bin/bash -e

name=reseek_bench3

mkdir -p ../sync_log

now=`date "+%Y-%m-%d-%H_%S"`

date "+%Y-%m-%d-%H_%S" \
	> ../sync_log/$now.started

aws s3 sync .. s3://serratus-rce-mirror/$name \
	--exclude ".git/*" \
	> /tmp/sync.$now.stdout \
	2> /tmp/sync.$now.stderr

mv -v /tmp/sync.$now.stderr \
	../sync_log/

grep upload /tmp/sync.$now.stdout \
	| sed "-es/Completed.*calculating.....//" \
	> ../sync_log/sync.$now.uploads

date "+%Y-%m-%d-%H_%S" \
	> ../sync_log/$now.finished

echo
echo
echo =============================
echo sync $now finished
echo =============================

ls -lh ../sync_log/*$now*
