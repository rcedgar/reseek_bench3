#!/bin/bash -e

foldseek=foldseek-v10

mkdir -p ../hits
mkdir -p ../time

tmpdir=../tmp_foldseek
rm -rf $tmpdir
mkdir -p $tmpdir
q=../foldseek_dbs/palms/palms
db=$c/data/pdb/foldseek_db/pdb

/bin/time -vo ../time/foldseek_palms_pdb \
$foldseek \
	  easy-search \
	  $db \
	  $db \
	  ../hits/foldseek_palms_pdb.hits \
	  $tmpdir \
	  --format-output query,target,evalue
