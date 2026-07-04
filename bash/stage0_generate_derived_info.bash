#!/bin/bash -e

mkdir -p ../derived_info

for db in  cath40 scop40 scop40x
do
	python ../py/build_derived_info.py \
		--lookup ../info/$db.lookup  \
		--output ../derived_info/$db.json
done
