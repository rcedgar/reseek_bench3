#!/bin/bash -e

rm -rf ../derived_info
mkdir -p ../derived_info

for ref in scop40 scop40x cath40
do
	../cpp/bin/build_derived_info \
		--lookup ../info/$ref.lookup \
		--output ../derived_info/$ref.json
done

ls -lh ../derived_info/*
