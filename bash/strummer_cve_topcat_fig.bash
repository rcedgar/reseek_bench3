#!/bin/bash -e

cd "$(dirname "$0")"

Q=/mnt/c/src/strummer_scop/lookup/test.lookup
T=/mnt/c/src/strummer_scop/lookup/train.lookup
MERGED=/mnt/c/src/strummer_scop/merged_hits/merged

mkdir -p ../topcat ../plot

python ../py/strummer_topcat.py \
	--hits "$MERGED" \
	--queries "$Q" \
	--targets "$T" \
	--output ../topcat/strummer.scop40t.topsf

python ../py/scop40t_topcat.py \
	--queries "$Q" \
	--targets "$T" \
	--big-hits ../big_hits \
	--output-dir ../topcat

python ../py/strummer_cve_topcat_fig.py \
	--topcat-dir ../topcat \
	--colors ../info/colors_paper.txt \
	--output ../plot/strummer_cve_topcat.svg

ls -lh ../plot/strummer_cve_topcat.svg
