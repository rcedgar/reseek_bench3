#!/bin/bash -e

cd ..

mkdir -p results

python py/scop40x_report.py \
	ana_log_scop40x/* \
	search_log_scop40x/* \
	| tee results/scop40x_report.txt

ls -lh results/scop40x_report.txt
