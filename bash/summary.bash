#!/bin/bash -e

mkdir -p ../results
python ../py/report_summary.py \
	../results/summary \
	../edf/* \
	../topcat/*

head ../results/summary.txt
ls -lh ../results/summary.txt
