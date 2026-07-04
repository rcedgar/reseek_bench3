#!/bin/bash -e

./stage1a_hits_to_edf.bash
./stage1b_hits_to_topcat.bash

mkdir -p ../results
python ../py/report_summary.py \
	../results/summary \
	../edf/*

head ../results/summary.txt
ls -lh ../results/summary.txt
echo SECONDS=$SECONDS