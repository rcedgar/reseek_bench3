#!/bin/bash -e

mkdir -p ../results
python ../py/report_summary.py \
	../results/summary \
	`ls ../edf/* | grep -v 28` \
	`ls ../topcat/* | grep -v 28`

head ../results/summary.txt
ls -lh ../results/summary.txt
