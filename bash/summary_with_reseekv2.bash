#!/bin/bash -e

mkdir -p ../results
python ../py/report_summary.py \
	../results/summary_with_reseekv2 \
	`ls ../edf/* | fgrep -v embed` \
	`ls ../topcat/* | fgrep -v embed`

head ../results/summary_with_reseekv2.txt
ls -lh ../results/summary_with_reseekv2.txt
