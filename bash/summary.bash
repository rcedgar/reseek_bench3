#!/bin/bash -e

mkdir -p ../results
python ../py/report_summary.py \
	../results/summary \
	`ls ../edf/* | fgrep -v 28` \
	`ls ../topcat/* | fgrep -v 28` \
	`ls ../topcat/* | grep -v 301`

head ../results/summary.txt
ls -lh ../results/summary.txt
