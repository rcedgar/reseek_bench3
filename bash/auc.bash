#!/bin/bash -e

python ../py/paper_auc.py \
	> ../results/auc.tsv
ls -lh ../results/auc.tsv

python ../py/paper_auc_top.py \
	> ../results/auc_top.tsv
ls -lh ../results/auc_top.tsv
