#!/bin/bash -e

outdir=../strumm
rm -rf $outdir/cdn*
mkdir -p $outdir
cd $outdir

fa=../data/true_cdn.fa
bcb=../data/true_cdn.bcb

muscle \
	-align $fa \
	-output cdn.muscleaa.afa

muscle \
	-align $bcb \
	-output cdn.muscle3d.afa

muscle \
	-strumm_build cdn.muscle3d.afa \
	-input $bcb \
	-output cdn.strumm \
	-map cdn.map \
	-seedmsaout cdn_strumm_seed.afa \
	-jalview_features cdn_strumm_seed_core_blocks.jalview

python ../py/prep_jalview.py \
	../data/cdn_motifs.jalview \
	cdn.muscle3d.afa \
	cdn.map \
	-o cdn_strumm_seed_motifs.jalview

python ../py/cdn_map_motifs.py \
	--motifs ../data/cdn_motifs.tsv \
	--orig cdn.muscle3d.afa \
	--map cdn.map \
	--seed cdn_strumm_seed.afa \
	-o cdn_strumm_seed_motif_cols.tsv
