#!/bin/bash -e

outdir=../strumm
rm -rf $outdir/rdrp*
mkdir -p $outdir
cd $outdir

fa=../data/rdrp_abc.fa
bcb=../data/rdrp_abc.bcb

muscle \
	-align $fa \
	-output rdrp_abc.muscleaa.afa

muscle \
	-align $bcb \
	-output rdrp_abc.muscle3d.afa

muscle \
	-strumm_build rdrp_abc.muscle3d.afa \
	-input $bcb \
	-output rdrp_abc.strumm \
	-map rdrp_abc.map \
	-seedmsaout rdrp_abc_strumm_seed.afa \
	-jalview_features rdrp_abc_strumm_seed_core_blocks.jalview

python ../py/prep_jalview.py \
	../data/rdrp_motifs.jalview \
	rdrp_abc.muscle3d.afa \
	rdrp_abc.map \
	-o rdrp_abc_strumm_seed_motifs.jalview
