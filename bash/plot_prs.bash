#!/bin/bash -e

mkdir -p ../plots
cd ../edf

###########################
# family
###########################
python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.family \
	--input dali.scop40x.family \
	--input tm.scop40x.family \
	--input reseek3_fast.scop40x.family \
	--input reseek3_sensitive.scop40x.family \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40c family" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_family.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40.family \
	--input dali.scop40.family \
	--input tm.scop40.family \
	--input reseek3_fast.scop40.family \
	--input reseek3_sensitive.scop40.family \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 family" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_family.svg

###########################
# superfamily
###########################
python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--ylim 0.7,1 \
	--xlim 3,0.5 \
	--title "PR SCOP40c superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_superfamily.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40.superfamily \
	--input dali.scop40.superfamily \
	--input tm.scop40.superfamily \
	--input reseek3_fast.scop40.superfamily \
	--input reseek3_sensitive.scop40.superfamily \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_superfamily.svg

###########################
# superfamilyx
###########################
python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.superfamilyx \
	--input dali.scop40x.superfamilyx \
	--input tm.scop40x.superfamilyx \
	--input reseek3_fast.scop40x.superfamilyx \
	--input reseek3_sensitive.scop40x.superfamilyx \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40c superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_superfamilyx.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40.superfamilyx \
	--input dali.scop40.superfamilyx \
	--input tm.scop40.superfamilyx \
	--input reseek3_fast.scop40.superfamilyx \
	--input reseek3_sensitive.scop40.superfamilyx \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_superfamilyx.svg

###########################
# fold
###########################
python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40c fold" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_fold.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40.fold \
	--input dali.scop40.fold \
	--input tm.scop40.fold \
	--input reseek3_fast.scop40.fold \
	--input reseek3_sensitive.scop40.fold \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 fold" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_fold.svg
