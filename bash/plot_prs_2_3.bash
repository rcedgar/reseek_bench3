#!/bin/bash -e

mkdir -p ../plots
cd ../edf

###########################
# family
###########################
python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40x.family \
	--input reseek28_sensitive.scop40x.family \
	--input reseek3_fast.scop40x.family \
	--input reseek3_sensitive.scop40x.family \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40c family" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_family_2_3.svg

python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40.family \
	--input reseek28_sensitive.scop40.family \
	--input reseek3_fast.scop40.family \
	--input reseek3_sensitive.scop40.family \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 family" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_family_2_3.svg

###########################
# superfamily
###########################
python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40x.superfamily \
	--input reseek28_sensitive.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--ylim 0.7,1 \
	--xlim 3,0.5 \
	--title "PR SCOP40c superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_superfamily_2_3.svg

python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40.superfamily \
	--input reseek28_sensitive.scop40.superfamily \
	--input reseek3_fast.scop40.superfamily \
	--input reseek3_sensitive.scop40.superfamily \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_superfamily_2_3.svg

###########################
# superfamilyx
###########################
python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40x.superfamilyx \
	--input reseek28_sensitive.scop40x.superfamilyx \
	--input reseek3_fast.scop40x.superfamilyx \
	--input reseek3_sensitive.scop40x.superfamilyx \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40c superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_superfamilyx_2_3.svg

python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40.superfamilyx \
	--input reseek28_sensitive.scop40.superfamilyx \
	--input reseek3_fast.scop40.superfamilyx \
	--input reseek3_sensitive.scop40.superfamilyx \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_superfamilyx_2_3.svg

###########################
# fold
###########################
python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40x.fold \
	--input reseek28_sensitive.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40c fold" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40x_fold_2_3.svg

python ../py/plot_curve.py \
	--type pr \
	--input reseek28_fast.scop40.fold \
	--input reseek28_sensitive.scop40.fold \
	--input reseek3_fast.scop40.fold \
	--input reseek3_sensitive.scop40.fold \
	--ylim 0.5,1 \
	--xlim 0,0.5 \
	--title "PR SCOP40 fold" \
	--colors ../info/colors.txt \
	--output ../plots/pr_scop40_fold_2_3.svg
