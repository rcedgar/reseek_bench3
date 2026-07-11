#!/bin/bash -e

mkdir -p ../plots
cd ../edf

###########################
# family
###########################
python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40x.family \
	--input reseek28_sensitive.scop40x.family \
	--input reseek3_fast.scop40x.family \
	--input reseek3_sensitive.scop40x.family \
	--yscale log \
	--xscale log \
	--ylim 1e-2,1 \
	--xlim 1e-7,1e-2 \
	--title "ROC SCOP40c family" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40x_family_2_3.svg

python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40.family \
	--input reseek28_sensitive.scop40.family \
	--input reseek3_fast.scop40.family \
	--input reseek3_sensitive.scop40.family \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--xlim 1e-7,1e-2 \
	--title "ROC SCOP40 family" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40_family_2_3.svg

###########################
# superfamily
###########################
python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40x.superfamily \
	--input reseek28_sensitive.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--yscale log \
	--xscale log \
	--ylim 1e-2,1 \
	--xlim 1e-7,1e-2 \
	--title "ROC SCOP40c superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40x_superfamily_2_3.svg

python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40.superfamily \
	--input reseek28_sensitive.scop40.superfamily \
	--input reseek3_fast.scop40.superfamily \
	--input reseek3_sensitive.scop40.superfamily \
	--yscale log \
	--xscale log \
	--ylim 1e-1,1 \
	--xlim 1e-7,1e-2 \
	--title "ROC SCOP40 superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40_superfamily_2_3.svg

###########################
# superfamilyx
###########################
python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40x.superfamilyx \
	--input reseek28_sensitive.scop40x.superfamilyx \
	--input reseek3_fast.scop40x.superfamilyx \
	--input reseek3_sensitive.scop40x.superfamilyx \
	--yscale log \
	--xscale log \
	--ylim 1e-2,1 \
	--xlim 1e-8,1e-3 \
	--title "ROC SCOP40c superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40x_superfamilyx_2_3.svg

python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40.superfamilyx \
	--input reseek28_sensitive.scop40.superfamilyx \
	--input reseek3_fast.scop40.superfamilyx \
	--input reseek3_sensitive.scop40.superfamilyx \
	--yscale log \
	--xscale log \
	--ylim 1e-2,1 \
	--xlim 1e-8,1e-3 \
	--title "ROC SCOP40 superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40_superfamilyx_2_3.svg

###########################
# fold
###########################
python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40x.fold \
	--input reseek28_sensitive.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--yscale log \
	--xscale log \
	--ylim 1e-2,1 \
	--xlim 1e-6,1e-2 \
	--title "ROC SCOP40c fold" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40x_fold_2_3.svg

python ../py/plot_curve.py \
	--type roc \
	--input reseek28_fast.scop40.fold \
	--input reseek28_sensitive.scop40.fold \
	--input reseek3_fast.scop40.fold \
	--input reseek3_sensitive.scop40.fold \
	--yscale log \
	--xscale log \
	--ylim 1e-2,1 \
	--xlim 1e-6,1e-2 \
	--title "ROC SCOP40 fold" \
	--colors ../info/colors.txt \
	--output ../plots/roc_scop40_fold_2_3.svg
