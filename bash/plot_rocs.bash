#!/bin/bash -e

#!/bin/bash -

mkdir -p ../plots
cd ../edf

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40x.family \
	--input tm.scop40x.family \
	--input foldseek.scop40x.family \
	--input reseek3.scop40x.family \
	--yscale log \
	--xscale log \
	--ylim 1e-4,1 \
	--xlim 1e-7,1 \
	--output ../plots/roc_scop40x_family.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input foldseek.scop40x.superfamily \
	--input reseek3.scop40x.superfamily \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--output ../plots/roc_scop40x_superfamily.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40x.superfamilyx \
	--input tm.scop40x.superfamilyx \
	--input foldseek.scop40x.superfamilyx \
	--input reseek3.scop40x.superfamilyx \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--output ../plots/roc_scop40x_superfamilyx.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input foldseek.scop40x.fold \
	--input reseek3.scop40x.fold \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--output ../plots/roc_scop40x_fold.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.family \
	--input tm.scop40.family \
	--input foldseek.scop40.family \
	--input reseek3.scop40.family \
	--yscale log \
	--xscale log \
	--ylim 1e-4,1 \
	--xlim 1e-7,1 \
	--output ../plots/roc_scop40_family.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.superfamily \
	--input tm.scop40.superfamily \
	--input foldseek.scop40.superfamily \
	--input reseek3.scop40.superfamily \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--output ../plots/roc_scop40_superfamily.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.superfamilyx \
	--input tm.scop40.superfamilyx \
	--input foldseek.scop40.superfamilyx \
	--input reseek3.scop40.superfamilyx \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--output ../plots/roc_scop40_superfamilyx.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.fold \
	--input tm.scop40.fold \
	--input foldseek.scop40.fold \
	--input reseek3.scop40.fold \
	--yscale log \
	--xscale log \
	--ylim 1e-3,1 \
	--output ../plots/roc_scop40_fold.svg
