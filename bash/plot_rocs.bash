#!/bin/bash -e

#!/bin/bash -

mkdir -p ../plots
cd ../edf

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.family \
	--input tm.scop40.family \
	--input foldseek.scop40.family \
	--input reseek3.scop40.family \
	--yscale log \
	--xscale log \
	--xlim 1e-4,1 \
	--ylim 1e-7,1 \
	--output ../plots/roc_scop40_family.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.superfamily \
	--input tm.scop40.superfamily \
	--input foldseek.scop40.superfamily \
	--input reseek3.scop40.superfamily \
	--yscale log \
	--xscale log \
	--xlim 1e-3,1 \
	--output ../plots/roc_scop40_superfamily.svg

python ../py/plot_curve.py \
	--type roc \
	--input dali.scop40.fold \
	--input tm.scop40.fold \
	--input foldseek.scop40.fold \
	--input reseek3.scop40.fold \
	--yscale log \
	--xscale log \
	--xlim 1e-3,1 \
	--output ../plots/roc_scop40_fold.svg
