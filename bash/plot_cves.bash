#!/bin/bash -e

#!/bin/bash -

mkdir -p ../plots
cd ../edf

python ../py/plot_curve.py \
	--type cve \
	--input dali.scop40.superfamily \
	--input tm.scop40.superfamily \
	--input foldseek.scop40.superfamily \
	--input reseek3.scop40.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--output ../plots/cve_scop40_superfamily.svg

python ../py/plot_curve.py \
	--type cve \
	--input dali.scop40.fold \
	--input tm.scop40.fold \
	--input foldseek.scop40.fold \
	--input reseek3.scop40.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--output ../plots/cve_scop40_fold.svg

python ../py/plot_curve.py \
	--type cve \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input foldseek.scop40x.superfamily \
	--input reseek3.scop40x.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--output ../plots/cve_scop40x_superfamily.svg

python ../py/plot_curve.py \
	--type cve \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input foldseek.scop40x.fold \
	--input reseek3.scop40x.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--output ../plots/cve_scop40x_fold.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.cath40.superfamily \
	--input reseek3.cath40.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--output ../plots/cve_cath40_superfamily.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.cath40.fold \
	--input reseek3.cath40.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--output ../plots/cve_cath40_fold.svg
