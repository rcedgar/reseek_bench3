#!/bin/bash -e

cd ../edf

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--input lol_prefilter.scop40x.superfamily \
	--input lol_exhaustive.scop40x.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--nolegend \
	--title "CVE SCOP40 sf (coll.)" \
	--colors ../info/colors.txt \
	--output ../paper_plots/lol_cve_scop40x_superfamily.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--input lol_prefilter.scop40x.superfamily \
	--input lol_exhaustive.scop40x.superfamily \
	--ylim 0.7,1 \
	--xlim 0,0.4 \
	--nolegend \
	--title "PR SCOP40 sf (coll.)" \
	--colors ../info/colors.txt \
	--output ../paper_plots/lol_pr_scop40x_superfamily.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--input lol_prefilter.scop40x.superfamily \
	--input lol_exhaustive.scop40x.superfamily \
	--xscale log \
	--ylim 0,0.4 \
	--xlim 1e-7,1e-4 \
	--nolegend \
	--title "ROC SCOP40 sf (coll.)" \
	--colors ../info/colors.txt \
	--output ../paper_plots/lol_roc_scop40x_superfamily.svg

cd ../paper_plots

svg_stack.py --direction=h \
	lol_cve_scop40x_superfamily.svg \
	lol_pr_scop40x_superfamily.svg \
	lol_roc_scop40x_superfamily.svg \
	> lol_paper_cve_pr_scop40x_superfamily.svg
ls -lh ../paper_plots/lol_paper_cve_pr_scop40x_superfamily.svg