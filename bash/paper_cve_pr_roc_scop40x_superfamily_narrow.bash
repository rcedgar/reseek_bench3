#!/bin/bash -e

cd ../edf

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--panelwidth 3 \
	--nolegend \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c superfamily" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/cve_scop40x_superfamily_narrow.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--panelwidth 3 \
	--nolegend \
	--ylim 0.7,1 \
	--xlim 0,0.4 \
	--title "PR SCOP40c superfamily" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/pr_scop40x_superfamily_narrow.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--panelwidth 3 \
	--nolegend \
	--xscale log \
	--ylim 0,0.4 \
	--xlim 1e-7,1e-4 \
	--title "ROC SCOP40c superfamily" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/roc_scop40x_superfamily_narrow.svg

cd ../paper_plots

svg_stack.py --direction=h \
	cve_scop40x_superfamily_narrow.svg \
	pr_scop40x_superfamily_narrow.svg \
	roc_scop40x_superfamily_narrow.svg \
	> paper_cve_pr_roc_scop40x_superfamily_narrow.svg
ls -lh ../paper_plots/paper_cve_pr_roc_scop40x_superfamily_narrow.svg
