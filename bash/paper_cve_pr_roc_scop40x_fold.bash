#!/bin/bash -e

cd ../edf

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c fold" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/cve_scop40x_fold.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--ylim 0.7,1 \
	--xlim 0,0.4 \
	--title "PR SCOP40c fold" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/pr_scop40x_fold.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--xscale log \
	--ylim 0,0.4 \
	--xlim 1e-7,1e-4 \
	--title "ROC SCOP40c fold" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/roc_scop40x_fold.svg

cd ../paper_plots

svg_stack.py --direction=h \
	cve_scop40x_fold.svg \
	pr_scop40x_fold.svg \
	roc_scop40x_fold.svg \
	> paper_cve_pr_roc_scop40x_fold.svg
ls -lh ../paper_plots/paper_cve_pr_roc_scop40x_fold.svg
