#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek28_sensitive.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--yscale log \
	--ylim 0.001,0.1 \
	--xlim 0,0.7 \
	--decimal-yticks \
	--nolegend \
	--ylabel "Errors per query" \
	--title "CVE SCOP40c topsf" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/cve_scop40x_topsf_with_reseekv2.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek28_sensitive.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--ylim 0.8,1 \
	--xlim 0,1 \
	--nolegend \
	--title "PR SCOP40c topsf" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/pr_scop40x_topsf_with_reseekv2.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek28_sensitive.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--xscale log \
	--ylim 0,0.7 \
	--xlim 0.001,0.1 \
	--decimal-xticks \
	--nolegend \
	--title "ROC SCOP40c topsf" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/roc_scop40x_topsf_with_reseekv2.svg

python ../py/plot_curve.py \
	--legend-only \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek28_sensitive.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/legend_scop40x_topsf_with_reseekv2.svg

cd ../paper_plots

svg_stack.py --direction=h \
	cve_scop40x_topsf_with_reseekv2.svg \
	pr_scop40x_topsf_with_reseekv2.svg \
	roc_scop40x_topsf_with_reseekv2.svg \
	> row_cve_pr_roc_scop40x_topsf_with_reseekv2.svg

svg_stack.py --direction=v \
	row_cve_pr_roc_scop40x_topsf_with_reseekv2.svg \
	legend_scop40x_topsf_with_reseekv2.svg \
	> paper_cve_pr_roc_scop40x_topsf_with_reseekv2.svg
ls -lh ../paper_plots/paper_cve_pr_roc_scop40x_topsf_with_reseekv2.svg
