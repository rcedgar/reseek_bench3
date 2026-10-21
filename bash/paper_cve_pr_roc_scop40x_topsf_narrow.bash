#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--panelwidth 2 \
	--panelheight 1.6 \
	--nolegend  --notitle --fontsize 6 \
	--yscale log \
	--ylim 0.001,0.1 \
	--xlim 0,0.7 \
	--decimal-yticks \
	--nolegend  --notitle --fontsize 6 \
	--ylabel "Errors per query" \
	--title "CVE SCOP40c topsf" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/cve_scop40x_topsf_narrow.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--panelwidth 2 \
	--panelheight 1.6 \
	--nolegend  --notitle --fontsize 6 \
	--ylim 0.8,1 \
	--xlim 0,0.7 \
	--nolegend  --notitle --fontsize 6 \
	--title "PR SCOP40c topsf" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/pr_scop40x_topsf_narrow.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--panelwidth 2 \
	--panelheight 1.6 \
	--nolegend  --notitle --fontsize 6 \
	--xscale log \
	--ylim 0,0.7 \
	--xlim 0.001,0.1 \
	--decimal-xticks \
	--nolegend  --notitle --fontsize 6 \
	--title "ROC SCOP40c topsf" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/roc_scop40x_topsf_narrow.svg

python ../py/plot_curve.py \
	--legend-only \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/legend_scop40x_topsf_narrow.svg

cd ../paper_plots

svg_stack.py --direction=h \
	cve_scop40x_topsf_narrow.svg \
	pr_scop40x_topsf_narrow.svg \
	roc_scop40x_topsf_narrow.svg \
	> row_cve_pr_roc_scop40x_topsf_narrow.svg

svg_stack.py --direction=v \
	row_cve_pr_roc_scop40x_topsf_narrow.svg \
	legend_scop40x_topsf_narrow.svg \
	> paper_cve_pr_roc_scop40x_topsf_narrow.svg
ls -lh ../paper_plots/paper_cve_pr_roc_scop40x_topsf_narrow.svg
