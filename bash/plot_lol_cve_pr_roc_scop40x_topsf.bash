#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--input lol_prob_prefilter.scop40x.topsf \
	--input lol_score_prefilter.scop40x.topsf \
	--yscale log \
	--ylim 0.001,0.1 \
	--xlim 0,0.7 \
	--decimal-yticks \
	--nolegend \
	--ylabel "Errors per query" \
	--title "CVE SCOP40x sf (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_cve_scop40x_topsf.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--input lol_prob_prefilter.scop40x.topsf \
	--input lol_score_prefilter.scop40x.topsf \
	--ylim 0.8,1 \
	--xlim 0,0.8 \
	--nolegend \
	--title "PR SCOP40x sf (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_pr_scop40x_topsf.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--input lol_prob_prefilter.scop40x.topsf \
	--input lol_score_prefilter.scop40x.topsf \
	--xscale log \
	--ylim 0,0.7 \
	--xlim 0.001,0.1 \
	--decimal-xticks \
	--nolegend \
	--title "ROC SCOP40x sf (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_roc_scop40x_topsf.svg

python ../py/plot_curve.py \
	--legend-only \
	--input foldseek.scop40x.topsf \
	--input dali.scop40x.topsf \
	--input tm.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--input lol_prob_prefilter.scop40x.topsf \
	--input lol_score_prefilter.scop40x.topsf \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_legend_scop40x_topsf.svg

cd ../plots

svg_stack.py --direction=h \
	lol_cve_scop40x_topsf.svg \
	lol_pr_scop40x_topsf.svg \
	lol_roc_scop40x_topsf.svg \
	> lol_row_cve_pr_roc_scop40x_topsf.svg

svg_stack.py --direction=v \
	lol_row_cve_pr_roc_scop40x_topsf.svg \
	lol_legend_scop40x_topsf.svg \
	> lol_paper_cve_pr_roc_scop40x_topsf.svg
ls -lh ../plots/lol_paper_cve_pr_roc_scop40x_topsf.svg
