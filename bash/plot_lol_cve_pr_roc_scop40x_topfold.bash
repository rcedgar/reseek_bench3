#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prob_prefilter.scop40x.topfold \
	--input lol_score_prefilter.scop40x.topfold \
	--yscale log \
	--ylim 0.001,0.1 \
	--xlim 0,0.3 \
	--decimal-yticks \
	--nolegend \
	--ylabel "Errors per query" \
	--title "CVE SCOP40x fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_cve_scop40x_topfold.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prob_prefilter.scop40x.topfold \
	--input lol_score_prefilter.scop40x.topfold \
	--ylim 0.8,1 \
	--xlim 0,0.6 \
	--nolegend \
	--title "PR SCOP40x fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_pr_scop40x_topfold.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prob_prefilter.scop40x.topfold \
	--input lol_score_prefilter.scop40x.topfold \
	--xscale log \
	--ylim 0,0.4 \
	--xlim 0.001,0.1 \
	--decimal-xticks \
	--nolegend \
	--title "ROC SCOP40x fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_roc_scop40x_topfold.svg

python ../py/plot_curve.py \
	--legend-only \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prob_prefilter.scop40x.topfold \
	--input lol_score_prefilter.scop40x.topfold \
	--colors ../info/colors_paper.txt \
	--output ../plots/lol_legend_scop40x_topfold.svg

cd ../plots

svg_stack.py --direction=h \
	lol_cve_scop40x_topfold.svg \
	lol_pr_scop40x_topfold.svg \
	lol_roc_scop40x_topfold.svg \
	> lol_row_cve_pr_roc_scop40x_topfold.svg

svg_stack.py --direction=v \
	lol_row_cve_pr_roc_scop40x_topfold.svg \
	lol_legend_scop40x_topfold.svg \
	> lol_paper_cve_pr_roc_scop40x_topfold.svg
ls -lh ../plots/lol_paper_cve_pr_roc_scop40x_topfold.svg
