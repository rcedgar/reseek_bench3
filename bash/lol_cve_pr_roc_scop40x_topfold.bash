#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prefilter.scop40x.topfold \
	--input lol_exhaustive.scop40x.topfold \
	--yscale log \
	--ylim 0.001,0.1 \
	--xlim 0,0.3 \
	--decimal-yticks \
	--nolegend \
	--ylabel "Errors per query" \
	--title "CVE SCOP40 fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/lol_cve_scop40x_topfold.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prefilter.scop40x.topfold \
	--input lol_exhaustive.scop40x.topfold \
	--ylim 0.8,1 \
	--xlim 0,0.6 \
	--nolegend \
	--title "PR SCOP40 fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/lol_pr_scop40x_topfold.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prefilter.scop40x.topfold \
	--input lol_exhaustive.scop40x.topfold \
	--xscale log \
	--ylim 0,0.4 \
	--xlim 0.001,0.1 \
	--decimal-xticks \
	--nolegend \
	--title "ROC SCOP40 fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/lol_roc_scop40x_topfold.svg

python ../py/plot_curve.py \
	--legend-only \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input lol_prefilter.scop40x.topfold \
	--input lol_exhaustive.scop40x.topfold \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/lol_legend_scop40x_topfold.svg

cd ../paper_plots

svg_stack.py --direction=h \
	lol_cve_scop40x_topfold.svg \
	lol_pr_scop40x_topfold.svg \
	lol_roc_scop40x_topfold.svg \
	> lol_row_cve_pr_roc_scop40x_topfold.svg

svg_stack.py --direction=v \
	lol_row_cve_pr_roc_scop40x_topfold.svg \
	lol_legend_scop40x_topfold.svg \
	> lol_paper_cve_pr_roc_scop40x_topfold.svg
ls -lh ../paper_plots/lol_paper_cve_pr_roc_scop40x_topfold.svg
