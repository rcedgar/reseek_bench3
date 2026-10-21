#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--yscale log \
	--ylim 0.001,0.1 \
	--xlim 0,0.4 \
	--decimal-yticks \
	--nolegend \
	--ylabel "Errors per query" \
	--title "CVE SCOP40 fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/cve_scop40x_topfold.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--ylim 0.8,1 \
	--xlim 0.0,0.4 \
	--nolegend \
	--title "PR SCOP40 fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/pr_scop40x_topfold.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--xscale log \
	--ylim 0,0.3 \
	--xlim 0.001,0.1 \
	--decimal-xticks \
	--nolegend \
	--title "ROC SCOP40 fold (cat.)" \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/roc_scop40x_topfold.svg

python ../py/plot_curve.py \
	--legend-only \
	--input foldseek.scop40x.topfold \
	--input dali.scop40x.topfold \
	--input tm.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--colors ../info/colors_paper.txt \
	--output ../paper_plots/legend_scop40x_topfold.svg

cd ../paper_plots

svg_stack.py --direction=h \
	cve_scop40x_topfold.svg \
	pr_scop40x_topfold.svg \
	roc_scop40x_topfold.svg \
	> row_cve_pr_roc_scop40x_topfold.svg

svg_stack.py --direction=v \
	row_cve_pr_roc_scop40x_topfold.svg \
	legend_scop40x_topfold.svg \
	> paper_cve_pr_roc_scop40x_topfold.svg
ls -lh ../paper_plots/paper_cve_pr_roc_scop40x_topfold.svg
