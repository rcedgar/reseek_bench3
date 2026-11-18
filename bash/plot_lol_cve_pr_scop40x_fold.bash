#!/bin/bash -e

cd ../edf

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--input lol_prob_prefilter.scop40x.fold \
	--input lol_score_prefilter.scop40x.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0.1,0.5 \
	--nolegend \
	--title "CVE SCOP40x sfx (coll.)" \
	--colors ../info/colors.txt \
	--output ../plots/lol_cve_scop40x_fold.svg

python ../py/plot_curve.py \
	--type pr \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--input lol_prob_prefilter.scop40x.fold \
	--input lol_score_prefilter.scop40x.fold \
	--ylim 0.7,1 \
	--xlim 0.2,0.5 \
	--nolegend \
	--title "PR SCOP40x sfx (coll.)" \
	--colors ../info/colors.txt \
	--output ../plots/lol_pr_scop40x_fold.svg

python ../py/plot_curve.py \
	--type roc \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--input lol_prob_prefilter.scop40x.fold \
	--input lol_score_prefilter.scop40x.fold \
	--xscale log \
	--ylim 0,0.5 \
	--xlim 1e-7,1e-4 \
	--nolegend \
	--title "ROC SCOP40x sfx (coll.)" \
	--colors ../info/colors.txt \
	--output ../plots/lol_roc_scop40x_fold.svg

cd ../plots

svg_stack.py --direction=h \
	lol_cve_scop40x_fold.svg \
	lol_pr_scop40x_fold.svg \
	lol_roc_scop40x_fold.svg \
	> lol_paper_cve_pr_scop40x_fold.svg
ls -lh ../plots/lol_paper_cve_pr_scop40x_fold.svg