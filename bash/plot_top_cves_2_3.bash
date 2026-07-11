#!/bin/bash -e

cd ../topcat

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_sensitive.scop40x.topsf \
	--input reseek28_fast.scop40x.topsf \
	--input reseek3_sensitive.scop40x.topsf \
	--input reseek3_fast.scop40x.topsf \
	--yscale log \
	--ylim 1e-2,1 \
	--xlim 0,1 \
	--title "CVE SCOP40c topsf" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_topsf_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_sensitive.scop40x.topfold \
	--input reseek28_fast.scop40x.topfold \
	--input reseek3_sensitive.scop40x.topfold \
	--input reseek3_fast.scop40x.topfold \
	--yscale log \
	--ylim 1e-2,1 \
	--xlim 0,0.8 \
	--title "CVE SCOP40c topfold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_topfold_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_sensitive.scop40.topsf \
	--input reseek28_fast.scop40.topsf \
	--input reseek3_sensitive.scop40.topsf \
	--input reseek3_fast.scop40.topsf \
	--yscale log \
	--ylim 1e-2,0.1 \
	--xlim 0,0.8 \
	--title "CVE SCOP40 topsf" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40_topsf_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek3_sensitive.scop40.topfold \
	--input reseek3_fast.scop40.topfold \
	--input reseek3_sensitive.scop40.topfold \
	--input reseek3_fast.scop40.topfold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,1 \
	--title "CVE SCOP40 topfold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40_topfold_2_3.svg

#################
# CATH topsf
# N_possible_tp=0
#################
# python ../py/plot_curve.py \
	# --type cve \
	# --input foldseek.cath40.topsf \
	# --input reseek3_sensitive.cath40.topsf \
	# --input reseek3_fast.cath40.topsf \
	# --yscale log \
	# --ylim 1e-2,10 \
	# --xlim 0,1 \
	# --title "CVE CATH40 topsf" \
	# --colors ../info/colors.txt \
	# --output ../plots/cve_cath40_topsf_2_3.svg
