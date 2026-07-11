#!/bin/bash -e

#!/bin/bash -

mkdir -p ../plots
cd ../edf

###########################
# family
###########################
python ../py/plot_curve.py \
	--type cve \
	--input foldseek.cath40.family \
	--input reseek3_fast.cath40.family \
	--input reseek3_sensitive.cath40.family \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE CATH40 family" \
	--colors ../info/colors.txt \
	--output ../plots/cve_cath40_family.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40.family \
	--input dali.scop40.family \
	--input tm.scop40.family \
	--input reseek3_fast.scop40.family \
	--input reseek3_sensitive.scop40.family \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40 family" \
	--output ../plots/cve_scop40_family.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.family \
	--input dali.scop40x.family \
	--input tm.scop40x.family \
	--input reseek3_fast.scop40x.family \
	--input reseek3_sensitive.scop40x.family \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c family" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_family.svg

###########################
# superfamily
###########################
python ../py/plot_curve.py \
	--type cve \
	--input foldseek.cath40.superfamily \
	--input reseek3_fast.cath40.superfamily \
	--input reseek3_sensitive.cath40.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE CATH40 superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/cve_cath40_superfamily.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40.superfamily \
	--input dali.scop40.superfamily \
	--input tm.scop40.superfamily \
	--input reseek3_fast.scop40.superfamily \
	--input reseek3_sensitive.scop40.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40 superfamily" \
	--output ../plots/cve_scop40_superfamily.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.superfamily \
	--input dali.scop40x.superfamily \
	--input tm.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_superfamily.svg

###########################
# superfamilyx
###########################
python ../py/plot_curve.py \
	--type cve \
	--input foldseek.cath40.superfamilyx \
	--input reseek3_fast.cath40.superfamilyx \
	--input reseek3_sensitive.cath40.superfamilyx \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE CATH40 superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/cve_cath40_superfamilyx.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40.superfamilyx \
	--input dali.scop40.superfamilyx \
	--input tm.scop40.superfamilyx \
	--input reseek3_fast.scop40.superfamilyx \
	--input reseek3_sensitive.scop40.superfamilyx \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40 superfamilyx" \
	--output ../plots/cve_scop40_superfamilyx.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.superfamilyx \
	--input dali.scop40x.superfamilyx \
	--input tm.scop40x.superfamilyx \
	--input reseek3_fast.scop40x.superfamilyx \
	--input reseek3_sensitive.scop40x.superfamilyx \
	--yscale log \
	--ylim 1e-3,10 \
	--xlim 0,0.6 \
	--title "CVE SCOP40c superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_superfamilyx.svg

###########################
# fold
###########################
python ../py/plot_curve.py \
	--type cve \
	--input foldseek.cath40.fold \
	--input reseek3_fast.cath40.fold \
	--input reseek3_sensitive.cath40.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.04 \
	--title "CVE CATH40 fold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_cath40_fold.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40.fold \
	--input dali.scop40.fold \
	--input tm.scop40.fold \
	--input reseek3_fast.scop40.fold \
	--input reseek3_sensitive.scop40.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.5 \
	--title "CVE SCOP40 fold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40_fold.svg

python ../py/plot_curve.py \
	--type cve \
	--input foldseek.scop40x.fold \
	--input dali.scop40x.fold \
	--input tm.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--title "CVE SCOP40c fold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_fold.svg
