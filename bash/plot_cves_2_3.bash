#!/bin/bash -e

#!/bin/bash -

mkdir -p ../plots
cd ../edf

###########################
# family
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40.family \
	--input reseek28_sensitive.scop40.family \
	--input reseek3_fast.scop40.family \
	--input reseek3_sensitive.scop40.family \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40 family" \
	--output ../plots/cve_scop40_family_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40x.family \
	--input reseek28_sensitive.scop40x.family \
	--input reseek3_fast.scop40x.family \
	--input reseek3_sensitive.scop40x.family \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c family" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_family_2_3.svg

###########################
# superfamily
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40.superfamily \
	--input reseek28_sensitive.scop40.superfamily \
	--input reseek3_fast.scop40.superfamily \
	--input reseek3_sensitive.scop40.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40 superfamily" \
	--output ../plots/cve_scop40_superfamily_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40x.superfamily \
	--input reseek28_sensitive.scop40x.superfamily \
	--input reseek3_fast.scop40x.superfamily \
	--input reseek3_sensitive.scop40x.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_superfamily_2_3.svg

###########################
# superfamilyx
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40.superfamilyx \
	--input reseek28_sensitive.scop40.superfamilyx \
	--input reseek3_fast.scop40.superfamilyx \
	--input reseek3_sensitive.scop40.superfamilyx \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40 superfamilyx" \
	--output ../plots/cve_scop40_superfamilyx_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40x.superfamilyx \
	--input reseek28_sensitive.scop40x.superfamilyx \
	--input reseek3_fast.scop40x.superfamilyx \
	--input reseek3_sensitive.scop40x.superfamilyx \
	--yscale log \
	--ylim 1e-3,10 \
	--xlim 0,0.6 \
	--title "CVE SCOP40c superfamilyx" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_superfamilyx_2_3.svg

###########################
# fold
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40.fold \
	--input reseek28_sensitive.scop40.fold \
	--input reseek3_fast.scop40.fold \
	--input reseek3_sensitive.scop40.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.5 \
	--title "CVE SCOP40 fold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40_fold_2_3.svg

python ../py/plot_curve.py \
	--type cve \
	--input reseek28_fast.scop40x.fold \
	--input reseek28_sensitive.scop40x.fold \
	--input reseek3_fast.scop40x.fold \
	--input reseek3_sensitive.scop40x.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--title "CVE SCOP40c fold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_fold_2_3.svg
