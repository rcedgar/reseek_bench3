#!/bin/bash -e

#!/bin/bash -

mkdir -p ../plots
cd ../cross_edf

###########################
# family
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek3_sensitive_family.scop40x.family \
	--input reseek3_sensitive_superfamily.scop40x.family \
	--input reseek3_sensitive_fold.scop40x.family \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--colors ../info/colors.txt \
	--title "CVE SCOP40c family" \
	--output ../plots/cve_scop40x_family_cross.svg

###########################
# superfamily
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek3_sensitive_family.scop40x.superfamily \
	--input reseek3_sensitive_superfamily.scop40x.superfamily \
	--input reseek3_sensitive_fold.scop40x.superfamily \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.4 \
	--title "CVE SCOP40c superfamily" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_superfamily_cross.svg

###########################
# fold
###########################
python ../py/plot_curve.py \
	--type cve \
	--input reseek3_sensitive_family.scop40x.fold \
	--input reseek3_sensitive_superfamily.scop40x.fold \
	--input reseek3_sensitive_fold.scop40x.fold \
	--yscale log \
	--ylim 1e-2,10 \
	--xlim 0,0.6 \
	--title "CVE SCOP40c fold" \
	--colors ../info/colors.txt \
	--output ../plots/cve_scop40x_fold_cross.svg

svg_stack.py --direction=h \
	../plots/cve_scop40x_family_cross.svg \
	../plots/cve_scop40x_superfamily_cross.svg \
	../plots/cve_scop40x_fold_cross.svg \
	> ../plots/cross.svg

ls -lh ../plots/cross.svg