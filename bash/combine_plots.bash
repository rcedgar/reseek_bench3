#!/bin/bash -e

cd ../plots

mkdir -p ../tmp
rm -f ../tmp/cve_row*.svg
rm -f ../tmp/topcve_row*.svg
rm -f ../tmp/cath_row*.svg
rm -f ../tmp/roc_row*.svg
rm -f ../tmp/pr_row*.svg
rm -f combined_*.svg

svg_stack.py --direction=h pr_scop40_family.svg		pr_scop40x_family.svg		> ../tmp/pr_row1.svg
svg_stack.py --direction=h pr_scop40_superfamily.svg   pr_scop40x_superfamily.svg > ../tmp/pr_row2.svg
svg_stack.py --direction=h pr_scop40_superfamilyx.svg  pr_scop40x_superfamilyx.svg > ../tmp/pr_row3.svg
svg_stack.py --direction=h pr_scop40_fold.svg			pr_scop40x_fold.svg		> ../tmp/pr_row4.svg

svg_stack.py --direction=h cve_scop40_family.svg		cve_scop40x_family.svg		> ../tmp/cve_row1.svg
svg_stack.py --direction=h cve_scop40_superfamily.svg   cve_scop40x_superfamily.svg > ../tmp/cve_row2.svg
svg_stack.py --direction=h cve_scop40_superfamilyx.svg  cve_scop40x_superfamilyx.svg > ../tmp/cve_row3.svg
svg_stack.py --direction=h cve_scop40_fold.svg			cve_scop40x_fold.svg		> ../tmp/cve_row4.svg

svg_stack.py --direction=h cve_scop40_topsf.svg   cve_scop40x_topsf.svg				> ../tmp/topcve_row1.svg
svg_stack.py --direction=h cve_scop40_topfold.svg			cve_scop40x_topfold.svg	> ../tmp/topcve_row2.svg

svg_stack.py --direction=h cve_cath40_superfamily.svg   cve_cath40_superfamilyx.svg	> ../tmp/cath_row1.svg
svg_stack.py --direction=h cve_cath40_fold.svg			cve_cath40_topfold.svg		> ../tmp/cath_row2.svg

svg_stack.py --direction=h roc_scop40_family.svg		roc_scop40x_family.svg		> ../tmp/roc_row1.svg
svg_stack.py --direction=h roc_scop40_superfamily.svg	roc_scop40x_superfamily.svg > ../tmp/roc_row2.svg
svg_stack.py --direction=h roc_scop40_superfamilyx.svg	roc_scop40x_superfamilyx.svg> ../tmp/roc_row3.svg
svg_stack.py --direction=h roc_scop40_fold.svg			roc_scop40x_fold.svg		> ../tmp/roc_row4.svg


svg_stack.py --direction=v ../tmp/cve_row*.svg > combined_cve.svg
svg_stack.py --direction=v ../tmp/topcve_row*.svg > combined_topcve.svg
svg_stack.py --direction=v ../tmp/cath_row*.svg > combined_cath.svg
svg_stack.py --direction=v ../tmp/roc_row*.svg > combined_roc.svg
svg_stack.py --direction=v ../tmp/pr_row*.svg > combined_pr.svg
