#!/bin/bash -e

cd ../plots

mkdir -p ../tmp
rm -f ../tmp/cve_row*_2_3.svg
rm -f ../tmp/topcve_row*_2_3.svg
rm -f ../tmp/roc_row*_2_3.svg
rm -f ../tmp/pr_row*_2_3.svg
rm -f combined_*_2_3.svg

svg_stack.py --direction=h pr_scop40_family_2_3.svg		pr_scop40x_family_2_3.svg		> ../tmp/pr_row1_2_3.svg
svg_stack.py --direction=h pr_scop40_superfamily_2_3.svg   pr_scop40x_superfamily_2_3.svg > ../tmp/pr_row2_2_3.svg
svg_stack.py --direction=h pr_scop40_superfamilyx_2_3.svg  pr_scop40x_superfamilyx_2_3.svg > ../tmp/pr_row3_2_3.svg
svg_stack.py --direction=h pr_scop40_fold_2_3.svg			pr_scop40x_fold_2_3.svg		> ../tmp/pr_row4_2_3.svg

svg_stack.py --direction=h cve_scop40_family_2_3.svg		cve_scop40x_family_2_3.svg		> ../tmp/cve_row1_2_3.svg
svg_stack.py --direction=h cve_scop40_superfamily_2_3.svg   cve_scop40x_superfamily_2_3.svg > ../tmp/cve_row2_2_3.svg
svg_stack.py --direction=h cve_scop40_superfamilyx_2_3.svg  cve_scop40x_superfamilyx_2_3.svg > ../tmp/cve_row3_2_3.svg
svg_stack.py --direction=h cve_scop40_fold_2_3.svg			cve_scop40x_fold_2_3.svg		> ../tmp/cve_row4_2_3.svg

svg_stack.py --direction=h cve_scop40_topsf_2_3.svg   cve_scop40x_topsf_2_3.svg				> ../tmp/topcve_row1_2_3.svg
svg_stack.py --direction=h cve_scop40_topfold_2_3.svg			cve_scop40x_topfold_2_3.svg	> ../tmp/topcve_row2_2_3.svg

svg_stack.py --direction=h roc_scop40_family_2_3.svg		roc_scop40x_family_2_3.svg		> ../tmp/roc_row1_2_3.svg
svg_stack.py --direction=h roc_scop40_superfamily_2_3.svg	roc_scop40x_superfamily_2_3.svg > ../tmp/roc_row2_2_3.svg
svg_stack.py --direction=h roc_scop40_superfamilyx_2_3.svg	roc_scop40x_superfamilyx_2_3.svg> ../tmp/roc_row3_2_3.svg
svg_stack.py --direction=h roc_scop40_fold_2_3.svg			roc_scop40x_fold_2_3.svg		> ../tmp/roc_row4_2_3.svg


svg_stack.py --direction=v ../tmp/cve_row*_2_3.svg > combined_cve_2_3.svg
svg_stack.py --direction=v ../tmp/topcve_row*_2_3.svg > combined_topcve_2_3.svg
svg_stack.py --direction=v ../tmp/roc_row*_2_3.svg > combined_roc_2_3.svg
svg_stack.py --direction=v ../tmp/pr_row*_2_3.svg > combined_pr_2_3.svg
