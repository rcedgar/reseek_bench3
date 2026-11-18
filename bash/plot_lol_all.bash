#!/bin/bash -e

./plot_lol_cve_pr_scop40x_fold.bash
./plot_lol_cve_pr_scop40_fold.bash
./plot_lol_cve_pr_roc_scop40x_topfold.bash

./plot_lol_cve_pr_scop40x_superfamily.bash
./plot_lol_cve_pr_scop40_superfamily.bash
./plot_lol_cve_pr_roc_scop40x_topsf.bash

./plot_lol_cve_pr_scop40x_superfamilyx.bash
./plot_lol_cve_pr_scop40_superfamilyx.bash
