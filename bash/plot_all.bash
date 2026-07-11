#!/bin/bash -e

./plot_cves.bash
./plot_prs.bash
./plot_rocs.bash
./plot_top_cves.bash
./combine_plots.bash