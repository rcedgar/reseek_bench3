cd ../paper_plots

svg_stack.py --direction=v \
	paper_cve_pr_roc_scop40x_superfamily.svg \
	paper_cve_pr_roc_scop40x_fold.svg \
	> paper_cve_pr_roc_scop40x_superfamily_fold.svg
ls -lh ../paper_plots/paper_cve_pr_roc_scop40x_superfamily_fold.svg
