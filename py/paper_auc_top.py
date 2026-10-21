#!/usr/bin/env python3
"""
paper_auc_top.py -- partial AUC for topsf, topsfx, and topfold.

References are scop40x and cath40. Summaries are
<topcat-dir>/<algo>.<reference>.<topsf|topsfx|topfold>.
Axis windows match paper_cve_pr_roc_scop40x_topsf_narrow.svg.
The integrals are the same as paper_auc.py. CVE is area over the curve
(larger is better), including the rectangle from xmin to where the curve
crosses the bottom of the window. PR and ROC stay area under the curve.

cath40 has no tm or dali summaries, and topsfx summaries are not in the
tree. Those rows are NA.

Defaults from the topsf narrow figure:
    CVE  coverage 0–0.7, errors per query 0.001–0.1 (log y)
    PR   recall 0–0.7, precision 0.8–1
    ROC  FPR 0.001–0.1 (log x), TPR 0–0.7

Examples:
    python paper_auc_top.py
    python paper_auc_top.py --cve-xlim 0,0.4 --pr-xlim 0,0.4 --roc-xlim 0.001,0.05
    python paper_auc_top.py --topcat-dir topcat --algo foldseek --algo reseek3_sensitive
"""

from __future__ import annotations

import paper_auc

TOP_TRUTHS = (
    paper_auc.TruthSpec("topsf", "topsf"),
    paper_auc.TruthSpec("topsfx", "topsfx"),
    paper_auc.TruthSpec("topfold", "topfold"),
)

TOPSF = paper_auc.FigureSpec(
    description=(
        "Partial AUC for topsf, topsfx, and topfold on scop40x and cath40. "
        "Axis windows match paper_cve_pr_roc_scop40x_topsf_narrow.svg."
    ),
    dir_flag="topcat-dir",
    dir_help=(
        "Directory of <algo>.<reference>.<truth> summaries (default: <repo>/topcat)"
    ),
    default_dir="topcat",
    algos=paper_auc.DEFAULT_ALGOS,
    references=paper_auc.DEFAULT_REFERENCES,
    truths=TOP_TRUTHS,
    cve_xlim="0,0.7",
    pr_xlim="0,0.7",
    roc_xlim="0.001,0.1",
    cve_ylim="0.001,0.1",
    pr_ylim="0.8,1",
    roc_ylim="0,0.7",
)


def main() -> int:
    return paper_auc.run(TOPSF)


if __name__ == "__main__":
    raise SystemExit(main())
