#!/usr/bin/env python3
"""
auc_rank.py -- mean AUC improvement of Reseek3 over each other algorithm.

Reads results/auc.tsv and results/auc_top.tsv. The comparison is
scop40x only, for sf, fold, topsf, and topfold. Fold is computed from
edf/ with the same windows as sf. reseek3_fast is omitted.
reseek3_sensitive is reported as Reseek3.

Each paired cell is one (panel, reference, truth) with a numeric area for
both algorithms. CVE and ROC use auc_log. PR uses auc_linear. The
percentage is 100 * (Reseek3 - other) / other. Pairs whose other area is
0, and NA cells, are left out.

Wilcoxon rank-sum (Mann-Whitney) compares the two area samples within
each panel, one-sided for Reseek3 being larger. Wilcoxon signed-rank and
a one-sample t-test are applied to the pooled percentages, also one-sided
for a positive improvement.

Example:
    python auc_rank.py
    python auc_rank.py --auc results/auc.tsv --auc-top results/auc_top.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import scipy.stats

import plot_auc_summary
import plot_curve

REPO = Path(__file__).resolve().parent.parent
MISSING = {"", "NA", "na"}
DROP_ALGOS = {"reseek3_fast"}
REFERENCE = "reseek3_sensitive"
REFERENCE_LABEL = "Reseek3"
PANELS = ("cve", "pr", "roc")
Cell = Tuple[str, str, str]
Pair = Tuple[str, float, float, float]


def read_tsv(path: Path) -> List[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def area_of(row: dict) -> Optional[float]:
    """CVE and ROC use the log integral. PR uses the linear integral."""
    panel = row["panel"]
    if panel == "cve" or panel == "roc":
        text = row["auc_log"]
    elif panel == "pr":
        text = row["auc_linear"]
    else:
        raise ValueError(f"unknown panel {panel!r}")
    if text.strip() in MISSING:
        return None
    return float(text)


def display_name(algo: str, styles: dict) -> str:
    if algo == REFERENCE:
        return REFERENCE_LABEL
    style = plot_curve.series_style(algo, 0, styles)
    if style.legend_label:
        return style.legend_label
    return algo


def load_cells(rows: Sequence[dict]) -> Tuple[Dict[Cell, Dict[str, float]], List[str]]:
    """Map each cell to algo -> area. Also return other algos in file order."""
    cells: Dict[Cell, Dict[str, float]] = {}
    others: List[str] = []
    for row in rows:
        algo = row["algo"]
        if algo in DROP_ALGOS:
            continue
        value = area_of(row)
        if value is None:
            continue
        cell = (row["panel"], row["reference"], row["truth"])
        slot = cells.setdefault(cell, {})
        if algo in slot:
            raise SystemExit(f"duplicate row for {algo} {cell}")
        slot[algo] = value
        if algo != REFERENCE and algo not in others:
            others.append(algo)
    return cells, others


def pairs_for(cells: Dict[Cell, Dict[str, float]], algo: str) -> List[Pair]:
    """(panel, Reseek3 area, other area, percent improvement)."""
    paired: List[Pair] = []
    for (panel, _reference, _truth), values in cells.items():
        if REFERENCE not in values or algo not in values:
            continue
        other = values[algo]
        if other == 0.0:
            continue
        reseek = values[REFERENCE]
        percent = 100.0 * (reseek - other) / other
        paired.append((panel, reseek, other, percent))
    return paired


def rank_sum_p(pairs: Sequence[Pair], panel: str) -> Optional[float]:
    """One-sided Mann-Whitney: Reseek3 areas greater than the other algorithm."""
    here = [pair for pair in pairs if pair[0] == panel]
    if not here:
        return None
    reseek = [pair[1] for pair in here]
    other = [pair[2] for pair in here]
    result = scipy.stats.mannwhitneyu(reseek, other, alternative="greater")
    return float(result.pvalue)


def signed_rank_p(percents: Sequence[float]) -> Optional[float]:
    """One-sided Wilcoxon signed-rank: percentage improvements greater than 0."""
    nonzero = [value for value in percents if value != 0.0]
    if not nonzero:
        return None
    result = scipy.stats.wilcoxon(nonzero, alternative="greater")
    return float(result.pvalue)


def ttest_p(percents: Sequence[float]) -> Optional[float]:
    """One-sided t-test: mean percentage improvement greater than 0."""
    if len(percents) < 2:
        return None
    if all(value == percents[0] for value in percents):
        return None
    result = scipy.stats.ttest_1samp(percents, 0.0, alternative="greater")
    return float(result.pvalue)


def fmt(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def report_rows(
    cells: Dict[Cell, Dict[str, float]],
    others: Sequence[str],
    styles: dict,
) -> List[str]:
    header = "\t".join(
        (
            "algo",
            "n",
            "mean_pct",
            "rank_sum_cve",
            "rank_sum_pr",
            "rank_sum_roc",
            "signed_rank_p",
            "ttest_p",
        )
    )
    lines = [header]
    for algo in others:
        paired = pairs_for(cells, algo)
        percents = [pair[3] for pair in paired]
        mean_pct = sum(percents) / len(percents) if percents else None
        fields = [
            display_name(algo, styles),
            str(len(paired)),
            fmt(mean_pct),
            fmt(rank_sum_p(paired, "cve")),
            fmt(rank_sum_p(paired, "pr")),
            fmt(rank_sum_p(paired, "roc")),
            fmt(signed_rank_p(percents) if percents else None),
            fmt(ttest_p(percents)),
        ]
        lines.append("\t".join(fields))
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Mean percentage AUC improvement of Reseek3 over each other algorithm, "
            "with rank-sum, signed-rank, and t-test p-values."
        )
    )
    ap.add_argument("--auc", type=Path, default=REPO / "results" / "auc.tsv")
    ap.add_argument("--auc-top", type=Path, default=REPO / "results" / "auc_top.tsv")
    ap.add_argument(
        "--colors",
        type=Path,
        default=REPO / "info" / "colors_paper.txt",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=REPO / "results" / "auc_ranks.txt",
        help="Also write the TSV here (default: results/auc_ranks.txt)",
    )
    args = ap.parse_args()

    for path in (args.auc, args.auc_top):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing or empty table: {path}")

    styles = plot_curve.load_series_styles(str(args.colors))
    wanted = set(plot_auc_summary.SHOW)
    rows = [
        row
        for row in read_tsv(args.auc) + read_tsv(args.auc_top)
        if (row["reference"], row["truth"]) in wanted
    ]
    rows.extend(plot_auc_summary.fold_auc_rows())
    cells, others = load_cells(rows)
    if REFERENCE not in {algo for values in cells.values() for algo in values}:
        raise SystemExit(f"no rows for {REFERENCE}")
    if not others:
        raise SystemExit("no algorithms to compare with Reseek3")

    lines = report_rows(cells, others, styles)
    text = "\n".join(lines) + "\n"
    sys.stdout.write(text)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
