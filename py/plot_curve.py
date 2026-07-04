#!/usr/bin/env python3
"""
plot_curve.py -- Stage 2: one panel from .edf or .tcat summary files.

Overlays multiple algorithms when each --input shares the same truth standard
and curve type.  Metric definitions are in truth_standards.md.

Examples:
  python plot_curve.py --type roc --input reseek.superfamily.edf \
      --input foldseek.superfamily.edf --output roc_superfamily.svg

  python plot_curve.py --type cve --input reseek.superfamily.edf \
      --xscale log --xlim 1e-4,0.1 --output cve_superfamily.svg

  python plot_curve.py --type cve --input reseek.topfold.tcat \
      --output cve_topfold.svg
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common
import hits_to_topcat

CurveType = str
SummaryKind = str

SERIES_COLORS = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
)


def parse_header_line(content: str) -> List[Tuple[str, str]]:
    """Parse '# key=value' or '# K1=v1 K2=v2 ...' into (key, value) pairs."""
    pairs: List[Tuple[str, str]] = []
    if content.startswith("SEPQ0.1=") or content.startswith("TEPQ0.001="):
        for part in content.split():
            k, v = part.split("=", 1)
            pairs.append((k, v))
    elif "=" in content:
        k, v = content.split("=", 1)
        pairs.append((k, v))
    return pairs


def read_summary(path: str) -> Tuple[dict, SummaryKind, List[str]]:
    """Load comment header, detect kind (edf|tcat), and return body lines."""
    hdr: dict[str, str] = {}
    kind: Optional[SummaryKind] = None
    body: List[str] = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                for k, v in parse_header_line(line[1:].strip()):
                    hdr[k] = v
            elif line.startswith("score\t"):
                kind = "edf"
            elif line.startswith("domain\t"):
                kind = "tcat"
            elif line.strip():
                body.append(line.strip())

    if kind is None:
        raise ValueError(f"{path}: cannot detect summary kind (missing column header)")
    return hdr, kind, body


def scores_are_evalues_from_header(hdr: dict) -> bool:
    direction = hdr.get("score_direction")
    if direction == "lower_better":
        return True
    if direction == "higher_better":
        return False
    raise ValueError(f"unknown score_direction: {direction!r}")


def parse_tcat_scores(body: List[str]) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    scores: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
    for line in body:
        dom, s_tp_s, s_xf_s = line.split("\t")
        s_tp = None if s_tp_s == "." else float(s_tp_s)
        s_xf = None if s_xf_s == "." else float(s_xf_s)
        scores[dom] = (s_tp, s_xf)
    return scores


def edf_curve_points(
    hdr: dict, body: List[str], curve_type: CurveType
) -> List[Tuple[float, float]]:
    ndom = int(hdr["ndom"])
    n_possible_tp = int(hdr["N_possible_tp"])
    n_possible_fp = int(hdr["N_possible_fp"])

    points: List[Tuple[float, float]] = []
    for line in body:
        _score, _n_tp, _n_fp, cum_tp_s, cum_fp_s = line.split("\t")
        cum_tp = int(cum_tp_s)
        cum_fp = int(cum_fp_s)
        coverage = cum_tp / n_possible_tp if n_possible_tp else 0.0

        if curve_type == "cve":
            x = coverage
            y = cum_fp / ndom if ndom else 0.0
        elif curve_type == "roc":
            x = coverage
            y = cum_fp / n_possible_fp if n_possible_fp else 0.0
        elif curve_type == "pr":
            denom = cum_tp + cum_fp
            if denom == 0:
                continue
            x = coverage
            y = cum_tp / denom
        else:
            raise ValueError(f"unknown curve type: {curve_type!r}")
        points.append((x, y))
    return points


def tcat_cve_curve_points(hdr: dict, body: List[str]) -> List[Tuple[float, float]]:
    n_q = int(hdr["N_possible_tp"])
    scores = parse_tcat_scores(body)
    scores_are_evalues = scores_are_evalues_from_header(hdr)
    return hits_to_topcat.topcat_cve_curve(scores, n_q, scores_are_evalues)


def curve_points_for_input(
    path: str, curve_type: CurveType
) -> Tuple[dict, SummaryKind, List[Tuple[float, float]]]:
    hdr, kind, body = read_summary(path)
    if kind == "edf":
        points = edf_curve_points(hdr, body, curve_type)
    else:
        if curve_type != "cve":
            raise ValueError(
                f"{path}: curve type {curve_type!r} is not defined for .tcat "
                "(only cve is supported)"
            )
        points = tcat_cve_curve_points(hdr, body)
    return hdr, kind, points


def parse_axis_limit(spec: Optional[str]) -> Optional[Tuple[float, float]]:
    if spec is None:
        return None
    lo_s, hi_s = spec.split(",", 1)
    return float(lo_s), float(hi_s)


def axis_labels(curve_type: CurveType, kind: SummaryKind) -> Tuple[str, str]:
    if curve_type == "cve":
        y = "Errors per query" if kind == "edf" else "Error (negative-control FPR)"
        return "Coverage", y
    if curve_type == "roc":
        return "TPR", "FPR"
    if curve_type == "pr":
        return "Recall", "Precision"
    raise ValueError(f"unknown curve type: {curve_type!r}")


def plot_title(curve_type: CurveType, truth: str) -> str:
    return f"{curve_type.upper()} ({truth})"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot one CVE / ROC / PR panel from homval summary files."
    )
    ap.add_argument(
        "--type",
        required=True,
        choices=("cve", "roc", "pr"),
        help="Curve type to plot",
    )
    ap.add_argument(
        "--input",
        action="append",
        required=True,
        metavar="PATH",
        help="Summary file (.edf or .tcat); repeat for algorithm overlay",
    )
    ap.add_argument("--output", required=True, help="Output SVG path")
    ap.add_argument(
        "--xscale",
        choices=("linear", "log"),
        default="linear",
        help="X axis scale (default linear)",
    )
    ap.add_argument(
        "--yscale",
        choices=("linear", "log"),
        default="linear",
        help="Y axis scale (default linear)",
    )
    ap.add_argument("--xlim", help="X axis limits as LO,HI")
    ap.add_argument("--ylim", help="Y axis limits as LO,HI")
    args = ap.parse_args()

    series: List[Tuple[str, List[Tuple[float, float]]]] = []
    truth: Optional[str] = None
    kind: Optional[SummaryKind] = None

    for path in args.input:
        hdr, file_kind, points = curve_points_for_input(path, args.type)
        file_truth = hdr.get("truth")
        if file_truth is None:
            raise ValueError(f"{path}: missing truth in header")
        if truth is None:
            truth = file_truth
            kind = file_kind
        elif file_truth != truth:
            raise ValueError(
                f"{path}: truth {file_truth!r} does not match {truth!r}"
            )
        elif file_kind != kind:
            raise ValueError(f"{path}: cannot mix .edf and .tcat inputs")

        label = common.algo_label_from_header(hdr, path)
        series.append((label, points))
        sys.stderr.write(f"loaded {path} ({len(points)} points, algo={label})\n")

    assert truth is not None and kind is not None

    fig, ax = plt.subplots(figsize=(8, 6))

    for idx, (label, points) in enumerate(series):
        if not points:
            continue
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        color = SERIES_COLORS[idx % len(SERIES_COLORS)]
        ax.plot(xs, ys, color=color, linewidth=1.5, label=label)

    ax.set_xscale(args.xscale)
    ax.set_yscale(args.yscale)

    xlim = parse_axis_limit(args.xlim)
    ylim = parse_axis_limit(args.ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    xlabel, ylabel = axis_labels(args.type, kind)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(plot_title(args.type, truth))
    ax.grid(True, which="major", linewidth=0.5, alpha=0.5)

    if len(series) > 1:
        ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(f"wrote {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
