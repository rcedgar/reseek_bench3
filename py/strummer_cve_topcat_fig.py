#!/usr/bin/env python3
"""
strummer_cve_topcat_fig.py -- CVE panel for Strummer vs structure methods.

Reads topsf .tcat summaries and writes one SVG. Axes match the paper
topsf CVE panel: coverage 0-0.7, errors per query on a log scale from
0.001 to 0.1.

Example:
  python strummer_cve_topcat_fig.py \\
      --topcat-dir ../topcat \\
      --colors ../info/colors_paper.txt \\
      --output ../plot/strummer_cve_topcat.svg
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common
import plot_curve

INPUTS = (
    "foldseek.scop40t.topcat",
    "dali.scop40t.topcat",
    "tm.scop40t.topcat",
    "reseek3_sensitive.scop40t.topcat",
    "strummer.scop40t.topsf",
)

STRUMMER_STYLE = plot_curve.SeriesStyle("blue", 1.5, "solid", "Strummer")


def main() -> int:
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    ap = argparse.ArgumentParser(
        description="Plot the scop40t topsf CVE comparison."
    )
    ap.add_argument(
        "--topcat-dir",
        default=os.path.join(root, "topcat"),
        help="Directory of topsf summary files",
    )
    ap.add_argument(
        "--colors",
        default=os.path.join(root, "info", "colors_paper.txt"),
        help="algo=color[,lw][,ls][,label] file",
    )
    ap.add_argument(
        "--output",
        default=os.path.join(root, "plot", "strummer_cve_topcat.svg"),
        help="Output SVG path",
    )
    args = ap.parse_args()

    styles = plot_curve.load_series_styles(args.colors)
    styles["strummer"] = STRUMMER_STYLE

    series: List[Tuple[str, List[Tuple[float, float]]]] = []
    truth = None
    for name in INPUTS:
        path = os.path.join(args.topcat_dir, name)
        hdr, _kind, points = plot_curve.curve_points_for_input(path, "cve")
        file_truth = hdr.get("truth")
        if file_truth is None:
            raise ValueError(f"{path}: missing truth in header")
        if truth is None:
            truth = file_truth
        elif file_truth != truth:
            raise ValueError(
                f"{path}: truth {file_truth!r} does not match {truth!r}"
            )
        algo = common.algo_label_from_header(hdr, path)
        series.append((algo, points))
        sys.stderr.write(f"loaded {path} ({len(points)} points, algo={algo})\n")

    fig, ax = plt.subplots(figsize=(4, 3))
    for idx, (label, points) in enumerate(series):
        if not points:
            continue
        style = plot_curve.series_style(label, idx, styles)
        legend_label = style.legend_label if style.legend_label is not None else label
        ax.plot(
            [p[0] for p in points],
            [p[1] for p in points],
            color=style.color,
            linewidth=style.linewidth,
            linestyle=style.linestyle,
            label=legend_label,
        )

    ax.set_xscale("linear")
    ax.set_yscale("log")
    ax.set_xlim(0, 0.7)
    ax.set_ylim(0.001, 0.1)
    plot_curve.apply_decimal_ticks(ax.yaxis, True)
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Errors per query")
    ax.set_title("CVE SCOP40t sf (cat.)")
    ax.grid(True, which="major", linewidth=0.5, alpha=0.5)
    ax.legend(fontsize=10)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(f"wrote {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
