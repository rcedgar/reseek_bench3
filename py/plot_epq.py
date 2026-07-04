#!/usr/bin/env python3
"""
plot_epq.py -- Plot FPEPQ (mean false positives per query) against E-value.

FPEPQ at a given score cutoff is cum_fp / ndom.  For an ideal E-value
calibration, E = FPEPQ; a dotted y=x reference line is drawn.

Input: a single .edf file whose scores are E-values (score_direction=lower_better).

Optional --dbsize N: scores are P-values; E-value = N * P.

Both axes use log scale by default.

Examples:
  python plot_epq.py --input reseek.superfamily.edf --output epq.svg

  python plot_epq.py --input reseek.superfamily.edf --output epq.svg \\
      --xlim 1e-10,10 --ylim 1e-10,10

  python plot_epq.py --input reseek.superfamily.edf --output epq.svg \\
      --dbsize 11211
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common
import plot_curve


def read_epq_points(
    path: str, dbsize: Optional[int]
) -> Tuple[dict, List[Tuple[float, float]]]:
    """Return (header, [(evalue, fpepq), ...]) from an .edf file."""
    hdr, kind, body = plot_curve.read_summary(path)
    if kind != "edf":
        raise ValueError(f"{path}: expected .edf, got .tcat")

    direction = hdr.get("score_direction")
    if direction != "lower_better" and dbsize is None:
        raise ValueError(
            f"{path}: score_direction={direction!r}; "
            "E-value calibration requires --evalue scores "
            "or --dbsize for P-value conversion"
        )

    ndom = int(hdr["ndom"])
    if ndom == 0:
        raise ValueError(f"{path}: ndom is 0")

    points: List[Tuple[float, float]] = []
    for line in body:
        score_s, _n_tp, _n_fp, _cum_tp, cum_fp_s = line.split("\t")
        score = float(score_s)
        evalue = score * dbsize if dbsize is not None else score
        fpepq = int(cum_fp_s) / ndom
        if evalue > 0 and fpepq > 0:
            points.append((evalue, fpepq))
    return hdr, points


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot FPEPQ (mean FPs per query) against E-value from an .edf file."
    )
    ap.add_argument("--input", required=True, help="Input .edf file")
    ap.add_argument("--output", required=True, help="Output path (svg, png, pdf)")
    ap.add_argument(
        "--dbsize",
        type=int,
        default=None,
        metavar="N",
        help="Database size; scores are P-values, E-value = N * P",
    )
    ap.add_argument("--xlim", help="X axis limits as LO,HI")
    ap.add_argument("--ylim", help="Y axis limits as LO,HI")
    args = ap.parse_args()

    hdr, points = read_epq_points(args.input, args.dbsize)
    if not points:
        print("error: no positive E-value / FPEPQ points", file=sys.stderr)
        return 1

    algo = common.algo_label_from_header(hdr, args.input)
    truth = hdr.get("truth", "")
    title = f"E-value calibration ({algo}, {truth})" if truth else f"E-value calibration ({algo})"

    xs = [x for x, _ in points]
    ys = [y for _, y in points]

    fig, ax = plt.subplots(figsize=(8, 6))

    xlim = plot_curve.parse_axis_limit(args.xlim)
    ylim = plot_curve.parse_axis_limit(args.ylim)

    if xlim is not None:
        lo, hi = xlim
    else:
        lo, hi = min(min(xs), min(ys)), max(max(xs), max(ys))

    ax.plot([lo, hi], [lo, hi], color="#999999", linewidth=1,
            linestyle=":", label="ideal (E = FPEPQ)", zorder=1)
    ax.plot(xs, ys, color=plot_curve.SERIES_COLORS[0], linewidth=1.5,
            label=algo, zorder=2)

    ax.set_xscale("log")
    ax.set_yscale("log")

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    xlabel = "E-value"
    if args.dbsize is not None:
        xlabel = f"E-value (dbsize={args.dbsize})"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("FPEPQ")
    ax.set_title(title)
    ax.grid(True, which="major", linewidth=0.5, alpha=0.5)
    ax.grid(True, which="minor", linewidth=0.3, alpha=0.3)
    ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(
        f"wrote {args.output} ({len(points)} points, algo={algo})\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
