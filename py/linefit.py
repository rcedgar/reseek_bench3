#!/usr/bin/env python3
"""
linefit.py -- Scatterplot with least-squares straight-line fit.

Input is a TSV with no header; first two columns are x and y.

Examples:
  python linefit.py --input points.tsv --output fit.svg

  python linefit.py --input points.tsv --output fit.svg --title "Scaling"
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_xy(path: str) -> List[Tuple[float, float]]:
    points: List[Tuple[float, float]] = []
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path}:{lineno}: need at least two columns")
            points.append((float(parts[0]), float(parts[1])))
    return points


def fit_line(points: List[Tuple[float, float]]) -> Tuple[float, float]:
    """Return (slope, intercept) for y = slope * x + intercept."""
    n = len(points)
    if n < 2:
        raise ValueError("need at least two points to fit a line")
    sum_x = sum(p[0] for p in points)
    sum_y = sum(p[1] for p in points)
    sum_xx = sum(p[0] * p[0] for p in points)
    sum_xy = sum(p[0] * p[1] for p in points)
    denom = n * sum_xx - sum_x * sum_x
    if denom == 0:
        raise ValueError("cannot fit line: all x values are identical")
    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n
    return slope, intercept


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Scatterplot with least-squares straight-line fit from TSV x,y."
    )
    ap.add_argument("--input", required=True, help="TSV file, no header; columns x, y")
    ap.add_argument("--output", required=True, help="Output path (svg, png, pdf)")
    ap.add_argument("--title", help="Plot title")
    ap.add_argument("--xlabel", default="x", help="X axis label (default: x)")
    ap.add_argument("--ylabel", default="y", help="Y axis label (default: y)")
    args = ap.parse_args()

    points = read_xy(args.input)
    slope, intercept = fit_line(points)
    if args.title:
        sys.stdout.write(f"{args.title}\n")
    sys.stdout.write(f"n={len(points)} slope={slope:.6g} intercept={intercept:.6g}\n")
    sys.stdout.write(f"seconds = {slope:.6g} * N + {intercept:.6g}\n")

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_lo = min(xs)
    x_hi = max(xs)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(xs, ys, s=36, zorder=2)
    ax.plot(
        [x_lo, x_hi],
        [slope * x_lo + intercept, slope * x_hi + intercept],
        color="#d62728",
        linewidth=1.5,
        zorder=3,
        label=f"y = {slope:.4g}x + {intercept:.4g}",
    )
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel(args.ylabel)
    if args.title is not None:
        ax.set_title(args.title)
    ax.grid(True, which="major", linewidth=0.5, alpha=0.5)
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(f"wrote {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
