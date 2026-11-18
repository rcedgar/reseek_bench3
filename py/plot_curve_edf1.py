#!/usr/bin/env python3
"""
plot_curve_edf1.py -- one panel from edf1 quantile tables.

Each --input is a qfract/sens table from hits_to_edf1.py.  The number of
points is taken from the file; # bins= is not used.  Axes are fixed to [0, 1].

Examples:
  python plot_curve_edf1.py --input reseek.superfamilyx.edf1 \
      --input foldseek.superfamilyx.edf1 --output edf1_superfamilyx.svg

  python plot_curve_edf1.py --input reseek.superfamilyx.edf1 \
      --colors algo_styles.txt --notitle --output edf1_superfamilyx.svg
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


def read_edf1(path: str) -> Tuple[dict, List[Tuple[float, float]]]:
    """Load an edf1 comment header and (qfract, sens) rows in file order."""
    hdr: dict[str, str] = {}
    points: List[Tuple[float, float]] = []
    saw_columns = False

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                for k, v in plot_curve.parse_header_line(line[1:].strip()):
                    hdr[k] = v
            elif line.startswith("qfract\t"):
                saw_columns = True
            elif line.strip():
                q_s, sens_s = line.split("\t")
                points.append((float(q_s), float(sens_s)))

    if not saw_columns:
        raise ValueError(f"{path}: missing qfract column header")
    return hdr, points


def fraction_tick_label(value: float) -> str:
    """Label 0 as '0'; every other tick with one decimal, including 1.0."""
    if value == 0:
        return "0"
    return f"{value:.1f}"


def apply_unit_axes(ax) -> None:
    """X ticks every 0.1 and y ticks every 0.2, both spanning 0 to 1."""
    xticks = [i / 10 for i in range(11)]
    yticks = [i / 5 for i in range(6)]
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([fraction_tick_label(x) for x in xticks])
    ax.set_yticklabels([fraction_tick_label(y) for y in yticks])


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot per-query sensitivity quantiles from edf1 tables."
    )
    ap.add_argument(
        "--input",
        action="append",
        required=True,
        metavar="PATH",
        help="edf1 table (qfract, sens); repeat for algorithm overlay",
    )
    ap.add_argument("--output", required=True, help="Output SVG path")
    ap.add_argument("--panelwidth", help="panel width inches", type=float, default=4)
    ap.add_argument("--panelheight", help="panel height inches", type=float, default=3)
    ap.add_argument(
        "--title",
        help="Plot title (default: Sensitivity (truth))",
    )
    ap.add_argument(
        "--notitle",
        action="store_true",
        help="Do not draw a plot title",
    )
    ap.add_argument(
        "--fontsize",
        type=float,
        help="Font size for axis labels and tick labels",
    )
    ap.add_argument("--ylabel", help="Y axis label (default: Sensitivity)")
    ap.add_argument(
        "--colors",
        metavar="PATH",
        help=(
            "Text file with algo_name=color[,lw][,ls][,label] per line; "
            "missing names use defaults"
        ),
    )
    ap.add_argument(
        "--nolegend",
        action="store_true",
        help="Do not draw a legend on the panel",
    )
    ap.add_argument(
        "--legend-only",
        action="store_true",
        help="Write a horizontal legend SVG from --input / --colors (no curves)",
    )
    args = ap.parse_args()

    styles = plot_curve.load_series_styles(args.colors)

    series: List[Tuple[str, List[Tuple[float, float]]]] = []
    truth: Optional[str] = None

    for path in args.input:
        hdr, points = read_edf1(path)
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

    if args.legend_only:
        plot_curve.write_legend_svg([algo for algo, _ in series], styles, args.output)
        sys.stderr.write(f"wrote {args.output}\n")
        return 0

    assert truth is not None

    fig, ax = plt.subplots(figsize=(args.panelwidth, args.panelheight))

    for idx, (label, points) in enumerate(series):
        if not points:
            continue
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        style = plot_curve.series_style(label, idx, styles)
        legend_label = style.legend_label if style.legend_label is not None else label
        ax.plot(
            xs,
            ys,
            color=style.color,
            linewidth=style.linewidth,
            linestyle=style.linestyle,
            label=legend_label,
        )

    apply_unit_axes(ax)

    xlabel = "Fraction of queries"
    ylabel = "Sensitivity" if args.ylabel is None else args.ylabel
    label_kw = {} if args.fontsize is None else {"fontsize": args.fontsize}
    ax.set_xlabel(xlabel, **label_kw)
    ax.set_ylabel(ylabel, **label_kw)
    if args.fontsize is not None:
        ax.tick_params(axis="both", which="both", labelsize=args.fontsize)
    if not args.notitle:
        title = args.title if args.title is not None else f"Sensitivity ({truth})"
        ax.set_title(title)
    ax.grid(True, which="major", linewidth=0.5, alpha=0.5)

    if not args.nolegend and len(series) > 1:
        ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(f"wrote {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
