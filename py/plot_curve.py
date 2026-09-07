#!/usr/bin/env python3
"""
plot_curve.py -- Stage 2: one panel from .edf or .tcat summary files.

Overlays multiple algorithms when each --input shares the same truth standard
and curve type.  ROC plots use FPR on the x-axis and TPR on the y-axis.
Metric definitions are in truth_standards.md.

Examples:
  python plot_curve.py --type roc --input reseek.superfamily.edf \
      --input foldseek.superfamily.edf --output roc_superfamily.svg

  python plot_curve.py --type cve --input reseek.superfamily.edf \
      --xscale log --xlim 1e-4,0.1 --output cve_superfamily.svg

  python plot_curve.py --type cve --input reseek.topfold.tcat \
      --output cve_topfold.svg

  python plot_curve.py --type roc --input reseek.topsf.tcat \
      --output roc_topsf.svg

  python plot_curve.py --type pr --input reseek.topsf.tcat \
      --output pr_topsf.svg

  python plot_curve.py --type roc --input reseek.superfamily.edf \
      --input foldseek.superfamily.edf --title "SCOP40 superfamily" \
      --colors algo_styles.txt --output roc_superfamily.svg

  # algo_styles.txt: reseek3=#d62728,2,dashed,ReSeek3
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, NamedTuple, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

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

DEFAULT_LINEWIDTH = 1.5


class SeriesStyle(NamedTuple):
    color: str
    linewidth: float
    linestyle: str
    legend_label: Optional[str]


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
            x = cum_fp / n_possible_fp if n_possible_fp else 0.0
            y = coverage
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


def tcat_curve_points(
    hdr: dict, body: List[str], curve_type: CurveType
) -> List[Tuple[float, float]]:
    n_q = int(hdr["N_possible_tp"])
    scores = parse_tcat_scores(body)
    scores_are_evalues = scores_are_evalues_from_header(hdr)
    if curve_type == "cve":
        return hits_to_topcat.topcat_cve_curve(scores, n_q, scores_are_evalues)
    if curve_type == "roc":
        return hits_to_topcat.topcat_roc_curve(scores, n_q, scores_are_evalues)
    if curve_type == "pr":
        return hits_to_topcat.topcat_pr_curve(scores, n_q, scores_are_evalues)
    raise ValueError(f"unknown curve type: {curve_type!r}")


def curve_points_for_input(
    path: str, curve_type: CurveType
) -> Tuple[dict, SummaryKind, List[Tuple[float, float]]]:
    hdr, kind, body = read_summary(path)
    if kind == "edf":
        points = edf_curve_points(hdr, body, curve_type)
    else:
        points = tcat_curve_points(hdr, body, curve_type)
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
        return "FPR", "TPR"
    if curve_type == "pr":
        return "Recall", "Precision"
    raise ValueError(f"unknown curve type: {curve_type!r}")


def plot_title(curve_type: CurveType, truth: str) -> str:
    return f"{curve_type.upper()} ({truth})"


def parse_series_style(spec: str) -> SeriesStyle:
    """Parse color[,lw][,ls][,label]; lw defaults to 1, ls to solid."""
    parts = [p.strip() for p in spec.split(",", 3)]
    color = parts[0]
    if not color:
        raise ValueError("expected color[,lw][,ls][,label]")
    lw_s = parts[1] if len(parts) > 1 else ""
    ls_s = parts[2] if len(parts) > 2 else ""
    legend_s = parts[3] if len(parts) > 3 else ""
    linewidth = float(lw_s) if lw_s else 1.0
    linestyle = ls_s if ls_s else "solid"
    legend_label = legend_s if legend_s else None
    return SeriesStyle(color, linewidth, linestyle, legend_label)


def load_series_styles(path: Optional[str]) -> Dict[str, SeriesStyle]:
    """Load algo_name=color[,lw][,ls][,label]; extra or missing names ignored."""
    styles: Dict[str, SeriesStyle] = {}
    if path is None:
        return styles
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                raise ValueError(
                    f"{path}:{lineno}: expected algo_name=color[,lw][,ls][,label]"
                )
            name, spec = line.split("=", 1)
            name = name.strip()
            spec = spec.strip()
            if not name or not spec:
                raise ValueError(
                    f"{path}:{lineno}: expected algo_name=color[,lw][,ls][,label]"
                )
            try:
                styles[name] = parse_series_style(spec)
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: {exc}") from exc
    return styles


def series_style(label: str, idx: int, styles: Dict[str, SeriesStyle]) -> SeriesStyle:
    if label in styles:
        return styles[label]
    return SeriesStyle(
        SERIES_COLORS[idx % len(SERIES_COLORS)],
        DEFAULT_LINEWIDTH,
        "solid",
        None,
    )


def format_decimal_tick(x: float, _pos: Optional[int] = None) -> str:
    """Tick label as a decimal (0.001, 0.01, 0.1), not scientific notation."""
    if x == 0:
        return "0"
    if abs(x) >= 1:
        if float(x).is_integer():
            return str(int(x))
        return f"{x:g}"
    return f"{x:.10f}".rstrip("0").rstrip(".")


def apply_decimal_ticks(axis, is_log: bool) -> None:
    axis.set_major_formatter(FuncFormatter(format_decimal_tick))
    if is_log:
        axis.set_major_locator(LogLocator(base=10))
        axis.set_minor_formatter(NullFormatter())


def display_label(algo: str, idx: int, styles: Dict[str, SeriesStyle]) -> str:
    style = series_style(algo, idx, styles)
    return style.legend_label if style.legend_label is not None else algo


def write_legend_svg(
    series_algos: List[str],
    styles: Dict[str, SeriesStyle],
    output: str,
) -> None:
    """Horizontal legend strip matching three (4 in, 3 in) panels stacked in a row."""
    handles: List[Line2D] = []
    labels: List[str] = []
    for idx, algo in enumerate(series_algos):
        style = series_style(algo, idx, styles)
        label = display_label(algo, idx, styles)
        handles.append(
            Line2D(
                [0],
                [0],
                color=style.color,
                linewidth=style.linewidth,
                linestyle=style.linestyle,
            )
        )
        labels.append(label)

    fig = plt.figure(figsize=(12, 0.55))
    fig.legend(
        handles,
        labels,
        loc="center",
        ncol=max(len(handles), 1),
        frameon=False,
        fontsize=10,
    )
    fig.savefig(output)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot one CVE / ROC / PR panel from homval summary files."
    )
    ap.add_argument(
        "--type",
        choices=("cve", "roc", "pr"),
        help="Curve type to plot (not required with --legend-only)",
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
    ap.add_argument(
        "--title",
        help="Plot title (default: CURVE_TYPE (truth standard))",
    )
    ap.add_argument("--ylabel", help="Y axis label (default depends on curve type)")
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
    ap.add_argument(
        "--decimal-xticks",
        action="store_true",
        help="Format x-axis ticks as decimals (0.001 not 10^-3)",
    )
    ap.add_argument(
        "--decimal-yticks",
        action="store_true",
        help="Format y-axis ticks as decimals (0.001 not 10^-3)",
    )
    args = ap.parse_args()

    if not args.legend_only and args.type is None:
        ap.error("--type is required unless --legend-only is set")

    styles = load_series_styles(args.colors)

    series: List[Tuple[str, List[Tuple[float, float]]]] = []
    truth: Optional[str] = None
    kind: Optional[SummaryKind] = None

    for path in args.input:
        if args.legend_only:
            hdr, file_kind, _body = read_summary(path)
            points: List[Tuple[float, float]] = []
        else:
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

        algo = common.algo_label_from_header(hdr, path)
        series.append((algo, points))
        sys.stderr.write(f"loaded {path} ({len(points)} points, algo={algo})\n")

    if args.legend_only:
        write_legend_svg([algo for algo, _ in series], styles, args.output)
        sys.stderr.write(f"wrote {args.output}\n")
        return 0

    assert truth is not None and kind is not None

    fig, ax = plt.subplots(figsize=(4, 3))

    for idx, (label, points) in enumerate(series):
        if not points:
            continue
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        style = series_style(label, idx, styles)
        legend_label = style.legend_label if style.legend_label is not None else label
        ax.plot(
            xs,
            ys,
            color=style.color,
            linewidth=style.linewidth,
            linestyle=style.linestyle,
            label=legend_label,
        )

    ax.set_xscale(args.xscale)
    ax.set_yscale(args.yscale)

    xlim = parse_axis_limit(args.xlim)
    ylim = parse_axis_limit(args.ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if args.decimal_xticks:
        apply_decimal_ticks(ax.xaxis, args.xscale == "log")
    if args.decimal_yticks:
        apply_decimal_ticks(ax.yaxis, args.yscale == "log")

    xlabel, ylabel = axis_labels(args.type, kind)
    if args.ylabel is not None:
        ylabel = args.ylabel
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    title = args.title if args.title is not None else plot_title(args.type, truth)
    ax.set_title(title)
    ax.grid(True, which="major", linewidth=0.5, alpha=0.5)

    if not args.nolegend and (len(series) > 1 or args.type == "roc"):
        ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(f"wrote {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
