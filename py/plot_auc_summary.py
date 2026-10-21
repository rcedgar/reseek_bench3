#!/usr/bin/env python3
"""
plot_auc_summary.py -- one horizontal bar chart per reference and truth.

Reads results/auc.tsv and results/auc_top.tsv. Every stored area is
higher-is-better. The figure is scop40x only: sf, fold, topsf, topfold.
Fold is computed from edf/ with the same windows as sf.

Reseek-fast is omitted. Each bar is stacked CVE, then PR, then ROC.
CVE uses auc_log only, PR uses auc_linear, ROC uses both integrals.
Before plotting, each of those three areas is divided by the best value
of that area in the same chart, then the three fractions are averaged so
the bar axis runs from 0 to 1.
Light, medium, and dark shades of the algorithm color mark the three
panels. Panel height follows the number of bars, so every bar has the
same thickness.

Each chart is sorted with the largest total at the top. The table gives
that total, how far the best score in the chart is ahead, and the mean
rank across the four charts. Numbers use three significant figures.
NA cells are left out of the rank.

Examples:
    python plot_auc_summary.py
    python plot_auc_summary.py --output paper_plots/auc_summary.svg
"""

from __future__ import annotations

import argparse
import colorsys
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

import common
import paper_auc
import plot_curve

REPO = Path(__file__).resolve().parent.parent
MISSING = {"", "NA", "na"}
DROP_ALGOS = {"reseek3_fast"}
PANEL_INDEX = {"cve": 0, "pr": 1, "roc": 2}
# Charts in display order. Fold is pairwise; topsf and topfold are top-hit.
SHOW = (
    ("scop40x", "sf"),
    ("scop40x", "fold"),
    ("scop40x", "topsf"),
    ("scop40x", "topfold"),
)
# HLS lightness and saturation scale for CVE, PR, ROC.
SHADE_LIGHTNESS = (0.90, 0.56, 0.20)
SHADE_SATURATION = (0.35, 0.72, 1.0)
FRACTION_MAX = 1.0
# Bar thickness is half the previous 0.144 in. The gap between bars is 0.018 in.
BAR_HEIGHT = 0.80
INCH_PER_BAR = 0.090
PANEL_GAP_IN = 0.15
TOP_IN = 0.26
BOTTOM_IN = 0.22
# Bar axis is half the previous 2.69 in. The table keeps its width.
FIG_W = 4.18
LABEL_IN = 0.58
BAR_AXIS_IN = 1.34
TABLE_IN = 2.15
FONT = 6.0
TITLE_FONT = 7.0
TABLE_COLUMNS = ((1 / 6, "Total"), (0.5, "Best ahead %"), (5 / 6, "Mean rank"))
Cell = Tuple[str, str]
Parts = Tuple[float, float, float]
ScoreKey = Tuple[str, str, str]


def read_tsv(path: Path) -> List[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fold_auc_rows() -> List[dict]:
    """scop40x fold areas, same axis windows as the superfamily figure."""
    spec = paper_auc.SUPERFAMILY
    args = argparse.Namespace(
        cve_xlim=spec.cve_xlim,
        pr_xlim=spec.pr_xlim,
        roc_xlim=spec.roc_xlim,
        cve_ylim=spec.cve_ylim,
        pr_ylim=spec.pr_ylim,
        roc_ylim=spec.roc_ylim,
    )
    panels = paper_auc.build_panels(args)
    summary_dir = REPO / "edf"
    rows: List[dict] = []
    for algo_name in spec.algos:
        if algo_name in DROP_ALGOS:
            continue
        path = summary_dir / f"{algo_name}.scop40x.fold"
        if not path.is_file():
            raise SystemExit(f"missing fold summary: {path}")
        curves = {}
        for panel in panels:
            hdr, _kind, points = plot_curve.curve_points_for_input(str(path), panel.name)
            algo = common.algo_label_from_header(hdr, algo_name)
            curves[panel.name] = (algo, points)
        for panel in panels:
            algo, points = curves[panel.name]
            linear, log_auc = paper_auc.auc_for_panel(panel, points, algo)
            xmin, xmax = panel.xlim
            ymin, ymax = panel.ylim
            rows.append(
                {
                    "panel": panel.name,
                    "reference": "scop40x",
                    "truth": "fold",
                    "algo": algo,
                    "xmin": paper_auc.fmt(xmin),
                    "xmax": paper_auc.fmt(xmax),
                    "ymin": paper_auc.fmt(ymin),
                    "ymax": paper_auc.fmt(ymax),
                    "auc_linear": paper_auc.fmt(linear),
                    "auc_log": log_auc,
                }
            )
    return rows


def rows_for_figure(pairwise: Sequence[dict], top: Sequence[dict]) -> List[dict]:
    """sf from the pairwise table, fold from edf, topsf and topfold from the top table."""
    wanted = set(SHOW)
    chosen = [
        row
        for row in list(pairwise) + list(top)
        if (row["reference"], row["truth"]) in wanted
    ]
    chosen.extend(fold_auc_rows())
    return chosen


def parse_auc(text: str) -> Optional[float]:
    if text.strip() in MISSING:
        return None
    return float(text)


def row_value(row: dict) -> Optional[float]:
    """Panel contribution. CVE is the log area only."""
    if row["algo"] in DROP_ALGOS:
        return None
    panel = row["panel"]
    if panel == "cve":
        return parse_auc(row["auc_log"])
    if panel == "pr":
        return parse_auc(row["auc_linear"])
    if panel == "roc":
        parts = [
            value
            for value in (parse_auc(row["auc_linear"]), parse_auc(row["auc_log"]))
            if value is not None
        ]
        if not parts:
            return None
        return sum(parts)
    raise ValueError(f"unknown panel {panel!r}")


def cell_parts(rows: Sequence[dict]) -> Dict[ScoreKey, Parts]:
    """CVE, PR, and ROC contributions for each reference, truth, and algorithm."""
    acc: Dict[ScoreKey, List[float]] = {}
    for row in rows:
        value = row_value(row)
        if value is None:
            continue
        key = (row["reference"], row["truth"], row["algo"])
        slot = acc.setdefault(key, [0.0, 0.0, 0.0])
        slot[PANEL_INDEX[row["panel"]]] += value
    return {key: (slot[0], slot[1], slot[2]) for key, slot in acc.items()}


def as_totals(parts: Dict[ScoreKey, Parts]) -> Dict[ScoreKey, float]:
    return {key: sum(value) for key, value in parts.items()}


def drop_zero_cells(scores: Dict[ScoreKey, float]) -> Dict[ScoreKey, float]:
    """Omit a reference/truth whose every score is zero (empty query set)."""
    cells = {(ref, truth) for ref, truth, _algo in scores}
    dropped = []
    for ref, truth in cells:
        values = [
            value
            for (cell_ref, cell_truth, _algo), value in scores.items()
            if cell_ref == ref and cell_truth == truth
        ]
        if values and max(values) == 0.0:
            dropped.append((ref, truth))
    if not dropped:
        return scores
    for ref, truth in dropped:
        sys.stderr.write(f"skipping {ref} {truth} (all scores are 0)\n")
    return {
        key: value
        for key, value in scores.items()
        if key[:2] not in set(dropped)
    }


def drop_zero_parts(parts: Dict[ScoreKey, Parts]) -> Dict[ScoreKey, Parts]:
    kept = drop_zero_cells(as_totals(parts))
    return {key: parts[key] for key in kept}


def fraction_of_best(parts: Dict[ScoreKey, Parts]) -> Dict[ScoreKey, Parts]:
    """Divide each CVE, PR, and ROC area by the best value in that chart."""
    cells = {(ref, truth) for ref, truth, _algo in parts}
    scaled: Dict[ScoreKey, Parts] = {}
    for ref, truth in cells:
        keys = [key for key in parts if key[0] == ref and key[1] == truth]
        best = [0.0, 0.0, 0.0]
        for key in keys:
            for index in range(3):
                if parts[key][index] > best[index]:
                    best[index] = parts[key][index]
        for key in keys:
            scaled[key] = tuple(
                parts[key][index] / best[index] if best[index] > 0.0 else 0.0
                for index in range(3)
            )
    return scaled


def cell_list(rows: Sequence[dict], scores: Dict[ScoreKey, float]) -> List[Cell]:
    """The four scop40x charts, in display order, that have scores."""
    del rows
    return [cell for cell in SHOW if any(key[:2] == cell for key in scores)]


def average_ranks(values: Sequence[float]) -> List[float]:
    """1 = largest. Ties share the mean of the ranks they span."""
    order = sorted(range(len(values)), key=lambda i: values[i], reverse=True)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        shared = (start + 1 + end) / 2.0
        for pos in range(start, end):
            ranks[order[pos]] = shared
        start = end
    return ranks


def mean_ranks(scores: Dict[ScoreKey, float]) -> Dict[str, Optional[float]]:
    algos = sorted({algo for _ref, _truth, algo in scores})
    cells = {(ref, truth) for ref, truth, _algo in scores}
    rank_lists: Dict[str, List[float]] = defaultdict(list)
    for ref, truth in cells:
        present = [algo for algo in algos if (ref, truth, algo) in scores]
        ranks = average_ranks([scores[(ref, truth, algo)] for algo in present])
        for algo, rank in zip(present, ranks):
            rank_lists[algo].append(rank)
    return {
        algo: (sum(rank_lists[algo]) / len(rank_lists[algo]) if rank_lists[algo] else None)
        for algo in algos
    }


def percent_ahead(scores_here: Dict[str, float]) -> Dict[str, Optional[float]]:
    if not scores_here:
        return {}
    best = max(scores_here.values())
    percent: Dict[str, Optional[float]] = {}
    for algo, score in scores_here.items():
        if score == best or score == 0.0:
            percent[algo] = None
        else:
            percent[algo] = (best - score) / score * 100.0
    return percent


def display_name(algo: str, styles: dict) -> str:
    if algo == "reseek3_sensitive":
        return "Reseek3"
    style = plot_curve.series_style(algo, 0, styles)
    if style.legend_label:
        return style.legend_label
    return algo


def bar_color(algo: str, styles: dict) -> str:
    return plot_curve.series_style(algo, 0, styles).color


def color_variants(color: str) -> Tuple[tuple, tuple, tuple]:
    """Light, medium, and dark shades of one paper color."""
    red, green, blue = mcolors.to_rgb(color)
    hue, _lightness, saturation = colorsys.rgb_to_hls(red, green, blue)
    shades = []
    for shade_l, sat_scale in zip(SHADE_LIGHTNESS, SHADE_SATURATION):
        shade_sat = 0.0 if saturation == 0.0 else min(1.0, saturation * sat_scale)
        shades.append(colorsys.hls_to_rgb(hue, shade_l, shade_sat))
    return shades[0], shades[1], shades[2]


def format_sig(value: float) -> str:
    """Two or three significant figures."""
    return f"{value:.3g}"


def draw_row(
    ax,
    tax,
    title: str,
    parts_here: Dict[str, Parts],
    percent: Dict[str, Optional[float]],
    mean_rank: Dict[str, Optional[float]],
    styles: dict,
    show_header: bool,
    show_xlabel: bool,
) -> None:
    """Stacked horizontal bars, largest total at the top, table aligned to the bars."""
    totals = {algo: sum(parts) for algo, parts in parts_here.items()}
    ordered = sorted(totals, key=lambda algo: (totals[algo], algo))
    n = len(ordered)
    ys = list(range(n))
    for y, algo in zip(ys, ordered):
        light, medium, dark = color_variants(bar_color(algo, styles))
        cve, pr, roc = parts_here[algo]
        ax.barh(y, cve, color=light, height=BAR_HEIGHT, linewidth=0)
        ax.barh(y, pr, left=cve, color=medium, height=BAR_HEIGHT, linewidth=0)
        ax.barh(y, roc, left=cve + pr, color=dark, height=BAR_HEIGHT, linewidth=0)
    ax.set_yticks(ys)
    ax.set_yticklabels([display_name(algo, styles) for algo in ordered])
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xlim(0, FRACTION_MAX)
    ax.set_xticks([0, 0.5, 1])
    if show_xlabel:
        ax.set_xlabel("Fraction of best", fontsize=FONT, labelpad=1)
    else:
        ax.tick_params(axis="x", labelbottom=False)
    ax.set_title(title, loc="left", fontsize=TITLE_FONT, pad=1)
    ax.grid(True, axis="x", linewidth=0.4, alpha=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=FONT, pad=1)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    tax.set_xlim(0, 1)
    tax.set_ylim(ax.get_ylim())
    tax.set_xticks([])
    tax.set_yticks([])
    for spine in tax.spines.values():
        spine.set_visible(False)
    trans = blended_transform_factory(tax.transAxes, tax.transData)
    if show_header:
        for x, label in TABLE_COLUMNS:
            tax.text(
                x,
                1.04,
                label,
                transform=tax.transAxes,
                ha="center",
                va="bottom",
                fontsize=FONT,
            )
    for y, algo in zip(ys, ordered):
        ahead = percent[algo]
        rank = mean_rank.get(algo)
        ahead_s = "" if ahead is None else format_sig(ahead)
        rank_s = "" if rank is None else format_sig(rank)
        values = (format_sig(totals[algo]), ahead_s, rank_s)
        for (x, _label), text in zip(TABLE_COLUMNS, values):
            tax.text(x, y, text, transform=trans, ha="center", va="center", fontsize=FONT)

    sys.stderr.write(f"{title}\n")
    for algo in reversed(ordered):
        ahead = percent[algo]
        rank = mean_rank.get(algo)
        ahead_s = "" if ahead is None else format_sig(ahead)
        rank_s = "" if rank is None else format_sig(rank)
        sys.stderr.write(
            f"  {algo}\t{format_sig(totals[algo])}\t{ahead_s}\t{rank_s}\n"
        )


def draw_panels(
    rows: Sequence[dict],
    parts: Dict[ScoreKey, Parts],
    styles: dict,
    axes_pairs: Sequence[Tuple],
) -> None:
    totals = as_totals(parts)
    ranks = mean_ranks(totals)
    cells = cell_list(rows, totals)
    last = len(cells) - 1
    for index, ((reference, truth), (bar_ax, table_ax)) in enumerate(zip(cells, axes_pairs)):
        here = {
            algo: value
            for (ref, tru, algo), value in parts.items()
            if ref == reference and tru == truth
        }
        scores = {algo: sum(value) for algo, value in here.items()}
        draw_row(
            bar_ax,
            table_ax,
            f"{reference} {truth}",
            here,
            percent_ahead(scores),
            ranks,
            styles,
            show_header=(index == 0),
            show_xlabel=(index == last),
        )


def panel_bar_counts(rows: Sequence[dict], parts: Dict[ScoreKey, Parts]) -> List[int]:
    totals = as_totals(parts)
    counts = []
    for reference, truth in cell_list(rows, totals):
        counts.append(
            sum(1 for ref, tru, _algo in parts if ref == reference and tru == truth)
        )
    return counts


def add_panel_axes(fig, counts: Sequence[int], fig_h: float) -> List[Tuple]:
    """Place each panel so its height is proportional to its bar count."""
    left = LABEL_IN / FIG_W
    bar_w = BAR_AXIS_IN / FIG_W
    table_left = (LABEL_IN + BAR_AXIS_IN + 0.06) / FIG_W
    table_w = TABLE_IN / FIG_W
    y = fig_h - TOP_IN
    pairs = []
    for index, n in enumerate(counts):
        height = INCH_PER_BAR * n
        y -= height
        bar_ax = fig.add_axes([left, y / fig_h, bar_w, height / fig_h])
        table_ax = fig.add_axes([table_left, y / fig_h, table_w, height / fig_h])
        pairs.append((bar_ax, table_ax))
        if index != len(counts) - 1:
            y -= PANEL_GAP_IN
    return pairs


def add_shade_legend(fig) -> None:
    light, medium, dark = color_variants("0.2")
    handles = [
        Patch(facecolor=light, edgecolor="none", label="CVE"),
        Patch(facecolor=medium, edgecolor="none", label="PR"),
        Patch(facecolor=dark, edgecolor="none", label="ROC"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        fontsize=FONT,
        bbox_to_anchor=((LABEL_IN + BAR_AXIS_IN / 2) / FIG_W, 0.995),
        borderpad=0.1,
        handlelength=1.2,
        handletextpad=0.4,
        columnspacing=1.0,
    )


def main() -> int:
    ap = argparse.ArgumentParser(
        description="One stacked AUC bar chart per reference and truth, with a summary table."
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
        default=REPO / "paper_plots" / "auc_summary.svg",
    )
    args = ap.parse_args()

    for path in (args.auc, args.auc_top):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing or empty table: {path}")

    styles = plot_curve.load_series_styles(str(args.colors))
    rows = rows_for_figure(read_tsv(args.auc), read_tsv(args.auc_top))
    parts = fraction_of_best(drop_zero_parts(cell_parts(rows)))
    parts = {
        key: tuple(value / 3.0 for value in triple)
        for key, triple in parts.items()
    }
    counts = panel_bar_counts(rows, parts)
    if len(counts) != len(SHOW):
        raise SystemExit(f"expected {len(SHOW)} charts, got {len(counts)}")

    n_panels = len(counts)
    fig_h = (
        TOP_IN
        + BOTTOM_IN
        + INCH_PER_BAR * sum(counts)
        + PANEL_GAP_IN * (n_panels - 1)
    )
    fig = plt.figure(figsize=(FIG_W, fig_h))
    axes_pairs = add_panel_axes(fig, counts, fig_h)
    draw_panels(rows, parts, styles, axes_pairs)
    add_shade_legend(fig)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    plt.close(fig)
    sys.stderr.write(f"wrote {args.output}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
