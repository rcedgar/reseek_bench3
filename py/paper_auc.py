#!/usr/bin/env python3
"""
paper_auc.py -- partial AUC for pairwise fam / sf / sfx curves.

References are scop40x and cath40. Summaries are
<edf-dir>/<algo>.<reference>.<family|superfamily|superfamilyx>.
Axis windows match paper_cve_pr_roc_scop40x_superfamily_narrow.svg.
Points are the polylines plot_curve.py draws. A vertical step (dx = 0)
contributes nothing. A missing summary (cath40 has no tm or dali) is
written as NA.

The integral does not extend past the first or last vertex, except CVE
area-over on the left. If the curve is on or below the bottom of the
window there, the shape is closed by the vertical line at xmin: the area
from xmin to the crossing of ymin is the rectangle
(x_intersect - xmin) * (ymax - ymin).

auc_linear
    Area inside the axes.
    CVE: trapezoidal ∫ (ymax - clip(y, ymin, ymax)) dx, area over the
    curve, so a lower error scores higher, plus that left rectangle.
    A curve above ymax contributes 0; a curve below ymin contributes
    ymax - ymin.
    PR and ROC: trapezoidal ∫ (clip(y, ymin, ymax) - ymin) dx, area under
    the curve, measured from the bottom of the plotted y range. A curve
    below ymin contributes 0; a curve above ymax contributes ymax - ymin.

auc_log
    CVE (log y): same area over the curve in log10(y), down from
    log10(ymax). The left rectangle has height log10(ymax) - log10(ymin).
    y <= 0 lies on that floor (log10 is undefined). Segments with positive
    y are straight in (coverage, log10(y)).
    ROC (log x): ∫ (clip(y, ymin, ymax) - ymin) d(log10(x)). Segments are
    straight in (log10(FPR), TPR). The x window must be positive.
    PR: na (both axes are linear).

Examples:
    python paper_auc.py
    python paper_auc.py --cve-xlim 0,0.2 --pr-xlim 0,0.2 --roc-xlim 1e-6,1e-5
    python paper_auc.py --cve-ylim 1e-2,10 --pr-ylim 0.7,1 --roc-ylim 0,0.4
    python paper_auc.py --edf-dir edf --algo foldseek --algo dali
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import common
import plot_curve

Point = Tuple[float, float]

# Algorithms and axis windows from paper_cve_pr_roc_scop40x_superfamily_narrow.bash.
DEFAULT_ALGOS = (
    "foldseek",
    "dali",
    "tm",
    "reseek3_fast",
    "reseek3_sensitive",
)
DEFAULT_REFERENCES = ("scop40x", "cath40")


class TruthSpec(NamedTuple):
    label: str
    filename: str


PAIRWISE_TRUTHS = (
    TruthSpec("fam", "family"),
    TruthSpec("sf", "superfamily"),
    TruthSpec("sfx", "superfamilyx"),
)

DEFINITIONS = """\
auc_linear CVE: area over the curve inside the axes, trapezoid of (ymax - clip(y, ymin, ymax)) vs x; left of the ymin crossing, rectangle (x_intersect - xmin) * (ymax - ymin)
auc_linear PR/ROC: area under the curve, trapezoid of (clip(y, ymin, ymax) - ymin) vs x
auc_log CVE: same area-over in log10(y), down from log10(ymax); left of the floor crossing, rectangle (x_intersect - xmin) * (log10(ymax) - log10(ymin)); y <= 0 lies on that floor
auc_log ROC: trapezoid of (clip(y, ymin, ymax) - ymin) vs log10(x)
auc_log PR: na (both axes linear)
"""


class Panel(NamedTuple):
    name: str
    xlim: Tuple[float, float]
    ylim: Tuple[float, float]
    log_axis: str  # "y", "x", or "none"


def parse_limit(spec: str, flag: str) -> Tuple[float, float]:
    try:
        lo_s, hi_s = spec.split(",", 1)
        lo, hi = float(lo_s), float(hi_s)
    except ValueError as exc:
        raise SystemExit(f"{flag}: expected LO,HI, got {spec!r}") from exc
    if not math.isfinite(lo) or not math.isfinite(hi) or not lo < hi:
        raise SystemExit(f"{flag}: need finite LO < HI, got {spec!r}")
    return lo, hi


def clip_polyline(points: Sequence[Point], xmin: float, xmax: float) -> List[Point]:
    """Portion of the polyline with x in [xmin, xmax], linear in (x, y)."""
    out: List[Point] = []

    def emit(x: float, y: float) -> None:
        if out and out[-1][0] == x and out[-1][1] == y:
            return
        out.append((x, y))

    if len(points) == 1:
        x, y = points[0]
        if xmin <= x <= xmax:
            emit(x, y)
        return out

    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 == x1:
            if xmin <= x0 <= xmax:
                emit(x0, y0)
                emit(x1, y1)
            continue

        dx = x1 - x0
        dy = y1 - y0
        t_lo = 0.0
        t_hi = 1.0
        for bound, is_min in ((xmin, True), (xmax, False)):
            t_bound = (bound - x0) / dx
            if dx > 0:
                if is_min:
                    t_lo = max(t_lo, t_bound)
                else:
                    t_hi = min(t_hi, t_bound)
            else:
                if is_min:
                    t_hi = min(t_hi, t_bound)
                else:
                    t_lo = max(t_lo, t_bound)
        if t_lo <= t_hi:
            emit(x0 + t_lo * dx, y0 + t_lo * dy)
            emit(x0 + t_hi * dx, y0 + t_hi * dy)
    return out


def segment_visible_area(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    ymin: float,
    ymax: float,
    from_ymax: bool = False,
) -> float:
    """Area of one segment inside ymin <= y <= ymax.

    Default height is clip(y) - ymin (area under). from_ymax uses
    ymax - clip(y) (area over). The segment is split where it crosses ymin
    or ymax so a linear crossing is exact. Vertical steps (dx = 0)
    contribute nothing.
    """
    if x0 == x1:
        return 0.0
    dy = y1 - y0
    cuts = [0.0, 1.0]
    if dy != 0.0:
        for bound in (ymin, ymax):
            t = (bound - y0) / dy
            if 0.0 < t < 1.0:
                cuts.append(t)
    cuts.sort()
    area = 0.0
    prev = cuts[0]
    dx = x1 - x0
    for t in cuts[1:]:
        if t <= prev:
            continue
        ya = y0 + prev * dy
        yb = y0 + t * dy
        clipped_a = min(max(ya, ymin), ymax)
        clipped_b = min(max(yb, ymin), ymax)
        if from_ymax:
            ha = ymax - clipped_a
            hb = ymax - clipped_b
        else:
            ha = clipped_a - ymin
            hb = clipped_b - ymin
        area += (t - prev) * dx * (ha + hb) / 2.0
        prev = t
    return area


def visible_auc(
    points: Sequence[Point],
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    from_ymax: bool = False,
) -> Tuple[float, Optional[Tuple[float, float]]]:
    """Area of the polyline inside the axis window.

    Height is measured up from ymin, or down from ymax when from_ymax is
    set. Returns (area, integrated x span). Span is None when the overlap
    has no width.
    """
    clipped = clip_polyline(points, xmin, xmax)
    if len(clipped) < 2:
        return 0.0, None
    span = (clipped[0][0], clipped[-1][0])
    if span[1] == span[0]:
        return 0.0, None
    area = 0.0
    for (x0, y0), (x1, y1) in zip(clipped, clipped[1:]):
        area += segment_visible_area(x0, y0, x1, y1, ymin, ymax, from_ymax)
    return area, span


def curve_x_span(points: Sequence[Point]) -> Optional[Tuple[float, float]]:
    if not points:
        return None
    xs = [x for x, _y in points]
    return min(xs), max(xs)


def covers_window(span: Optional[Tuple[float, float]], xmin: float, xmax: float) -> bool:
    if span is None:
        return False
    return span[0] <= xmin and span[1] >= xmax


def fmt(value: float) -> str:
    return f"{value:.8g}"


def fmt_span(span: Optional[Tuple[float, float]]) -> str:
    if span is None:
        return "none"
    return f"[{fmt(span[0])}, {fmt(span[1])}]"


def extend_floor_to_xmin(
    points: Sequence[Point], xmin: float, xmax: float, ymin: float
) -> List[Point]:
    """Close area-over to xmin when the curve is on or below ymin there.

    The added segment lies on ymin, so its area-over is the rectangle from
    xmin to the first in-window vertex. Vertices at or below ymin keep full
    height until the polyline crosses ymin, so the total left of that
    crossing is (x_intersect - xmin) * (ymax - ymin). A curve that enters
    the window from the left, or is already above ymin, is left unchanged.
    """
    if not points:
        return []
    first_inside: Optional[int] = None
    for i, (x, _y) in enumerate(points):
        if x >= xmin:
            first_inside = i
            break
    if first_inside is None or first_inside > 0:
        return list(points)
    x0, y0 = points[0]
    if x0 > xmax or x0 <= xmin or y0 > ymin:
        return list(points)
    return [(xmin, ymin), *points]


def warn_if_partial(
    panel: str,
    algo: str,
    which: str,
    points: Sequence[Point],
    xmin: float,
    xmax: float,
    integrated: Optional[Tuple[float, float]],
) -> None:
    if covers_window(integrated, xmin, xmax):
        return
    span = curve_x_span(points)
    if covers_window(span, xmin, xmax):
        return
    sys.stderr.write(
        f"warning: {panel} {algo} {which} curve x spans {fmt_span(span)}; "
        f"window is [{fmt(xmin)}, {fmt(xmax)}]; "
        f"integrated {fmt_span(integrated)} only\n"
    )


def cve_log_auc(
    points: Sequence[Point],
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
) -> float:
    """Area over the curve in (x, log10(y)), down from log10(ymax).

    y <= 0 lies on the floor log10(ymin); the left side is closed to xmin
    when the curve is on or below that floor. Positive-y segments are
    straight in (x, log10(y)).
    """
    log_ymin = math.log10(ymin)
    log_ymax = math.log10(ymax)
    logged = [
        (x, log_ymin if y <= 0.0 else math.log10(y)) for x, y in points
    ]
    logged = extend_floor_to_xmin(logged, xmin, xmax, log_ymin)
    area, _span = visible_auc(
        logged, xmin, xmax, log_ymin, log_ymax, from_ymax=True
    )
    return area


def roc_log_auc(
    points: Sequence[Point],
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    panel: str,
    algo: str,
) -> float:
    logged = [(math.log10(x), y) for x, y in points if x > 0.0]
    log_lo = math.log10(xmin)
    log_hi = math.log10(xmax)
    area, integrated = visible_auc(logged, log_lo, log_hi, ymin, ymax)
    # Report the warning in FPR, not log10(FPR).
    positive = [(x, y) for x, y in points if x > 0.0]
    pos_span = curve_x_span(positive)
    if not covers_window(pos_span, xmin, xmax):
        if integrated is None:
            shown = "none"
        else:
            shown = f"[{fmt(10 ** integrated[0])}, {fmt(10 ** integrated[1])}]"
        sys.stderr.write(
            f"warning: {panel} {algo} auc_log curve x spans {fmt_span(pos_span)}; "
            f"window is [{fmt(xmin)}, {fmt(xmax)}]; "
            f"integrated {shown} only\n"
        )
    return area


def auc_for_panel(
    panel: Panel, points: Sequence[Point], algo: str
) -> Tuple[float, str]:
    xmin, xmax = panel.xlim
    ymin, ymax = panel.ylim
    over = panel.name == "cve"
    linear_points = (
        extend_floor_to_xmin(points, xmin, xmax, ymin) if over else points
    )
    linear, integrated = visible_auc(
        linear_points, xmin, xmax, ymin, ymax, from_ymax=over
    )
    warn_if_partial(panel.name, algo, "auc_linear", points, xmin, xmax, integrated)
    if panel.log_axis == "y":
        log_auc = fmt(cve_log_auc(points, xmin, xmax, ymin, ymax))
    elif panel.log_axis == "x":
        log_auc = fmt(roc_log_auc(points, xmin, xmax, ymin, ymax, panel.name, algo))
    else:
        log_auc = "na"
    return linear, log_auc


class FigureSpec(NamedTuple):
    description: str
    dir_flag: str
    dir_help: str
    default_dir: str
    algos: Tuple[str, ...]
    references: Tuple[str, ...]
    truths: Tuple[TruthSpec, ...]
    cve_xlim: str
    pr_xlim: str
    roc_xlim: str
    cve_ylim: str
    pr_ylim: str
    roc_ylim: str


SUPERFAMILY = FigureSpec(
    description=(
        "Partial AUC for pairwise fam, sf, and sfx on scop40x and cath40. "
        "Axis windows match paper_cve_pr_roc_scop40x_superfamily_narrow.svg."
    ),
    dir_flag="edf-dir",
    dir_help=(
        "Directory of <algo>.<reference>.<truth> summaries (default: <repo>/edf)"
    ),
    default_dir="edf",
    algos=DEFAULT_ALGOS,
    references=DEFAULT_REFERENCES,
    truths=PAIRWISE_TRUTHS,
    cve_xlim="0,0.4",
    pr_xlim="0,0.4",
    roc_xlim="1e-7,1e-4",
    cve_ylim="1e-2,10",
    pr_ylim="0.7,1",
    roc_ylim="0,0.4",
)


def build_panels(args: argparse.Namespace) -> List[Panel]:
    cve_x = parse_limit(args.cve_xlim, "--cve-xlim")
    pr_x = parse_limit(args.pr_xlim, "--pr-xlim")
    roc_x = parse_limit(args.roc_xlim, "--roc-xlim")
    cve_y = parse_limit(args.cve_ylim, "--cve-ylim")
    pr_y = parse_limit(args.pr_ylim, "--pr-ylim")
    roc_y = parse_limit(args.roc_ylim, "--roc-ylim")
    if roc_x[0] <= 0.0 or roc_x[1] <= 0.0:
        raise SystemExit(
            f"--roc-xlim: log-x window must be positive, got {args.roc_xlim!r}"
        )
    if cve_y[0] <= 0.0:
        raise SystemExit(
            f"--cve-ylim: log-y window must be positive, got {args.cve_ylim!r}"
        )
    return [
        Panel("cve", cve_x, cve_y, "y"),
        Panel("pr", pr_x, pr_y, "none"),
        Panel("roc", roc_x, roc_y, "x"),
    ]


def run(spec: FigureSpec, argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=spec.description)
    ap.add_argument(
        "--algo",
        action="append",
        dest="algos",
        metavar="NAME",
        help=(
            "Algorithm to include (repeatable). "
            f"Default: {', '.join(spec.algos)}. "
            f"Reads <{spec.dir_flag}>/ALGO.REFERENCE.TRUTH"
        ),
    )
    ap.add_argument(
        f"--{spec.dir_flag}",
        dest="summary_dir",
        type=Path,
        default=None,
        help=spec.dir_help,
    )
    ap.add_argument(
        "--cve-xlim",
        default=spec.cve_xlim,
        help=f"CVE x window LO,HI (default {spec.cve_xlim})",
    )
    ap.add_argument(
        "--pr-xlim",
        default=spec.pr_xlim,
        help=f"PR x window LO,HI (default {spec.pr_xlim})",
    )
    ap.add_argument(
        "--roc-xlim",
        default=spec.roc_xlim,
        help=f"ROC x window LO,HI (default {spec.roc_xlim}); must be positive",
    )
    ap.add_argument(
        "--cve-ylim",
        default=spec.cve_ylim,
        help=f"CVE y window LO,HI (default {spec.cve_ylim}); must be positive",
    )
    ap.add_argument(
        "--pr-ylim",
        default=spec.pr_ylim,
        help=f"PR y window LO,HI (default {spec.pr_ylim})",
    )
    ap.add_argument(
        "--roc-ylim",
        default=spec.roc_ylim,
        help=f"ROC y window LO,HI (default {spec.roc_ylim})",
    )
    args = ap.parse_args(argv)

    algos = list(args.algos) if args.algos else list(spec.algos)
    if args.summary_dir is not None:
        summary_dir = args.summary_dir
    else:
        summary_dir = Path(__file__).resolve().parent.parent / spec.default_dir
    panels = build_panels(args)

    sys.stderr.write(DEFINITIONS)

    rows: List[str] = []
    for reference in spec.references:
        for truth in spec.truths:
            for algo_name in algos:
                path = summary_dir / f"{algo_name}.{reference}.{truth.filename}"
                curves: Dict[str, Tuple[str, List[Point]]] = {}
                if path.is_file():
                    for panel in panels:
                        hdr, _kind, points = plot_curve.curve_points_for_input(
                            str(path), panel.name
                        )
                        algo = common.algo_label_from_header(hdr, algo_name)
                        curves[panel.name] = (algo, points)
                        sys.stderr.write(
                            f"loaded {path} ({len(points)} points, "
                            f"algo={algo}, panel={panel.name})\n"
                        )
                else:
                    sys.stderr.write(f"warning: missing {path}; writing NA\n")

                for panel in panels:
                    xmin, xmax = panel.xlim
                    ymin, ymax = panel.ylim
                    if panel.name in curves:
                        algo, points = curves[panel.name]
                        linear, log_auc = auc_for_panel(panel, points, algo)
                        linear_s = fmt(linear)
                    else:
                        algo = algo_name
                        linear_s = "NA"
                        log_auc = "NA"
                    rows.append(
                        "\t".join(
                            (
                                panel.name,
                                reference,
                                truth.label,
                                algo,
                                fmt(xmin),
                                fmt(xmax),
                                fmt(ymin),
                                fmt(ymax),
                                linear_s,
                                log_auc,
                            )
                        )
                    )

    header = "\t".join(
        (
            "panel",
            "reference",
            "truth",
            "algo",
            "xmin",
            "xmax",
            "ymin",
            "ymax",
            "auc_linear",
            "auc_log",
        )
    )
    sys.stdout.write(header + "\n")
    for row in rows:
        sys.stdout.write(row + "\n")
    return 0


def main() -> int:
    return run(SUPERFAMILY)


if __name__ == "__main__":
    raise SystemExit(main())
