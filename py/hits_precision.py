#!/usr/bin/env python3
"""
hits_precision.py -- Score-to-precision lookup table from a hits file.

Columns 1 and 2 are query and target, each labeled domain/family
(class.fold.superfamily.family). Self-hits are ignored. A hit is a true
positive when query and target share --truth: fam, sf, or fold.

Cumulative precision at a score is cum_tp / (cum_tp + cum_fp). The table has
nr_bins knots spaced evenly in precision from the best score's precision down
to min_precision, or to the last score's precision if the curve never falls
that far. Each knot's score is the first downward crossing of that precision,
interpolated in log10(score) when --evalue is set and in raw score otherwise.
Scores better than the first knot clamp to the first precision; scores worse
than the last knot are below the table.

The score column (--field, default 3) is higher-better unless --evalue.

Example:
  python hits_precision.py \\
      --input ../big_reseek3_search_scop40x_raw/sf.hits --truth sf
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


TRUTHS = ("fam", "sf", "fold")
Point = Tuple[float, float]


def split_label(label: str) -> Tuple[str, str]:
    """Return (domain, family_id) from 'domain/family_id'."""
    dom, sep, fam = label.partition("/")
    if not sep or not fam or not dom:
        raise ValueError(f"expected domain/family, got {label!r}")
    return dom, fam


def category(fam: str, truth: str) -> str:
    """Fold, superfamily, or family id from a 4-field family id."""
    flds = fam.split(".")
    if len(flds) < 4:
        raise ValueError(f"family_id must have 4 dot fields, got {fam!r}")
    if truth == "fold":
        return flds[0] + "." + flds[1]
    if truth == "sf":
        return flds[0] + "." + flds[1] + "." + flds[2]
    if truth == "fam":
        return fam
    raise ValueError(f"unknown truth {truth!r}")


def precision_curve(
    hits_path: str,
    truth: str,
    score_field: int,
    scores_are_evalues: bool,
) -> List[Point]:
    """
    Cumulative precision at each distinct score, best-first.

    Self-hits (same domain) are ignored. Precision is cum_tp / (cum_tp + cum_fp).
    """
    if score_field < 1:
        raise ValueError(f"field must be >= 1, got {score_field}")
    s_idx = score_field - 1
    q_idx, t_idx = 0, 1
    need = max(q_idx, t_idx, s_idx)

    hist: Dict[float, List[int]] = defaultdict(lambda: [0, 0])
    n_self = 0
    n_bad = 0
    n_considered = 0

    with open(hits_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            flds = line.split()
            if len(flds) <= need:
                n_bad += 1
                continue
            try:
                q_dom, q_fam = split_label(flds[q_idx])
                t_dom, t_fam = split_label(flds[t_idx])
                score = float(flds[s_idx])
                q_cat = category(q_fam, truth)
                t_cat = category(t_fam, truth)
            except ValueError:
                n_bad += 1
                continue

            if q_dom == t_dom:
                n_self += 1
                continue

            n_considered += 1
            if q_cat == t_cat:
                hist[score][0] += 1
            else:
                hist[score][1] += 1

    sys.stderr.write(
        f"ignored {n_self} self-hits, skipped {n_bad} bad lines, "
        f"considered {n_considered}\n"
    )
    if n_considered == 0:
        raise ValueError(f"{hits_path}: no hits to score")

    ordered = sorted(hist, reverse=not scores_are_evalues)
    cum_tp = 0
    cum_fp = 0
    points: List[Point] = []
    for score in ordered:
        n_tp, n_fp = hist[score]
        cum_tp += n_tp
        cum_fp += n_fp
        denom = cum_tp + cum_fp
        if denom == 0:
            continue
        points.append((score, cum_tp / denom))
    if not points:
        raise ValueError(f"{hits_path}: no hits to score")
    return points


def points_through_threshold(
    points: List[Point], min_precision: float
) -> Tuple[List[Point], float, float]:
    """
    Keep rows through the first precision below min_precision.

    Return (points, p_hi, p_bottom). p_hi is the precision at the best score.
    p_bottom is min_precision when the curve crosses it, otherwise the precision
    of the last row.
    """
    if not 0.0 <= min_precision <= 1.0:
        raise ValueError(f"min_precision must be in [0, 1], got {min_precision}")

    p_hi = points[0][1]
    if p_hi < min_precision:
        raise ValueError(
            f"best precision {p_hi:.6g} is below min_precision {min_precision:.6g}"
        )

    kept: List[Point] = []
    crossed = False
    for score, prec in points:
        kept.append((score, prec))
        if prec < min_precision:
            crossed = True
            break

    p_bottom = min_precision if crossed else kept[-1][1]
    return kept, p_hi, p_bottom


def precision_levels(p_hi: float, p_bottom: float, nr_bins: int) -> List[float]:
    """Even steps from p_hi down to p_bottom, including both ends."""
    if nr_bins < 1:
        raise ValueError(f"nr_bins must be >= 1, got {nr_bins}")
    if nr_bins == 1:
        return [p_hi]
    levels = [
        p_hi + (p_bottom - p_hi) * i / (nr_bins - 1)
        for i in range(nr_bins)
    ]
    levels[0] = p_hi
    levels[-1] = p_bottom
    return levels


def interpolate_score(
    left: Point, right: Point, target: float, use_log: bool
) -> float:
    """Score where precision crosses target, between left (above) and right (below)."""
    s0, p0 = left
    s1, p1 = right
    if target == p0 or p0 == p1:
        return s0
    if target == p1:
        return s1
    t = (p0 - target) / (p0 - p1)
    if use_log:
        if s0 <= 0.0 or s1 <= 0.0:
            raise ValueError("log10 interpolation requires scores > 0")
        log_s = math.log10(s0) + t * (math.log10(s1) - math.log10(s0))
        return 10.0 ** log_s
    return s0 + t * (s1 - s0)


def score_at_precision(points: List[Point], target: float, use_log: bool) -> float:
    """
    Score at the first downward crossing of target.

    Later recoveries above target are ignored. If precision never falls below
    target, return the last score.
    """
    last_ok: int | None = None
    for i, (score, prec) in enumerate(points):
        if prec >= target:
            last_ok = i
            continue
        if last_ok is None:
            raise ValueError(
                f"precision drops below {target:.6g} before any point reaches it"
            )
        return interpolate_score(points[last_ok], (score, prec), target, use_log)
    if last_ok is None:
        raise ValueError(f"no point with precision >= {target:.6g}")
    return points[last_ok][0]


def build_knots(
    points: List[Point],
    min_precision: float,
    nr_bins: int,
    use_log: bool,
) -> List[Point]:
    kept, p_hi, p_bottom = points_through_threshold(points, min_precision)
    knots: List[Point] = []
    for level in precision_levels(p_hi, p_bottom, nr_bins):
        knots.append((score_at_precision(kept, level, use_log), level))
    return knots


def format_float(value: float) -> str:
    """Shortest general-format decimal that round-trips through float()."""
    for digits in range(1, 18):
        text = format(value, f".{digits}g")
        if float(text) == value:
            return text
    return repr(value)


def write_table(
    path: str,
    hdr: Dict[str, str],
    knots: List[Point],
    min_precision: float,
    nr_bins: int,
    use_log: bool,
) -> None:
    interp = "log10_score" if use_log else "linear_score"
    lines: List[str] = []
    if "truth" in hdr:
        lines.append(f"# truth={hdr['truth']}")
    if "score_direction" in hdr:
        lines.append(f"# score_direction={hdr['score_direction']}")
    lines.append(f"# min_precision={format_float(min_precision)}")
    lines.append(f"# nr_bins={nr_bins}")
    lines.append(f"# interp={interp}")
    lines.append("# score\tprecision")
    for score, prec in knots:
        lines.append(f"{format_float(score)}\t{format_float(prec)}")
    text = "\n".join(lines) + "\n"
    if path == "-":
        sys.stdout.write(text)
        return
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build a score-to-precision interpolation table from a hits file."
    )
    ap.add_argument("--input", required=True, help="Hits file (query, target, score)")
    ap.add_argument(
        "--truth",
        required=True,
        choices=TRUTHS,
        help="True positive when query and target share this level",
    )
    ap.add_argument(
        "--field",
        type=int,
        default=3,
        help="1-based score column (default: 3). Query and target are columns 1 and 2",
    )
    ap.add_argument(
        "--output",
        default=None,
        help="Output TSV (default: <input>.precision.tsv)",
    )
    ap.add_argument(
        "--min-precision",
        type=float,
        default=0.01,
        help="Lowest precision to tabulate; keep precision >= this (default: 0.01)",
    )
    ap.add_argument(
        "--nr-bins",
        type=int,
        default=32,
        help="Number of output rows, including both ends (default: 32)",
    )
    ap.add_argument(
        "--evalue",
        action="store_true",
        help="Lower score is better (default: higher is better, as for raw scores)",
    )
    args = ap.parse_args()

    output = args.output if args.output is not None else args.input + ".precision.tsv"
    use_log = bool(args.evalue)

    try:
        points = precision_curve(args.input, args.truth, args.field, use_log)
        knots = build_knots(points, args.min_precision, args.nr_bins, use_log)
        hdr = {
            "truth": args.truth,
            "score_direction": "lower_better" if use_log else "higher_better",
        }
        write_table(output, hdr, knots, args.min_precision, args.nr_bins, use_log)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
