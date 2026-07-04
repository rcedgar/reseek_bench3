#!/usr/bin/env python3
"""
hits_to_topcat.py -- Stage 1b: condense hits to leave-category-out .tcat summary.

Single streaming pass.  For each query in the truth query set, tracks S_tp and
S_xf per truth_standards.md.  Top3 is computed from the implied CVE curve.

Example:
  python hits_to_topcat.py --hits algo.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json \
      --score --truth topsf --output algo.topsf.tcat
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Optional, Tuple

import common
import lookup


def is_tp_at_threshold(
    s_tp: Optional[float],
    s_xf: Optional[float],
    threshold: float,
    scores_are_evalues: bool,
) -> bool:
    """Positive-control TP at threshold T (truth_standards.md section 8)."""
    if not common.passes_threshold(s_tp, threshold, scores_are_evalues):
        return False
    if not common.passes_threshold(s_xf, threshold, scores_are_evalues):
        return True
    return common.at_least_as_good(s_tp, s_xf, scores_are_evalues)


def classify_at_threshold(
    scores: Dict[str, Tuple[Optional[float], Optional[float]]],
    threshold: float,
    n_q: int,
    scores_are_evalues: bool,
) -> Tuple[float, float]:
    """
    Return (coverage, error) at threshold T.

    coverage = TP / |Q|   (positive control, fixed denominator)
    error    = FP_neg / |Q|  (negative control FPR)
    """
    if n_q <= 0:
        return 0.0, 0.0

    n_tp = 0
    n_fp_neg = 0
    for s_tp, s_xf in scores.values():
        if is_tp_at_threshold(s_tp, s_xf, threshold, scores_are_evalues):
            n_tp += 1
        if common.passes_threshold(s_xf, threshold, scores_are_evalues):
            n_fp_neg += 1

    return n_tp / n_q, n_fp_neg / n_q


def topcat_cve_curve(
    scores: Dict[str, Tuple[Optional[float], Optional[float]]],
    n_q: int,
    scores_are_evalues: bool,
) -> List[Tuple[float, float]]:
    """Return (coverage, error) points sweeping thresholds strict-first."""
    taus = sort_thresholds_strict_first(
        collect_thresholds(scores), scores_are_evalues
    )
    taus.append(loose_threshold(scores_are_evalues))
    points: List[Tuple[float, float]] = []
    for tau in taus:
        cov, err = classify_at_threshold(scores, tau, n_q, scores_are_evalues)
        points.append((cov, err))
    return points


def collect_thresholds(
    scores: Dict[str, Tuple[Optional[float], Optional[float]]],
) -> List[float]:
    """Distinct non-missing S_tp and S_xf values across the query set."""
    seen: set[float] = set()
    for s_tp, s_xf in scores.values():
        if s_tp is not None:
            seen.add(s_tp)
        if s_xf is not None:
            seen.add(s_xf)
    return list(seen)


def sort_thresholds_strict_first(
    thresholds: List[float], scores_are_evalues: bool
) -> List[float]:
    return sorted(thresholds, reverse=not scores_are_evalues)


def loose_threshold(scores_are_evalues: bool) -> float:
    """Threshold accepting every reported score."""
    return float("inf") if scores_are_evalues else float("-inf")


def compute_top3(
    scores: Dict[str, Tuple[Optional[float], Optional[float]]],
    n_q: int,
    scores_are_evalues: bool,
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Compute TEPQ0.001, TEPQ0.01, TEPQ0.1, Top3.

    Sweep thresholds strict-first; record coverage at the first point where
    error >= each FPR threshold.  Unreached thresholds back-fill with final
    coverage (truth_standards.md section 9).
    """
    if n_q <= 0:
        return None, None, None, None

    taus = sort_thresholds_strict_first(
        collect_thresholds(scores), scores_are_evalues
    )
    taus.append(loose_threshold(scores_are_evalues))

    tepq0001: Optional[float] = None
    tepq001: Optional[float] = None
    tepq01: Optional[float] = None
    final_cov = 0.0

    for tau in taus:
        cov, err = classify_at_threshold(scores, tau, n_q, scores_are_evalues)
        final_cov = cov
        if err >= 0.001 and tepq0001 is None:
            tepq0001 = cov
        if err >= 0.01 and tepq001 is None:
            tepq001 = cov
        if err >= 0.1 and tepq01 is None:
            tepq01 = cov

    if tepq0001 is None:
        tepq0001 = final_cov
    if tepq001 is None:
        tepq001 = final_cov
    if tepq01 is None:
        tepq01 = final_cov

    top3 = tepq0001 * 2 + tepq001 * 1.5 + tepq01
    return tepq0001, tepq001, tepq01, top3


def condense_hits(
    hits_path: str,
    tax: lookup.Taxonomy,
    truth: str,
    query_set: frozenset[str],
    q_idx: int,
    t_idx: int,
    s_idx: int,
    scores_are_evalues: bool,
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """Stream hits; return domain -> (S_tp, S_xf) for queries in query_set."""
    scores: Dict[str, Tuple[Optional[float], Optional[float]]] = {
        dom: (None, None) for dom in sorted(query_set)
    }
    n_ignored = 0
    n_considered = 0

    with open(hits_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            flds = line.split()
            if len(flds) <= max(q_idx, t_idx, s_idx):
                continue
            q = lookup.parse_domain(flds[q_idx])
            t = lookup.parse_domain(flds[t_idx])
            try:
                score = float(flds[s_idx])
            except ValueError:
                continue

            if q not in query_set:
                continue

            track = lookup.topcat_track(q, t, truth, tax)
            if track is None:
                n_ignored += 1
                continue

            n_considered += 1
            s_tp, s_xf = scores[q]
            if track == "tp":
                if common.better(score, s_tp, scores_are_evalues):
                    s_tp = score
            else:
                if common.better(score, s_xf, scores_are_evalues):
                    s_xf = score
            scores[q] = (s_tp, s_xf)

    sys.stderr.write(f"ignored {n_ignored} hits, considered {n_considered}\n")
    return scores


def format_score(score: Optional[float]) -> str:
    return "." if score is None else f"{score:g}"


def write_tcat(
    out_path: str,
    scores: Dict[str, Tuple[Optional[float], Optional[float]]],
    *,
    algo: str,
    truth: str,
    fields_spec: str,
    scores_are_evalues: bool,
    n_possible_tp: int,
) -> None:
    """Write .tcat file with header and per-query S_tp / S_xf rows."""
    score_dir = "lower_better" if scores_are_evalues else "higher_better"
    tepq0001, tepq001, tepq01, top3 = compute_top3(
        scores, n_possible_tp, scores_are_evalues
    )

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(f"# algo={algo}\n")
        f.write(f"# truth={truth}\n")
        f.write(f"# fields={fields_spec}\n")
        f.write(f"# score_direction={score_dir}\n")
        f.write(f"# N_possible_tp={n_possible_tp}\n")
        if tepq0001 is not None:
            s = f"# TEPQ0.001={tepq0001:.3f} TEPQ0.01={tepq001:.3f} TEPQ0.1={tepq01:.3f} Top3={top3:.3f}\n"
            sys.stderr.write(s)
            f.write(s)
        f.write("domain\tS_tp\tS_xf\n")
        for dom in sorted(scores):
            s_tp, s_xf = scores[dom]
            f.write(f"{dom}\t{format_score(s_tp)}\t{format_score(s_xf)}\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Condense hits to homval leave-category-out .tcat summary."
    )
    ap.add_argument("--hits", required=True, help="Hits TSV path")
    ap.add_argument("--lookup", required=True, help="Lookup TSV path")
    ap.add_argument(
        "--derived-info",
        required=True,
        help="derived_info JSON from build_derived_info.py",
    )
    ap.add_argument(
        "--truth",
        required=True,
        choices=lookup.TOPCAT_TRUTHS,
        help="Leave-category-out truth standard (topfold or topsf)",
    )
    ap.add_argument("--output", required=True, help="Output .tcat path")
    ap.add_argument(
        "--fields",
        default="1,2,3",
        help="1-based column indices query,target,score (default 1,2,3)",
    )
    ap.add_argument(
        "--algo",
        default="",
        help="Algorithm label for .tcat header (default: hits basename)",
    )
    common.add_score_direction_args(ap)
    args = ap.parse_args()

    scores_are_evalues = common.scores_are_evalues_from_args(args)
    q_idx, t_idx, s_idx = common.parse_fields(args.fields)
    fields_spec = common.fields_spec_1based(q_idx, t_idx, s_idx)

    info = common.load_derived_info(args.derived_info)
    n_possible_tp = common.topcat_denominators(info, args.truth)
    query_set = common.topcat_query_set(info, args.truth)

    algo = args.algo or args.hits.replace("\\", "/").split("/")[-1]

    sys.stderr.write(f"\n\n")
    sys.stderr.write(f"loading lookup {args.lookup}...\n")
    tax = lookup.load_lookup(args.lookup)

    sys.stderr.write(
        f"reading hits {args.hits} truth={args.truth} "
        f"(|Q|={n_possible_tp})...\n"
    )
    scores = condense_hits(
        args.hits,
        tax,
        args.truth,
        query_set,
        q_idx,
        t_idx,
        s_idx,
        scores_are_evalues,
    )

    sys.stderr.write(f"writing {args.output}...\n")
    write_tcat(
        args.output,
        scores,
        algo=algo,
        truth=args.truth,
        fields_spec=fields_spec,
        scores_are_evalues=scores_are_evalues,
        n_possible_tp=n_possible_tp,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
