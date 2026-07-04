#!/usr/bin/env python3
"""
hits_to_edf.py -- Stage 1a: condense pair-list hits to an .edf summary file.

Single streaming pass over a hits TSV.  Non-ignored hits are classified TP/FP
per the truth standard; counts are accumulated per distinct score value.

CVE coverage = cum_tp / N_possible_tp
CVE error    = cum_fp / ndom   (NOT ROC FPR)

Sum3 = 2*SEPQ0.1 + 1.5*SEPQ1 + SEPQ10  (coverage at error <= 0.1, 1, 10)
SFFP = fraction of all possible TPs with score strictly better than first FP

Example:
  python hits_to_edf.py --hits algo.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json \
      --truth superfamily --score --output algo.superfamily.edf
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import common
import lookup


def compute_sum3(
    rows: List[Tuple[float, int, int]],
    ndom: int,
    n_possible_tp: int,
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Compute SEPQ0.1, SEPQ1, SEPQ10, Sum3 from EDF rows (score, cum_tp, cum_fp).

    rows must be sorted best-first.  Matches FastBench::BenchPairCVE: record
    coverage at the last point before each score change where error <= threshold.
    Unreached thresholds back-fill with final coverage.
    """
    if n_possible_tp <= 0:
        return None, None, None, None

    sepq01: Optional[float] = None
    sepq1: Optional[float] = None
    sepq10: Optional[float] = None
    last_score: Optional[float] = None

    for score, cum_tp, cum_fp in rows:
        if last_score is not None and score != last_score:
            cov = cum_tp / n_possible_tp
            err = cum_fp / ndom
            if err <= 0.1:
                sepq01 = cov
            if err <= 1:
                sepq1 = cov
            if err <= 10:
                sepq10 = cov
        last_score = score

    if not rows:
        return None, None, None, None

    final_cov = rows[-1][1] / n_possible_tp
    if sepq01 is None:
        sepq01 = final_cov
    if sepq1 is None:
        sepq1 = final_cov
    if sepq10 is None:
        sepq10 = final_cov

    sum3 = sepq01 * 2 + sepq1 * 1.5 + sepq10
    return sepq01, sepq1, sepq10, sum3


def compute_sffp(
    per_query_hits: Dict[str, List[Tuple[float, bool]]],
    n_possible_tp: int,
    scores_are_evalues: bool,
) -> Optional[float]:
    """
    Sensitivity above First False Positive (per-query).

    For each query, find its first FP (best-scoring FP), then count TPs for
    that query with score strictly better than that FP.  Sum across queries
    and divide by N_possible_tp.

    per_query_hits: query -> list of (score, is_tp).
    If a query has no FPs, all its TPs count.
    """
    if n_possible_tp <= 0 or not per_query_hits:
        return None

    total_tp_above = 0
    for _query, hits in per_query_hits.items():
        tp_scores: List[float] = []
        best_fp_score: Optional[float] = None
        for score, is_tp in hits:
            if is_tp:
                tp_scores.append(score)
            else:
                if best_fp_score is None or common.better(score, best_fp_score, scores_are_evalues):
                    best_fp_score = score

        if best_fp_score is None:
            total_tp_above += len(tp_scores)
        else:
            for s in tp_scores:
                if common.better(s, best_fp_score, scores_are_evalues):
                    total_tp_above += 1

    return total_tp_above / n_possible_tp


def sort_scores_best_first(
    scores: List[float], scores_are_evalues: bool
) -> List[float]:
    return sorted(scores, reverse=not scores_are_evalues)


def condense_hits(
    hits_path: str,
    tax: lookup.Taxonomy,
    truth: str,
    q_idx: int,
    t_idx: int,
    s_idx: int,
    scores_are_evalues: bool,
) -> Tuple[Dict[float, List[int]], Dict[str, List[Tuple[float, bool]]]]:
    """Stream hits file; return (histogram, per_query_hits).

    histogram: score -> [n_tp, n_fp]
    per_query_hits: query -> list of (score, is_tp) for SFFP calculation
    """
    hist: Dict[float, List[int]] = defaultdict(lambda: [0, 0])
    per_query: Dict[str, List[Tuple[float, bool]]] = defaultdict(list)
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

            if lookup.pair_ignore(q, t, truth, tax):
                n_ignored += 1
                continue

            n_considered += 1
            is_tp = lookup.pair_is_tp(q, t, truth, tax)
            if is_tp:
                hist[score][0] += 1
            else:
                hist[score][1] += 1
            per_query[q].append((score, is_tp))

    sys.stderr.write(f"ignored {n_ignored} hits, considered {n_considered}\n")
    return hist, per_query


def write_edf(
    out_path: str,
    hist: Dict[float, List[int]],
    per_query_hits: Dict[str, List[Tuple[float, bool]]],
    *,
    algo: str,
    reference: str,
    truth: str,
    fields_spec: str,
    scores_are_evalues: bool,
    ndom: int,
    n_possible_tp: int,
    n_possible_fp: int,
) -> None:
    """Write .edf file with header, cumulative columns, and Sum3 stats."""
    score_dir = "lower_better" if scores_are_evalues else "higher_better"
    ordered = sort_scores_best_first(list(hist.keys()), scores_are_evalues)

    cum_tp = 0
    cum_fp = 0
    body_rows: List[str] = []
    curve_rows: List[Tuple[float, int, int]] = []

    for score in ordered:
        n_tp, n_fp = hist[score]
        cum_tp += n_tp
        cum_fp += n_fp
        curve_rows.append((score, cum_tp, cum_fp))
        body_rows.append(f"{score:g}\t{n_tp}\t{n_fp}\t{cum_tp}\t{cum_fp}")

    sepq01, sepq1, sepq10, sum3 = compute_sum3(curve_rows, ndom, n_possible_tp)
    sffp = compute_sffp(per_query_hits, n_possible_tp, scores_are_evalues)

    with open(out_path, "w", encoding="utf-8") as f:
        common.write_algo_reference_header(f, algo, reference)
        f.write(f"# truth={truth}\n")
        f.write(f"# fields={fields_spec}\n")
        f.write(f"# score_direction={score_dir}\n")
        f.write(f"# ndom={ndom}\n")
        f.write(f"# N_possible_tp={n_possible_tp}\n")
        f.write(f"# N_possible_fp={n_possible_fp}\n")
        if sepq01 is not None:
            s = f"# SEPQ0.1={sepq01:.3f} SEPQ1={sepq1:.3f} SEPQ10={sepq10:.3f} Sum3={sum3:.3f} SFFP={sffp:.3f}\n"
            sys.stderr.write(s)
            f.write(s)
        f.write("score\tn_tp\tn_fp\tcum_tp\tcum_fp\n")
        for line in body_rows:
            f.write(line + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Condense pair-list hits to homval .edf summary."
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
        choices=lookup.TRUTH_STANDARDS,
        help="Truth standard for TP/FP classification",
    )
    ap.add_argument("--output", required=True, help="Output .edf path")
    ap.add_argument(
        "--fields",
        default="1,2,3",
        help="1-based column indices query,target,score (default 1,2,3)",
    )
    ap.add_argument(
        "--algo",
        default="",
        help="Algorithm label for .edf header (default: from hits basename)",
    )
    ap.add_argument(
        "--reference",
        default=None,
        help="Reference DB label for .edf header (default: from hits basename)",
    )
    common.add_score_direction_args(ap)
    args = ap.parse_args()

    scores_are_evalues = common.scores_are_evalues_from_args(args)
    q_idx, t_idx, s_idx = common.parse_fields(args.fields)
    fields_spec = common.fields_spec_1based(q_idx, t_idx, s_idx)

    info = common.load_derived_info(args.derived_info)
    ndom, n_possible_tp, n_possible_fp = common.pair_denominators(info, args.truth)

    algo, reference = common.resolve_algo_reference(
        args.hits, algo=args.algo, reference=args.reference
    )

    sys.stderr.write(f"\n\n")
    sys.stderr.write(f"loading lookup {args.lookup}...\n")
    tax = lookup.load_lookup(args.lookup)

    sys.stderr.write(f"reading hits {args.hits} truth={args.truth}...\n")
    hist, per_query_hits = condense_hits(
        args.hits, tax, args.truth, q_idx, t_idx, s_idx, scores_are_evalues
    )

    sys.stderr.write(f"writing {args.output} ({len(hist)} distinct scores)...\n")
    write_edf(
        args.output,
        hist,
        per_query_hits,
        algo=algo,
        reference=reference,
        truth=args.truth,
        fields_spec=fields_spec,
        scores_are_evalues=scores_are_evalues,
        ndom=ndom,
        n_possible_tp=n_possible_tp,
        n_possible_fp=n_possible_fp,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
