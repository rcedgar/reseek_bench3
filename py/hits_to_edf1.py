#!/usr/bin/env python3
"""
hits_to_edf1.py -- per-query sensitivity quantiles for plotting.

Same pair-list truth standards as hits_to_edf.py.  For each lookup domain
with at least one possible true partner, sensitivity is the fraction of those
partners reported with a score strictly better than that query's first false
positive.  Self-hits are excluded.  Queries are sorted best-first and sampled
at --bins evenly spaced fractions of the query set (default 16).

Output columns:
  qfract   fraction of queries (k / bins)
  sens     sensitivity of the query at that quantile

Example:
  python hits_to_edf1.py --hits algo.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json \
      --truth superfamily --score --output algo.superfamily.edf1
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import Dict, List, Optional, Tuple

import common
import lookup


def possible_tp_count(q: str, truth: str, tax: lookup.Taxonomy) -> int:
    """Number of true-positive targets for query q, excluding the self-hit."""
    if truth == "fold":
        return len(tax.fold2doms[tax.dom2fold[q]]) - 1
    if truth in ("superfamily", "superfamilyx"):
        return len(tax.sf2doms[tax.dom2sf[q]]) - 1
    if truth == "family":
        return len(tax.fam2doms[tax.dom2fam[q]]) - 1
    raise ValueError(f"unknown truth standard: {truth!r}")


def collect_query_hits(
    hits_path: str,
    tax: lookup.Taxonomy,
    truth: str,
    q_idx: int,
    t_idx: int,
    s_idx: int,
    scores_are_evalues: bool,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, float]]:
    """Stream hits. Return (best TP score per query/target, best FP score per query)."""
    tp_best: Dict[str, Dict[str, float]] = {}
    best_fp: Dict[str, float] = {}
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
            if lookup.pair_is_tp(q, t, truth, tax):
                scores = tp_best.get(q)
                if scores is None:
                    scores = {}
                    tp_best[q] = scores
                prev = scores.get(t)
                if prev is None or common.better(score, prev, scores_are_evalues):
                    scores[t] = score
            else:
                prev_fp = best_fp.get(q)
                if prev_fp is None or common.better(score, prev_fp, scores_are_evalues):
                    best_fp[q] = score

    sys.stderr.write(f"ignored {n_ignored} hits, considered {n_considered}\n")
    return tp_best, best_fp


def query_sensitivities(
    tax: lookup.Taxonomy,
    truth: str,
    tp_best: Dict[str, Dict[str, float]],
    best_fp: Dict[str, float],
    scores_are_evalues: bool,
    n_possible_tp: int,
) -> List[float]:
    """
    Sensitivity of each eligible query, best-first.

    Eligible queries have at least one possible TP.  A missing hit list scores 0.
    The sum of per-query denominators must match derived_info N_possible_tp.
    """
    sensitivities: List[float] = []
    denom_sum = 0

    for q in tax.dom2fam:
        n_pos = possible_tp_count(q, truth, tax)
        denom_sum += n_pos
        if n_pos <= 0:
            continue

        fp: Optional[float] = best_fp.get(q)
        n_above = 0
        for score in tp_best.get(q, {}).values():
            if fp is None or common.better(score, fp, scores_are_evalues):
                n_above += 1
        sensitivities.append(n_above / n_pos)

    if denom_sum != n_possible_tp:
        raise ValueError(
            f"per-query TP denominators sum to {denom_sum}, "
            f"derived_info N_possible_tp={n_possible_tp}"
        )

    sensitivities.sort(reverse=True)
    return sensitivities


def quantile_rows(sensitivities: List[float], bins: int) -> List[Tuple[float, float]]:
    """Sample the best-first sensitivity list at qfract = k/bins."""
    if bins < 1:
        raise ValueError(f"bins must be >= 1, got {bins}")
    n = len(sensitivities)
    if n == 0:
        return []

    rows: List[Tuple[float, float]] = []
    for k in range(1, bins + 1):
        qfract = k / bins
        rank = (k * n + bins - 1) // bins  # ceil(k * n / bins)
        idx = rank - 1
        rows.append((qfract, sensitivities[idx]))
    return rows


def write_edf1(
    out_path: str,
    rows: List[Tuple[float, float]],
    *,
    algo: str,
    reference: str,
    truth: str,
    fields_spec: str,
    scores_are_evalues: bool,
    ndom: int,
    n_queries: int,
    n_possible_tp: int,
    bins: int,
) -> None:
    """Write the qfract/sens table with a comment header."""
    score_dir = "lower_better" if scores_are_evalues else "higher_better"
    with open(out_path, "w", encoding="utf-8") as f:
        common.write_algo_reference_header(f, algo, reference)
        f.write(f"# truth={truth}\n")
        f.write(f"# fields={fields_spec}\n")
        f.write(f"# score_direction={score_dir}\n")
        f.write(f"# ndom={ndom}\n")
        f.write(f"# N_queries={n_queries}\n")
        f.write(f"# N_possible_tp={n_possible_tp}\n")
        f.write(f"# bins={bins}\n")
        f.write("qfract\tsens\n")
        for qfract, sens in rows:
            f.write(f"{qfract:g}\t{sens:g}\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Per-query sensitivity quantiles from pair-list hits."
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
    ap.add_argument("--output", required=True, help="Output table path")
    ap.add_argument(
        "--fields",
        default="1,2,3",
        help="1-based column indices query,target,score (default 1,2,3)",
    )
    ap.add_argument(
        "--algo",
        default="",
        help="Algorithm label for the header (default: from hits basename)",
    )
    ap.add_argument(
        "--reference",
        default=None,
        help="Reference DB label for the header (default: from hits basename)",
    )
    ap.add_argument(
        "--bins",
        type=int,
        default=16,
        help="Number of query-fraction quantiles (default 16)",
    )
    common.add_score_direction_args(ap)
    args = ap.parse_args()

    if args.bins < 1:
        ap.error(f"--bins must be >= 1, got {args.bins}")

    scores_are_evalues = common.scores_are_evalues_from_args(args)
    q_idx, t_idx, s_idx = common.parse_fields(args.fields)
    fields_spec = common.fields_spec_1based(q_idx, t_idx, s_idx)

    info = common.load_derived_info(args.derived_info)
    ndom, n_possible_tp, _n_possible_fp = common.pair_denominators(info, args.truth)

    algo, reference = common.resolve_algo_reference(
        args.hits, algo=args.algo, reference=args.reference
    )

    sys.stderr.write("\n\n")
    sys.stderr.write(f"loading lookup {args.lookup}...\n")
    tax = lookup.load_lookup(args.lookup)

    sys.stderr.write(f"reading hits {args.hits} truth={args.truth}...\n")
    tp_best, best_fp = collect_query_hits(
        args.hits, tax, args.truth, q_idx, t_idx, s_idx, scores_are_evalues
    )

    sensitivities = query_sensitivities(
        tax, args.truth, tp_best, best_fp, scores_are_evalues, n_possible_tp
    )
    if not sensitivities:
        sys.stderr.write("no queries with a possible TP\n")
    rows = quantile_rows(sensitivities, args.bins)

    sys.stderr.write(
        f"writing {args.output} ({len(rows)} quantiles, "
        f"{len(sensitivities)} queries)...\n"
    )
    write_edf1(
        args.output,
        rows,
        algo=algo,
        reference=reference,
        truth=args.truth,
        fields_spec=fields_spec,
        scores_are_evalues=scores_are_evalues,
        ndom=ndom,
        n_queries=len(sensitivities),
        n_possible_tp=n_possible_tp,
        bins=args.bins,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
