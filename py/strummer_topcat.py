#!/usr/bin/env python3
"""
strummer_topcat.py -- topsf summary from Strummer top hits on scop40t.

Queries are scop40q (test.lookup). The database is scop40t (train.lookup).
A query is included only when scop40t has a different family of the same
superfamily, so a correct superfamily call is reachable.

The merged hit file is already the top Strummer call per query:
  1 domain/family
  2 predicted superfamily id
  3 P-value (lower is better)

Ties at that P-value are all kept. The predicted superfamily is column 2;
other hits are not in this file.

Example:
  python strummer_topcat.py \\
      --hits /mnt/c/src/strummer_scop/merged_hits/merged \\
      --queries /mnt/c/src/strummer_scop/lookup/test.lookup \\
      --targets /mnt/c/src/strummer_scop/lookup/train.lookup \\
      --output ../topcat/strummer.scop40t.topsf
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Optional, Tuple

import common
import hits_to_topcat
import lookup


def reachable_query_set(
    q_tax: lookup.Taxonomy, t_tax: lookup.Taxonomy
) -> frozenset[str]:
    """scop40q domains with a different family of the same superfamily in scop40t."""
    sf_fams: Dict[str, set[str]] = {}
    for sf, doms in t_tax.sf2doms.items():
        sf_fams[sf] = {t_tax.dom2fam[dom] for dom in doms}
    out: List[str] = []
    for dom, fam in q_tax.dom2fam.items():
        fams = sf_fams.get(q_tax.dom2sf[dom])
        if fams is not None and any(other != fam for other in fams):
            out.append(dom)
    return frozenset(out)


def note_top(
    best: Dict[str, List],
    q: str,
    score: float,
    same: bool,
    scores_are_evalues: bool,
) -> None:
    """Record a candidate top hit. rec is [score, saw_tp, saw_xf]."""
    rec = best.get(q)
    if rec is None or common.better(score, rec[0], scores_are_evalues):
        best[q] = [score, same, not same]
        return
    if score == rec[0]:
        if same:
            rec[1] = True
        else:
            rec[2] = True


def scores_from_best(
    query_set: frozenset[str],
    best: Dict[str, List],
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """Map each query to (S_tp, S_xf) from its top hit, or missing."""
    scores: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
    for dom in query_set:
        rec = best.get(dom)
        if rec is None:
            scores[dom] = (None, None)
            continue
        score, is_tp, is_xf = rec
        scores[dom] = (
            score if is_tp else None,
            score if is_xf else None,
        )
    return scores


def condense_strummer(
    hits_path: str,
    query_set: frozenset[str],
    q_tax: lookup.Taxonomy,
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """Top Strummer superfamily call per reachable query. P-values, lower better."""
    best: Dict[str, List] = {}
    q_sf = q_tax.dom2sf
    n_rows = 0
    n_kept = 0
    with open(hits_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            flds = line.split()
            if len(flds) < 3:
                continue
            n_rows += 1
            q = lookup.parse_domain(flds[0])
            if q not in query_set:
                continue
            try:
                score = float(flds[2])
            except ValueError:
                continue
            note_top(best, q, score, flds[1] == q_sf[q], True)
            n_kept += 1
    sys.stderr.write(
        f"strummer rows {n_rows}, kept {n_kept}, queries with a hit {len(best)}\n"
    )
    return scores_from_best(query_set, best)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Strummer top hits to a topsf .tcat summary on scop40t."
    )
    ap.add_argument("--hits", required=True, help="Merged Strummer hits TSV")
    ap.add_argument("--queries", required=True, help="scop40q lookup (test.lookup)")
    ap.add_argument("--targets", required=True, help="scop40t lookup (train.lookup)")
    ap.add_argument("--output", required=True, help="Output .topsf path")
    ap.add_argument("--algo", default="strummer", help="Algorithm label")
    args = ap.parse_args()

    sys.stderr.write(f"loading queries {args.queries}\n")
    q_tax = lookup.load_lookup(args.queries)
    sys.stderr.write(f"loading targets {args.targets}\n")
    t_tax = lookup.load_lookup(args.targets)
    query_set = reachable_query_set(q_tax, t_tax)
    sys.stderr.write(f"|Q|={len(query_set)}\n")

    scores = condense_strummer(args.hits, query_set, q_tax)
    sys.stderr.write(f"writing {args.output}\n")
    hits_to_topcat.write_tcat(
        args.output,
        scores,
        algo=args.algo,
        truth="topsf",
        fields_spec="1,2,3",
        scores_are_evalues=True,
        n_possible_tp=len(query_set),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
