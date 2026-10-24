#!/usr/bin/env python3
"""
scop40t_topcat.py -- topsf summary for structure methods on the scop40t database.

For each query in scop40q whose superfamily has another family in scop40t,
the superfamily call is the superfamily of the single best hit to a domain
in scop40t. Other hits are ignored. Same-family targets are ignored.

One streaming pass per hits file. Only the current best score per query is kept.

Example:
  python scop40t_topcat.py \\
      --queries /mnt/c/src/strummer_scop/lookup/test.lookup \\
      --targets /mnt/c/src/strummer_scop/lookup/train.lookup \\
      --big-hits ../big_hits \\
      --output-dir ../topcat
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import hits_to_topcat
import lookup
import strummer_topcat

# algo, hits filename, scores_are_evalues (lower better)
METHODS = (
    ("reseek3_sensitive", "reseek3_sensitive_superfamily.scop40.hits", True),
    ("foldseek", "foldseek.scop40.hits", True),
    ("dali", "dali.scop40.hits", False),
    ("tm", "tm.scop40.hits", False),
)

PROGRESS_EVERY = 5_000_000


def scan_hits(
    hits_path: str,
    query_set: frozenset[str],
    q_tax: lookup.Taxonomy,
    t_tax: lookup.Taxonomy,
    scores_are_evalues: bool,
) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    """Single pass: best scop40t target per scop40q query."""
    best: Dict[str, List] = {}
    q_sf = q_tax.dom2sf
    q_fam = q_tax.dom2fam
    t_sf = t_tax.dom2sf
    t_fam = t_tax.dom2fam
    parse_domain = lookup.parse_domain
    note_top = strummer_topcat.note_top
    n_lines = 0
    n_kept = 0

    with open(hits_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            flds = line.split()
            if len(flds) < 3:
                continue
            n_lines += 1
            if n_lines % PROGRESS_EVERY == 0:
                sys.stderr.write(f"{hits_path}: {n_lines} lines\n")
            q = parse_domain(flds[0])
            if q not in query_set:
                continue
            t = parse_domain(flds[1])
            tsf = t_sf.get(t)
            if tsf is None or t == q or t_fam[t] == q_fam[q]:
                continue
            try:
                score = float(flds[2])
            except ValueError:
                continue
            note_top(best, q, score, tsf == q_sf[q], scores_are_evalues)
            n_kept += 1

    sys.stderr.write(
        f"{hits_path}: lines {n_lines}, kept {n_kept}, "
        f"queries with a hit {len(best)}\n"
    )
    return strummer_topcat.scores_from_best(query_set, best)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Top-hit topsf summaries for methods searched against scop40t."
    )
    ap.add_argument("--queries", required=True, help="scop40q lookup (test.lookup)")
    ap.add_argument("--targets", required=True, help="scop40t lookup (train.lookup)")
    ap.add_argument("--big-hits", required=True, help="Directory of method hit files")
    ap.add_argument("--output-dir", required=True, help="Directory for METHOD.scop40t.topcat")
    args = ap.parse_args()

    sys.stderr.write(f"loading queries {args.queries}\n")
    q_tax = lookup.load_lookup(args.queries)
    sys.stderr.write(f"loading targets {args.targets}\n")
    t_tax = lookup.load_lookup(args.targets)
    query_set = strummer_topcat.reachable_query_set(q_tax, t_tax)
    n_q = len(query_set)
    sys.stderr.write(f"|Q|={n_q}\n")

    os.makedirs(args.output_dir, exist_ok=True)
    for algo, filename, scores_are_evalues in METHODS:
        hits_path = os.path.join(args.big_hits, filename)
        out_path = os.path.join(args.output_dir, f"{algo}.scop40t.topcat")
        sys.stderr.write(f"reading {hits_path}\n")
        scores = scan_hits(
            hits_path, query_set, q_tax, t_tax, scores_are_evalues
        )
        sys.stderr.write(f"writing {out_path}\n")
        hits_to_topcat.write_tcat(
            out_path,
            scores,
            algo=algo,
            truth="topsf",
            fields_spec="1,2,3",
            scores_are_evalues=scores_are_evalues,
            n_possible_tp=n_q,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
