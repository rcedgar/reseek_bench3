#!/usr/bin/env python3
"""
build_derived_info.py -- Stage 0: precompute denominators and query sets.

Reads a lookup TSV and writes a derived_info JSON file used by Stage 1 scripts.

Example:
  python build_derived_info.py --lookup ../info/scop40x.lookup \
      --output ../derived_info/scop40x.json
"""

from __future__ import annotations

import argparse
import json
import sys

import lookup


TRUTH_STANDARDS = lookup.TRUTH_STANDARDS


def _ordered_pairs(n: int) -> int:
    """Ordered pairs among n items, excluding self: n * (n - 1)."""
    return n * (n - 1)


def count_pair_denominators(tax: lookup.Taxonomy, truth: str) -> tuple[int, int]:
    """
    Count N_possible_tp and N_possible_fp for ordered pairs (q, t), q != t.

    Matches the pair-list sweep in hits_to_edf.py.
    """
    n = tax.ndom

    if truth in ("fold", "superfamily", "family"):
        if truth == "fold":
            groups = tax.fold2doms
        elif truth == "superfamily":
            groups = tax.sf2doms
        else:
            groups = tax.fam2doms
        n_tp = sum(_ordered_pairs(len(doms)) for doms in groups.values())
        n_fp = _ordered_pairs(n) - n_tp
        return n_tp, n_fp

    if truth == "superfamilyx":
        n_tp = sum(_ordered_pairs(len(doms)) for doms in tax.sf2doms.values())
        n_ignore = 0
        for fold, sfs in tax.fold2sfs.items():
            n_fold = len(tax.fold2doms[fold])
            n_tp_in_fold = sum(
                _ordered_pairs(len(tax.sf2doms[sf])) for sf in sfs
            )
            n_ignore += _ordered_pairs(n_fold) - n_tp_in_fold
        n_fp = _ordered_pairs(n) - n_tp - n_ignore
        return n_tp, n_fp

    raise ValueError(f"unknown truth standard: {truth!r}")


def build_summary(tax: lookup.Taxonomy) -> dict:
    """Assemble derived_info JSON content."""
    qs_topfold = lookup.query_set_topfold(tax)
    qs_topsf = lookup.query_set_topsf(tax)

    out = {
        "ndom": tax.ndom,
        "query_set_topfold": qs_topfold,
        "query_set_topsf": qs_topsf,
        "N_possible_tp": {
            "topfold": len(qs_topfold),
            "topsf": len(qs_topsf),
        },
        "N_possible_fp": {},
    }

    for truth in TRUTH_STANDARDS:
        n_tp, n_fp = count_pair_denominators(tax, truth)
        out["N_possible_tp"][truth] = n_tp
        out["N_possible_fp"][truth] = n_fp

    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Precompute homval denominators and query sets from a lookup file."
    )
    ap.add_argument("--lookup", required=True, help="Path to lookup TSV")
    ap.add_argument("--output", required=True, help="Output derived_info JSON path")
    args = ap.parse_args()

    sys.stderr.write(f"reading {args.lookup}...\n")
    tax = lookup.load_lookup(args.lookup)
    summary = build_summary(tax)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    sys.stderr.write(f"wrote {args.output}\n")
    sys.stderr.write(f"  ndom={summary['ndom']}\n")
    sys.stderr.write(f"  query_set_topfold={len(summary['query_set_topfold'])}\n")
    sys.stderr.write(f"  query_set_topsf={len(summary['query_set_topsf'])}\n")
    for truth in TRUTH_STANDARDS:
        sys.stderr.write(
            f"  {truth}: N_possible_tp={summary['N_possible_tp'][truth]}"
            f" N_possible_fp={summary['N_possible_fp'][truth]}\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
