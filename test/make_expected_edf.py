#!/usr/bin/env python3
"""
make_expected_edf.py -- brute-force reference .edf generator (no homval imports).

Reads test.lookup and test_hits.tsv from the same directory, classifies each
hit under the README truth rules, and writes expected.superfamily.edf and
expected.superfamilyx.edf for comparison with hits_to_edf.py.

All logic is inlined and transparent; data is tiny so everything is O(n^2)
or smaller.  If this agrees with hits_to_edf.py, both implement the same spec.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOOKUP = os.path.join(HERE, "test.lookup")
HITS = os.path.join(HERE, "test_hits.tsv")

# Higher score = stronger hit (--score mode)
SCORES_ARE_EVALUES = False


# ---------------------------------------------------------------------------
# Lookup parsing (family_id = class.fold.superfamily.family)
# ---------------------------------------------------------------------------

def parse_domain(label):
    """Strip suffix after '/' from hit-file label."""
    return label.split("/")[0]


def load_lookup(path):
    dom2fam = {}
    dom2sf = {}
    dom2fold = {}
    with open(path, encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_no}: need 2 tab fields")
            dom, fam = parts[0], parts[1]
            flds = fam.split(".")
            if len(flds) < 4:
                raise ValueError(f"{path}:{line_no}: bad family_id {fam!r}")
            fold = flds[0] + "." + flds[1]
            sf = fold + "." + flds[2]
            dom2fam[dom] = fam
            dom2sf[dom] = sf
            dom2fold[dom] = fold
    return dom2fam, dom2sf, dom2fold


# ---------------------------------------------------------------------------
# Truth rules (verbatim from README)
# ---------------------------------------------------------------------------

def pair_ignore(q, t, truth, dom2fam, dom2sf, dom2fold):
    if q == t:
        return True
    if q not in dom2fam or t not in dom2fam:
        return True
    if truth == "superfamilyx":
        if dom2fold[q] == dom2fold[t] and dom2sf[q] != dom2sf[t]:
            return True
    return False


def pair_is_tp(q, t, truth, dom2fam, dom2sf, dom2fold):
    if truth in ("superfamily", "superfamilyx"):
        return dom2sf[q] == dom2sf[t]
    if truth == "fold":
        return dom2fold[q] == dom2fold[t]
    if truth == "family":
        return dom2fam[q] == dom2fam[t]
    raise ValueError(truth)


# ---------------------------------------------------------------------------
# Brute-force denominators: every ordered pair (q, t), q != t
# ---------------------------------------------------------------------------

def count_pair_denominators(dom2fam, dom2sf, dom2fold, truth):
    doms = sorted(dom2fam.keys())
    n_tp = 0
    n_fp = 0
    for q in doms:
        for t in doms:
            if pair_ignore(q, t, truth, dom2fam, dom2sf, dom2fold):
                continue
            if pair_is_tp(q, t, truth, dom2fam, dom2sf, dom2fold):
                n_tp += 1
            else:
                n_fp += 1
    return len(doms), n_tp, n_fp


# ---------------------------------------------------------------------------
# Read hits; bucket into histogram[score] = [n_tp, n_fp]
# ---------------------------------------------------------------------------

def read_hits_histogram(path, truth, dom2fam, dom2sf, dom2fold):
    hist = {}  # score -> [n_tp, n_fp]
    per_query = {}  # query -> list of (score, is_tp)
    n_ignored = 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            flds = line.split()
            if len(flds) < 3:
                continue
            q = parse_domain(flds[0])
            t = parse_domain(flds[1])
            score = float(flds[2])
            if pair_ignore(q, t, truth, dom2fam, dom2sf, dom2fold):
                n_ignored += 1
                continue
            if score not in hist:
                hist[score] = [0, 0]
            is_tp = pair_is_tp(q, t, truth, dom2fam, dom2sf, dom2fold)
            if is_tp:
                hist[score][0] += 1
            else:
                hist[score][1] += 1
            per_query.setdefault(q, []).append((score, is_tp))
    return hist, per_query, n_ignored


def sort_scores_best_first(scores):
    return sorted(scores, reverse=not SCORES_ARE_EVALUES)


# ---------------------------------------------------------------------------
# Sum3 (FastBench::BenchPairCVE): coverage when error <= threshold, recorded
# at last point before each score change
# ---------------------------------------------------------------------------

def compute_sum3(curve_rows, ndom, n_possible_tp):
    """curve_rows: list of (score, cum_tp, cum_fp) best-first."""
    if n_possible_tp <= 0 or not curve_rows:
        return None, None, None, None

    sepq01 = sepq1 = sepq10 = None
    last_score = None
    for score, cum_tp, cum_fp in curve_rows:
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

    final_cov = curve_rows[-1][1] / n_possible_tp
    if sepq01 is None:
        sepq01 = final_cov
    if sepq1 is None:
        sepq1 = final_cov
    if sepq10 is None:
        sepq10 = final_cov
    sum3 = sepq01 * 2 + sepq1 * 1.5 + sepq10
    return sepq01, sepq1, sepq10, sum3


def compute_sffp(per_query, n_possible_tp):
    """Sensitivity above first false positive (per-query).

    For each query, find its best FP score, then count TPs with strictly
    better score.  Sum across queries, divide by N_possible_tp.
    """
    if n_possible_tp <= 0 or not per_query:
        return None
    total_tp_above = 0
    for _q, hits in per_query.items():
        tp_scores = []
        best_fp_score = None
        for score, is_tp in hits:
            if is_tp:
                tp_scores.append(score)
            else:
                if best_fp_score is None or score > best_fp_score:
                    best_fp_score = score
        if best_fp_score is None:
            total_tp_above += len(tp_scores)
        else:
            for s in tp_scores:
                if s > best_fp_score:
                    total_tp_above += 1
    return total_tp_above / n_possible_tp


def build_edf(truth, dom2fam, dom2sf, dom2fold, hist, per_query):
    ndom, n_possible_tp, n_possible_fp = count_pair_denominators(
        dom2fam, dom2sf, dom2fold, truth
    )
    ordered = sort_scores_best_first(list(hist.keys()))
    cum_tp = 0
    cum_fp = 0
    body_lines = []
    curve_rows = []
    for score in ordered:
        n_tp, n_fp = hist[score]
        cum_tp += n_tp
        cum_fp += n_fp
        curve_rows.append((score, cum_tp, cum_fp))
        body_lines.append(f"{score:g}\t{n_tp}\t{n_fp}\t{cum_tp}\t{cum_fp}")

    sepq01, sepq1, sepq10, sum3 = compute_sum3(curve_rows, ndom, n_possible_tp)
    sffp = compute_sffp(per_query, n_possible_tp)
    score_dir = "lower_better" if SCORES_ARE_EVALUES else "higher_better"

    lines = [
        f"# algo=brute_force_reference",
        f"# truth={truth}",
        f"# fields=1,2,3",
        f"# score_direction={score_dir}",
        f"# ndom={ndom}",
        f"# N_possible_tp={n_possible_tp}",
        f"# N_possible_fp={n_possible_fp}",
        f"# SEPQ0.1={sepq01:.3f} SEPQ1={sepq1:.3f} SEPQ10={sepq10:.3f} Sum3={sum3:.3f} SFFP={sffp:.3f}",
        "score\tn_tp\tn_fp\tcum_tp\tcum_fp",
    ]
    lines.extend(body_lines)
    return "\n".join(lines) + "\n"


def main():
    dom2fam, dom2sf, dom2fold = load_lookup(LOOKUP)
    for truth in ("superfamily", "superfamilyx"):
        hist, per_query, n_ignored = read_hits_histogram(
            HITS, truth, dom2fam, dom2sf, dom2fold
        )
        out_path = os.path.join(HERE, f"expected.{truth}.edf")
        content = build_edf(truth, dom2fam, dom2sf, dom2fold, hist, per_query)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(content)
        sys.stderr.write(
            f"wrote {out_path} ({len(hist)} scores, {n_ignored} ignored)\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
