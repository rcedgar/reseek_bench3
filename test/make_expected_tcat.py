#!/usr/bin/env python3
"""
make_expected_tcat.py -- brute-force reference .tcat generator (no homval imports).

Reads test.lookup, test_hits.tsv, and derived_info.json from the same directory,
computes S_tp/S_xf and Top3 per truth_standards.md, and writes expected.topsf.tcat
and expected.topfold.tcat for comparison with hits_to_topcat.py.
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOOKUP = os.path.join(HERE, "test.lookup")
HITS = os.path.join(HERE, "test_hits.tsv")
DERIVED = os.path.join(HERE, "derived_info.json")

SCORES_ARE_EVALUES = False


def parse_domain(label):
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


def topcat_track(q, t, truth, dom2fam, dom2sf, dom2fold):
    if q == t:
        return None
    if q not in dom2fam or t not in dom2fam:
        return None
    if truth == "topsf":
        if dom2fam[q] == dom2fam[t]:
            return None
        if dom2sf[q] == dom2sf[t]:
            return "tp"
        return "xf"
    if truth == "topfold":
        if dom2sf[q] == dom2sf[t]:
            return None
        if dom2fold[q] == dom2fold[t]:
            return "tp"
        return "xf"
    raise ValueError(truth)


def better(score, best):
    if best is None:
        return True
    if SCORES_ARE_EVALUES:
        return score < best
    return score > best


def passes_threshold(score, threshold):
    if score is None:
        return False
    if SCORES_ARE_EVALUES:
        return score <= threshold
    return score >= threshold


def at_least_as_good(score, other):
    if SCORES_ARE_EVALUES:
        return score <= other
    return score >= other


def is_tp_at_threshold(s_tp, s_xf, threshold):
    if not passes_threshold(s_tp, threshold):
        return False
    if not passes_threshold(s_xf, threshold):
        return True
    return at_least_as_good(s_tp, s_xf)


def classify_at_threshold(scores, threshold, n_q):
    n_tp = 0
    n_fp_neg = 0
    for s_tp, s_xf in scores.values():
        if is_tp_at_threshold(s_tp, s_xf, threshold):
            n_tp += 1
        if passes_threshold(s_xf, threshold):
            n_fp_neg += 1
    return n_tp / n_q, n_fp_neg / n_q


def collect_thresholds(scores):
    seen = set()
    for s_tp, s_xf in scores.values():
        if s_tp is not None:
            seen.add(s_tp)
        if s_xf is not None:
            seen.add(s_xf)
    return list(seen)


def sort_thresholds_strict_first(thresholds):
    return sorted(thresholds, reverse=not SCORES_ARE_EVALUES)


def loose_threshold():
    return float("inf") if SCORES_ARE_EVALUES else float("-inf")


def compute_top3(scores, n_q):
    taus = sort_thresholds_strict_first(collect_thresholds(scores))
    taus.append(loose_threshold())

    tepq0001 = tepq001 = tepq01 = None
    final_cov = 0.0
    for tau in taus:
        cov, err = classify_at_threshold(scores, tau, n_q)
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


def read_hits_scores(path, truth, query_set, dom2fam, dom2sf, dom2fold):
    scores = {dom: (None, None) for dom in sorted(query_set)}
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
            if q not in query_set:
                continue
            track = topcat_track(q, t, truth, dom2fam, dom2sf, dom2fold)
            if track is None:
                continue
            s_tp, s_xf = scores[q]
            if track == "tp":
                if better(score, s_tp):
                    s_tp = score
            elif better(score, s_xf):
                s_xf = score
            scores[q] = (s_tp, s_xf)
    return scores


def format_score(score):
    return "." if score is None else f"{score:g}"


def build_tcat(truth, scores, n_possible_tp):
    tepq0001, tepq001, tepq01, top3 = compute_top3(scores, n_possible_tp)
    score_dir = "lower_better" if SCORES_ARE_EVALUES else "higher_better"
    lines = [
        "# algo=brute_force_reference",
        f"# truth={truth}",
        "# fields=1,2,3",
        f"# score_direction={score_dir}",
        f"# N_possible_tp={n_possible_tp}",
        (
            f"# TEPQ0.001={tepq0001:.3f} TEPQ0.01={tepq001:.3f} "
            f"TEPQ0.1={tepq01:.3f} Top3={top3:.3f}"
        ),
        "domain\tS_tp\tS_xf",
    ]
    for dom in sorted(scores):
        s_tp, s_xf = scores[dom]
        lines.append(f"{dom}\t{format_score(s_tp)}\t{format_score(s_xf)}")
    return "\n".join(lines) + "\n"


def main():
    dom2fam, dom2sf, dom2fold = load_lookup(LOOKUP)
    with open(DERIVED, encoding="utf-8") as f:
        info = json.load(f)

    for truth in ("topsf", "topfold"):
        query_set = frozenset(info[f"query_set_{truth}"])
        n_possible_tp = int(info["N_possible_tp"][truth])
        scores = read_hits_scores(
            HITS, truth, query_set, dom2fam, dom2sf, dom2fold
        )
        out_path = os.path.join(HERE, f"expected.{truth}.tcat")
        content = build_tcat(truth, scores, n_possible_tp)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(content)
        sys.stderr.write(f"wrote {out_path} (|Q|={n_possible_tp})\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
