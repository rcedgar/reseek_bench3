"""
common.py -- small shared helpers for homval scripts.

Kept minimal so the homval directory can ship as self-contained supplementary
material with lookup.py as the other shared module.
"""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Score direction (--evalue or --score, mutually exclusive; no auto-detection)
# ---------------------------------------------------------------------------

def add_score_direction_args(ap: argparse.ArgumentParser) -> argparse._MutuallyExclusiveGroup:
    """Add --evalue and --score flags; exactly one is required."""
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument(
        "--evalue",
        action="store_true",
        help="Lower score is better (E-values, P-values)",
    )
    g.add_argument(
        "--score",
        action="store_true",
        help="Higher score is better (bit scores, similarity)",
    )
    return g


def scores_are_evalues_from_args(args: argparse.Namespace) -> bool:
    return bool(args.evalue)


def better(score: float, best: Optional[float], scores_are_evalues: bool) -> bool:
    """
    Return True if score is strictly better than best.

    best is None means no prior score; any score beats None.
    """
    if best is None:
        return True
    if scores_are_evalues:
        return score < best
    return score > best


def passes_threshold(
    score: Optional[float], threshold: float, scores_are_evalues: bool
) -> bool:
    """Return True if score reaches the acceptance threshold."""
    if score is None:
        return False
    if scores_are_evalues:
        return score <= threshold
    return score >= threshold


def at_least_as_good(
    score: float, other: float, scores_are_evalues: bool
) -> bool:
    """Return True if score is at least as good as other (ties count as good)."""
    if scores_are_evalues:
        return score <= other
    return score >= other


# ---------------------------------------------------------------------------
# Hit TSV field indices (1-based in user docs, 0-based internally)
# ---------------------------------------------------------------------------

def parse_fields(spec: str) -> Tuple[int, int, int]:
    """
    Parse 'query,target,score' field spec (1-based) -> 0-based indices.

    Example: '1,2,3' -> (0, 1, 2)
    """
    parts = [p.strip() for p in spec.split(",")]
    if len(parts) != 3:
        raise ValueError(f"fields must have 3 comma-separated 1-based indices, got {spec!r}")
    q, t, s = (int(p) - 1 for p in parts)
    if min(q, t, s) < 0:
        raise ValueError(f"field indices must be >= 1, got {spec!r}")
    return q, t, s


def fields_spec_1based(q_idx: int, t_idx: int, s_idx: int) -> str:
    """Format 0-based indices as 1-based fields spec for .edf header."""
    return f"{q_idx + 1},{t_idx + 1},{s_idx + 1}"


# ---------------------------------------------------------------------------
# Summary file algo / reference headers (.edf, .tcat)
# ---------------------------------------------------------------------------

def parse_algo_reference(label: str) -> Tuple[str, str]:
    """
    Split a hits basename or 'algo.reference' label into (algo, reference).

    Strips a trailing .tsv or .txt extension before splitting on the first dot.
    """
    base = label.replace("\\", "/").split("/")[-1]
    for ext in (".tsv", ".txt"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    if "." in base:
        algo, reference = base.split(".", 1)
        return algo, reference
    return base, ""


def resolve_algo_reference(
    hits_path: str,
    algo: str = "",
    reference: Optional[str] = None,
) -> Tuple[str, str]:
    """Resolve algo and reference for summary file headers."""
    default_algo, default_ref = parse_algo_reference(
        hits_path.replace("\\", "/").split("/")[-1]
    )
    resolved_algo = algo or default_algo
    resolved_ref = default_ref if reference is None else reference
    return resolved_algo, resolved_ref


def write_algo_reference_header(f, algo: str, reference: str) -> None:
    """Write # algo= and # reference= comment lines."""
    f.write(f"# algo={algo}\n")
    if reference:
        f.write(f"# reference={reference}\n")


def algo_label_from_header(hdr: Dict[str, str], fallback: str = "") -> str:
    """Legend label from summary header; handles legacy combined algo.ref."""
    algo = hdr.get("algo", fallback)
    if "reference" not in hdr and "." in algo:
        return algo.split(".", 1)[0]
    return algo


# ---------------------------------------------------------------------------
# derived_info JSON (Stage 0 output, read by Stage 1)
# ---------------------------------------------------------------------------

def load_derived_info(path: str) -> Dict[str, Any]:
    """Load derived_info JSON written by build_derived_info.py."""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def pair_denominators(info: Dict[str, Any], truth: str) -> Tuple[int, int, int]:
    """
    Return (ndom, N_possible_tp, N_possible_fp) for pair-list truth standard.
    """
    if truth not in info.get("N_possible_tp", {}):
        raise KeyError(f"truth {truth!r} not in derived_info N_possible_tp")
    if truth not in info.get("N_possible_fp", {}):
        raise KeyError(f"truth {truth!r} not in derived_info N_possible_fp")
    return (
        int(info["ndom"]),
        int(info["N_possible_tp"][truth]),
        int(info["N_possible_fp"][truth]),
    )


def topcat_denominators(info: Dict[str, Any], truth: str) -> int:
    """Return N_possible_tp (query set size |Q|) for a topcat truth standard."""
    if truth not in info.get("N_possible_tp", {}):
        raise KeyError(f"truth {truth!r} not in derived_info N_possible_tp")
    return int(info["N_possible_tp"][truth])


def topcat_query_set(info: Dict[str, Any], truth: str) -> frozenset[str]:
    """Return domain ids in the topcat query set for this truth standard."""
    key = f"query_set_{truth}"
    if key not in info:
        raise KeyError(f"{key!r} not in derived_info")
    return frozenset(info[key])
