"""
lookup.py -- SCOP/CATH domain taxonomy from a lookup TSV file.

Lookup file format (no header, tab-separated):
  column 1: domain_id     e.g. d1asha_
  column 2: family_id     e.g. a.1.1.2  (class.fold.superfamily.family)

The family_id string is split to derive superfamily (first three dot fields)
and fold (first two dot fields).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Set


def parse_domain(label: str) -> str:
    """Strip annotation suffix after '/' from a hit-file domain label."""
    return label.split("/")[0]


def _split_fam_id(fam: str) -> tuple[str, str, str]:
    """Return (fold, superfamily, family) from a dotted family_id string."""
    flds = fam.split(".")
    if len(flds) < 4:
        raise ValueError(f"family_id must have 4 dot fields, got: {fam!r}")
    fold = flds[0] + "." + flds[1]
    sf = fold + "." + flds[2]
    return fold, sf, fam


@dataclass
class Taxonomy:
    """Domain annotations and inverse maps built from one lookup file."""

    dom2fam: Dict[str, str] = field(default_factory=dict)
    dom2sf: Dict[str, str] = field(default_factory=dict)
    dom2fold: Dict[str, str] = field(default_factory=dict)

    # Inverse maps: category -> set of domain ids
    fold2doms: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    sf2doms: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    fam2doms: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))

    # For query-set construction
    fold2sfs: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    sf2fams: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))

    @property
    def domains(self) -> FrozenSet[str]:
        return frozenset(self.dom2fam.keys())

    @property
    def ndom(self) -> int:
        return len(self.dom2fam)


def load_lookup(path: str) -> Taxonomy:
    """
    Read a lookup TSV and populate a Taxonomy.

    Duplicate domain ids are rejected.  Blank lines are skipped.
    """
    tax = Taxonomy()
    seen: Set[str] = set()

    with open(path, encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_no}: expected 2 tab fields")
            dom, fam = parts[0], parts[1]
            if dom in seen:
                raise ValueError(f"{path}:{line_no}: duplicate domain {dom!r}")
            seen.add(dom)

            fold, sf, fam = _split_fam_id(fam)
            tax.dom2fam[dom] = fam
            tax.dom2sf[dom] = sf
            tax.dom2fold[dom] = fold
            tax.fold2doms[fold].add(dom)
            tax.sf2doms[sf].add(dom)
            tax.fam2doms[fam].add(dom)
            tax.fold2sfs[fold].add(sf)
            tax.sf2fams[sf].add(fam)

    return tax


def query_set_topfold(tax: Taxonomy) -> List[str]:
    """Domains in folds containing >1 distinct superfamily."""
    out = []
    for dom in sorted(tax.dom2fam):
        fold = tax.dom2fold[dom]
        if len(tax.fold2sfs[fold]) > 1:
            out.append(dom)
    return out


def query_set_topsf(tax: Taxonomy) -> List[str]:
    """Domains in superfamilies containing >1 distinct family."""
    out = []
    for dom in sorted(tax.dom2fam):
        sf = tax.dom2sf[dom]
        if len(tax.sf2fams[sf]) > 1:
            out.append(dom)
    return out


# ---------------------------------------------------------------------------
# Pair-list truth classification (shared by build_derived_info and hits_to_edf)
# ---------------------------------------------------------------------------

TRUTH_STANDARDS = ("fold", "superfamily", "family", "superfamilyx")


def pair_ignore(q: str, t: str, truth: str, tax: Taxonomy) -> bool:
    """Return True if ordered pair (q, t) is excluded from pair-list evaluation."""
    if q == t:
        return True
    if q not in tax.dom2fam or t not in tax.dom2fam:
        return True
    if truth == "superfamilyx":
        qfold, tfold = tax.dom2fold[q], tax.dom2fold[t]
        qsf, tsf = tax.dom2sf[q], tax.dom2sf[t]
        if qfold == tfold and qsf != tsf:
            return True
    return False


def pair_is_tp(q: str, t: str, truth: str, tax: Taxonomy) -> bool:
    """Return True if (q, t) is a true positive under the truth standard."""
    if truth == "fold":
        return tax.dom2fold[q] == tax.dom2fold[t]
    if truth == "superfamily":
        return tax.dom2sf[q] == tax.dom2sf[t]
    if truth == "family":
        return tax.dom2fam[q] == tax.dom2fam[t]
    if truth == "superfamilyx":
        return tax.dom2sf[q] == tax.dom2sf[t]
    raise ValueError(f"unknown truth standard: {truth!r}")


# ---------------------------------------------------------------------------
# Leave-category-out (topcat) hit classification
# ---------------------------------------------------------------------------

TOPCAT_TRUTHS = ("topfold", "topsf")


def topcat_track(
    q: str, t: str, truth: str, tax: Taxonomy
) -> Optional[str]:
    """
    Classify directed hit (q, t) for topcat score tracking.

    Returns None if the hit is ignored (self, unknown domain, or held out).
    Otherwise returns "tp" (updates S_tp) or "xf" (updates S_xf).
    """
    if q == t:
        return None
    if q not in tax.dom2fam or t not in tax.dom2fam:
        return None

    if truth == "topsf":
        if tax.dom2fam[q] == tax.dom2fam[t]:
            return None
        if tax.dom2sf[q] == tax.dom2sf[t]:
            return "tp"
        return "xf"

    if truth == "topfold":
        if tax.dom2sf[q] == tax.dom2sf[t]:
            return None
        if tax.dom2fold[q] == tax.dom2fold[t]:
            return "tp"
        return "xf"

    raise ValueError(f"unknown topcat truth: {truth!r}")
