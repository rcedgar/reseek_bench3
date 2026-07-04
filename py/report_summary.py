#!/usr/bin/env python3
"""
report_summary.py -- tabular report of Sum3 / Top3 statistics from summary files.

Reads .edf (SEPQ/Sum3) and .tcat (TEPQ/Top3) comment headers and writes
NAME.tsv plus a column-aligned NAME.txt.

Each table is headed by metric, truth standard, and reference, e.g.
Sum3>superfamily>scop40 or Top3>topsf>scop40x.

Example:
  python report_summary.py bench_superfamily \\
      ../edf/reseek28.scop40.superfamily ../edf/foldseek.scop40.superfamily \\
      ../edf/dali.scop40.superfamily ../edf/tm.scop40.superfamily
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import common
import lookup

ALL_TRUTHS = frozenset(lookup.TRUTH_STANDARDS + lookup.TOPCAT_TRUTHS)

SUM3_COLUMNS = ("SEPQ0.1", "SEPQ1", "SEPQ10", "Sum3", "algo")
SUM3_KEYS = ("SEPQ0.1", "SEPQ1", "SEPQ10", "Sum3")

SFFP_COLUMNS = ("SFFP", "algo")
SFFP_KEYS = ("SFFP",)

TOP3_COLUMNS = ("TEPQ0.001", "TEPQ0.01", "TEPQ0.1", "Top3", "algo")
TOP3_KEYS = ("TEPQ0.001", "TEPQ0.01", "TEPQ0.1", "Top3")

GroupKey = Tuple[str, str]


@dataclass(frozen=True)
class Row:
    truth: str
    reference: str
    algo: str
    values: Tuple[str, ...]
    sort_key: float


def parse_header_line(content: str) -> List[Tuple[str, str]]:
    """Parse '# key=value' or '# K1=v1 K2=v2 ...' into (key, value) pairs."""
    pairs: List[Tuple[str, str]] = []
    if content.startswith("SEPQ0.1=") or content.startswith("TEPQ0.001="):
        for part in content.split():
            k, v = part.split("=", 1)
            pairs.append((k, v))
    elif "=" in content:
        k, v = content.split("=", 1)
        pairs.append((k, v))
    return pairs


def read_summary_header(path: str) -> dict[str, str]:
    hdr: dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                for k, v in parse_header_line(line[1:].strip()):
                    hdr[k] = v
            elif line.strip():
                break
    return hdr


def basename_parts(path: str) -> List[str]:
    base = Path(path).name
    for ext in (".edf", ".tcat"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base.split(".")


def truth_from_path(path: str, hdr: dict[str, str]) -> str:
    if hdr.get("truth"):
        return hdr["truth"]
    parts = basename_parts(path)
    if parts and parts[-1] in ALL_TRUTHS:
        return parts[-1]
    return ""


def reference_from_path(path: str, hdr: dict[str, str]) -> str:
    if hdr.get("reference"):
        return hdr["reference"]
    parts = basename_parts(path)
    if len(parts) >= 3 and parts[-1] in ALL_TRUTHS:
        return parts[-2]
    return ""


def algo_from_header(hdr: dict[str, str], path: str) -> str:
    fallback = Path(path).name.split(".", 1)[0]
    return common.algo_label_from_header(hdr, fallback=fallback)


def section_subhead(section: str, truth: str, reference: str) -> str:
    return f"{section}>{truth}>{reference}"


def load_row(path: str) -> Optional[List[Tuple[str, Row]]]:
    """Return list of ('sum3'|'sffp'|'top3', Row) or None if no summary stats."""
    hdr = read_summary_header(path)
    truth = truth_from_path(path, hdr)
    reference = reference_from_path(path, hdr)
    algo = algo_from_header(hdr, path)

    results: List[Tuple[str, Row]] = []

    if all(k in hdr for k in SUM3_KEYS):
        values = tuple(hdr[k] for k in SUM3_KEYS)
        results.append(("sum3", Row(
            truth=truth,
            reference=reference,
            algo=algo,
            values=values,
            sort_key=float(hdr["Sum3"]),
        )))

    if all(k in hdr for k in SFFP_KEYS):
        values = tuple(hdr[k] for k in SFFP_KEYS)
        results.append(("sffp", Row(
            truth=truth,
            reference=reference,
            algo=algo,
            values=values,
            sort_key=float(hdr["SFFP"]),
        )))

    if all(k in hdr for k in TOP3_KEYS):
        values = tuple(hdr[k] for k in TOP3_KEYS)
        results.append(("top3", Row(
            truth=truth,
            reference=reference,
            algo=algo,
            values=values,
            sort_key=float(hdr["Top3"]),
        )))

    return results if results else None


def group_rows(
    paths: Sequence[str],
) -> Tuple[Dict[GroupKey, List[Row]], Dict[GroupKey, List[Row]], Dict[GroupKey, List[Row]]]:
    sum3: Dict[GroupKey, List[Row]] = {}
    sffp: Dict[GroupKey, List[Row]] = {}
    top3: Dict[GroupKey, List[Row]] = {}
    buckets = {"sum3": sum3, "sffp": sffp, "top3": top3}
    for path in paths:
        loaded = load_row(path)
        if loaded is None:
            print(f"warning: no Sum3/Top3 stats in {path}", file=sys.stderr)
            continue
        for kind, row in loaded:
            bucket = buckets[kind]
            key = (row.truth, row.reference)
            bucket.setdefault(key, []).append(row)
    for bucket in buckets.values():
        for rows in bucket.values():
            rows.sort(key=lambda r: (r.sort_key, r.algo))
    return sum3, sffp, top3


def sorted_group_keys(groups: Dict[GroupKey, List[Row]]) -> List[GroupKey]:
    return sorted(groups)


def write_section_tsv(
    f,
    section: str,
    columns: Sequence[str],
    groups: Dict[GroupKey, List[Row]],
) -> None:
    for truth, reference in sorted_group_keys(groups):
        f.write(f"# {section_subhead(section, truth, reference)}\n")
        f.write("\t".join(columns) + "\n")
        for row in groups[(truth, reference)]:
            f.write("\t".join((*row.values, row.algo)) + "\n")
        f.write("\n")


def write_tsv(
    out_path: Path,
    sum3: Dict[GroupKey, List[Row]],
    sffp: Dict[GroupKey, List[Row]],
    top3: Dict[GroupKey, List[Row]],
) -> None:
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        if sum3:
            write_section_tsv(f, "Sum3", SUM3_COLUMNS, sum3)
        if sffp:
            write_section_tsv(f, "SFFP", SFFP_COLUMNS, sffp)
        if top3:
            write_section_tsv(f, "Top3", TOP3_COLUMNS, top3)


def column_widths(
    columns: Sequence[str], rows: Sequence[Row]
) -> List[int]:
    widths = [len(c) for c in columns]
    for row in rows:
        for i, value in enumerate(row.values):
            widths[i] = max(widths[i], len(value))
        widths[-1] = max(widths[-1], len(row.algo))
    return widths


def format_table(
    columns: Sequence[str],
    rows: Sequence[Row],
    widths: Sequence[int],
) -> List[str]:
    lines: List[str] = []
    header = "  ".join(c.ljust(widths[i]) for i, c in enumerate(columns))
    lines.append(header)
    for row in rows:
        cells = [
            *(v.rjust(widths[i]) for i, v in enumerate(row.values)),
            row.algo.ljust(widths[-1]),
        ]
        lines.append("  ".join(cells))
    return lines


def write_section_txt(
    lines: List[str],
    section: str,
    columns: Sequence[str],
    groups: Dict[GroupKey, List[Row]],
) -> None:
    for truth, reference in sorted_group_keys(groups):
        rows = groups[(truth, reference)]
        widths = column_widths(columns, rows)
        lines.append(section_subhead(section, truth, reference))
        lines.extend(format_table(columns, rows, widths))
        lines.append("")


def write_txt(
    out_path: Path,
    sum3: Dict[GroupKey, List[Row]],
    sffp: Dict[GroupKey, List[Row]],
    top3: Dict[GroupKey, List[Row]],
) -> None:
    lines: List[str] = []
    if sum3:
        write_section_txt(lines, "Sum3", SUM3_COLUMNS, sum3)
    if sffp:
        write_section_txt(lines, "SFFP", SFFP_COLUMNS, sffp)
    if top3:
        write_section_txt(lines, "Top3", TOP3_COLUMNS, top3)
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines).rstrip() + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Report Sum3 / Top3 statistics from homval summary files."
    )
    ap.add_argument(
        "name",
        help="Output basename (writes NAME.tsv and NAME.txt)",
    )
    ap.add_argument(
        "inputs",
        nargs="+",
        metavar="SUMMARY",
        help="Summary files (.edf or .tcat) as the last arguments",
    )
    args = ap.parse_args()

    if not args.inputs:
        ap.error("at least one summary file is required")

    sum3, sffp, top3 = group_rows(args.inputs)
    if not sum3 and not sffp and not top3:
        print("error: no summary statistics found in any input file", file=sys.stderr)
        return 1

    base = Path(args.name)
    tsv_path = base.with_suffix(".tsv")
    txt_path = base.with_suffix(".txt")
    write_tsv(tsv_path, sum3, sffp, top3)
    write_txt(txt_path, sum3, sffp, top3)
    print(f"wrote {tsv_path} and {txt_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
