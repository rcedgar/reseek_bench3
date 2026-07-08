#!/usr/bin/env python3
"""
report_summary.py -- tabular report of Sum3 / Top3 statistics from summary files.

Reads .edf (SEPQ/Sum3) and .tcat (TEPQ/Top3) comment headers and writes
NAME.tsv plus a column-aligned NAME.txt.

Each table is headed by metric, truth standard, and reference, e.g.
Sum3>superfamily>scop40 or Top3>topsf>scop40x.

A final Rank section lists every Sum3, SFFP, PR90, and Top3 benchmark as a row,
algos as columns (sorted by decreasing median rank), rank 1=best with
competition tie-breaking (1, 2, 2, 4, ...), and a median footer row per algo.
RankNoSFFP is the same but omits SFFP rows.

RankPValues and RankNoSFFPPValues give two-sided Wilcoxon signed-rank p-values
for each unordered algo pair on paired competition ranks (same benchmark rows
as the corresponding Rank table; only rows where both algos appear). The .txt
tables are symmetric off-diagonal matrices: row vs column with > better,
< worse, ~ if p>0.05; mirrored cells share the p-value with > and < reversed.

Example:
  python report_summary.py bench_superfamily \\
      ../edf/reseek28.scop40.superfamily ../edf/foldseek.scop40.superfamily \\
      ../edf/dali.scop40.superfamily ../edf/tm.scop40.superfamily
"""

from __future__ import annotations

import argparse
import math
import statistics
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

PR90_COLUMNS = ("PR90", "algo")
PR90_KEYS = ("PR90",)

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


def header_n_possible_tp(hdr: dict[str, str]) -> Optional[int]:
    raw = hdr.get("N_possible_tp")
    if raw is None:
        return None
    try:
        return int(float(raw))
    except ValueError:
        return None


def load_row(path: str) -> Optional[List[Tuple[str, Row]]]:
    """Return list of ('sum3'|'sffp'|'top3', Row) or None if no summary stats."""
    hdr = read_summary_header(path)
    truth = truth_from_path(path, hdr)
    reference = reference_from_path(path, hdr)
    algo = algo_from_header(hdr, path)
    n_possible_tp = header_n_possible_tp(hdr)

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

    if all(k in hdr for k in PR90_KEYS):
        values = tuple(hdr[k] for k in PR90_KEYS)
        results.append(("pr90", Row(
            truth=truth,
            reference=reference,
            algo=algo,
            values=values,
            sort_key=float(hdr["PR90"]),
        )))

    if all(k in hdr for k in TOP3_KEYS) and n_possible_tp != 0:
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
) -> Tuple[
    Dict[GroupKey, List[Row]],
    Dict[GroupKey, List[Row]],
    Dict[GroupKey, List[Row]],
    Dict[GroupKey, List[Row]],
]:
    sum3: Dict[GroupKey, List[Row]] = {}
    sffp: Dict[GroupKey, List[Row]] = {}
    pr90: Dict[GroupKey, List[Row]] = {}
    top3: Dict[GroupKey, List[Row]] = {}
    buckets = {"sum3": sum3, "sffp": sffp, "pr90": pr90, "top3": top3}
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
    return sum3, sffp, pr90, top3


def sorted_group_keys(groups: Dict[GroupKey, List[Row]]) -> List[GroupKey]:
    return sorted(groups)


def competition_ranks(rows: Sequence[Row]) -> Dict[str, int]:
    """Rank algos 1=best; ties share rank, next rank skips (1, 2, 2, 4, ...)."""
    ordered = sorted(rows, key=lambda r: (-r.sort_key, r.algo))
    ranks: Dict[str, int] = {}
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j].sort_key == ordered[i].sort_key:
            j += 1
        rank = i + 1
        for row in ordered[i:j]:
            ranks[row.algo] = rank
        i = j
    return ranks


@dataclass(frozen=True)
class RankTable:
    row_labels: Tuple[str, ...]
    algos: Tuple[str, ...]
    ranks: Tuple[Tuple[str, ...], ...]
    medians: Tuple[str, ...]


@dataclass(frozen=True)
class RankPairPValue:
    algo_a: str
    algo_b: str
    n_cells: int
    p_value: Optional[float]


@dataclass(frozen=True)
class RankPValueTable:
    pairs: Tuple[RankPairPValue, ...]


def collect_rank_rows(
    sections: Sequence[Tuple[str, Dict[GroupKey, List[Row]]]],
) -> Tuple[List[str], List[Dict[str, int]]]:
    row_labels: List[str] = []
    rank_rows: List[Dict[str, int]] = []
    for section, groups in sections:
        for truth, reference in sorted_group_keys(groups):
            rows = groups[(truth, reference)]
            row_labels.append(section_subhead(section, truth, reference))
            rank_rows.append(competition_ranks(rows))
    return row_labels, rank_rows


def build_rank_table(
    sections: Sequence[Tuple[str, Dict[GroupKey, List[Row]]]],
) -> Optional[RankTable]:
    row_labels, rank_rows = collect_rank_rows(sections)
    if not row_labels:
        return None

    algo_ranks: Dict[str, List[int]] = {}
    for rank_row in rank_rows:
        for algo, rank in rank_row.items():
            algo_ranks.setdefault(algo, []).append(rank)

    algos = sorted(
        algo_ranks,
        key=lambda a: (-statistics.median(algo_ranks[a]), a),
    )
    medians = tuple(f"{statistics.median(algo_ranks[a]):g}" for a in algos)
    ranks = tuple(
        tuple(
            str(rank_row[algo]) if algo in rank_row else ""
            for algo in algos
        )
        for rank_row in rank_rows
    )
    return RankTable(
        row_labels=tuple(row_labels),
        algos=tuple(algos),
        ranks=ranks,
        medians=medians,
    )


def _average_ranks(values: Sequence[float]) -> List[float]:
    """Average ranks for ties (1-based)."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i + 1
        while j < n and values[order[j]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    return ranks


def _normal_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def wilcoxon_signed_rank_pvalue(
    xs: Sequence[float], ys: Sequence[float]
) -> Optional[float]:
    """Two-sided Wilcoxon signed-rank test; discards zero differences."""
    if len(xs) != len(ys) or len(xs) == 0:
        return None
    diffs = [float(x) - float(y) for x, y in zip(xs, ys)]
    nonzero = [(d) for d in diffs if d != 0.0]
    n = len(nonzero)
    if n == 0:
        return 1.0
    if n == 1:
        return 1.0

    abs_diffs = [abs(d) for d in nonzero]
    ranks = _average_ranks(abs_diffs)
    w_plus = sum(r for r, d in zip(ranks, nonzero) if d > 0.0)

    tie_counts: Dict[float, int] = {}
    for v in abs_diffs:
        tie_counts[v] = tie_counts.get(v, 0) + 1
    tie_term = sum(t ** 3 - t for t in tie_counts.values() if t > 1)

    mu = n * (n + 1) / 4.0
    var = n * (n + 1) * (2 * n + 1) / 24.0 - tie_term / 48.0
    if var <= 0.0:
        return 1.0
    z = (abs(w_plus - mu) - 0.5) / math.sqrt(var)
    return min(1.0, 2.0 * _normal_sf(z))


def build_rank_pvalue_table(
    sections: Sequence[Tuple[str, Dict[GroupKey, List[Row]]]],
) -> Optional[RankPValueTable]:
    _, rank_rows = collect_rank_rows(sections)
    if not rank_rows:
        return None
    algos = sorted({algo for row in rank_rows for algo in row})
    if len(algos) < 2:
        return None

    pairs: List[RankPairPValue] = []
    for i, algo_a in enumerate(algos):
        for algo_b in algos[i + 1:]:
            xs: List[float] = []
            ys: List[float] = []
            for row in rank_rows:
                if algo_a in row and algo_b in row:
                    xs.append(float(row[algo_a]))
                    ys.append(float(row[algo_b]))
            p_value = (
                wilcoxon_signed_rank_pvalue(xs, ys)
                if len(xs) >= 2
                else None
            )
            pairs.append(RankPairPValue(
                algo_a=algo_a,
                algo_b=algo_b,
                n_cells=len(xs),
                p_value=p_value,
            ))
    return RankPValueTable(pairs=tuple(pairs))


def _fmt_pvalue(p_value: Optional[float]) -> str:
    if p_value is None:
        return "."
    return f"{p_value:.4g}"


def write_rank_pvalue_tsv(f, section: str, table: RankPValueTable) -> None:
    f.write(f"# {section}\n")
    f.write("\t".join(("algo_a", "algo_b", "n_cells", "p_value")) + "\n")
    for row in table.pairs:
        f.write("\t".join((
            row.algo_a,
            row.algo_b,
            str(row.n_cells),
            _fmt_pvalue(row.p_value),
        )) + "\n")
    f.write("\n")


def _rank_pvalue_lookup(
    table: RankPValueTable,
) -> Dict[Tuple[str, str], Optional[float]]:
    return {(p.algo_a, p.algo_b): p.p_value for p in table.pairs}


def _directional_pvalue_cell(
    row_algo: str,
    col_algo: str,
    pmap: Dict[Tuple[str, str], Optional[float]],
    rank_rows: Sequence[Dict[str, int]],
) -> str:
    """Prefix p-value: row>col better, row<col worse, ~ if p>0.05."""
    a, b = (
        (row_algo, col_algo)
        if row_algo < col_algo
        else (col_algo, row_algo)
    )
    p_value = pmap.get((a, b))
    if p_value is None:
        return "."

    diffs: List[float] = []
    for row in rank_rows:
        if row_algo in row and col_algo in row:
            diffs.append(float(row[row_algo]) - float(row[col_algo]))
    if not diffs:
        return "."

    if p_value > 0.05:
        prefix = "~"
    else:
        med_diff = statistics.median(diffs)
        if med_diff < 0:
            prefix = ">"
        elif med_diff > 0:
            prefix = "<"
        else:
            prefix = "~"

    return f"{prefix}{_fmt_pvalue(p_value)}"


def format_rank_pvalue_table(
    table: RankPValueTable,
    algos: Sequence[str],
    rank_rows: Sequence[Dict[str, int]],
) -> List[str]:
    """Symmetric off-diagonal: row vs col; mirror cell swaps > and <."""
    if len(algos) < 2:
        return []

    pmap = _rank_pvalue_lookup(table)
    columns = ("", *algos)
    widths = [len(c) for c in columns]
    matrix: List[Tuple[str, ...]] = []

    for i, row_algo in enumerate(algos):
        cells: List[str] = [row_algo]
        for j, col_algo in enumerate(algos):
            if i == j:
                cells.append("")
            else:
                cells.append(_directional_pvalue_cell(
                    row_algo, col_algo, pmap, rank_rows
                ))
        matrix.append(tuple(cells))

    for row in matrix:
        for i, cell in enumerate(row):
            if cell:
                widths[i] = max(widths[i], len(cell))
    for i, algo in enumerate(algos, start=1):
        widths[i] = max(widths[i], len(algo))

    lines = [
        "  ".join(
            columns[i].ljust(widths[i]) if i == 0 else columns[i].rjust(widths[i])
            for i in range(len(columns))
        ),
    ]
    for row in matrix:
        cells = [row[0].ljust(widths[0])]
        cells.extend(
            cell.rjust(widths[i + 1]) if cell else " " * widths[i + 1]
            for i, cell in enumerate(row[1:])
        )
        lines.append("  ".join(cells))
    return lines


def rank_sections(
    sum3: Dict[GroupKey, List[Row]],
    sffp: Dict[GroupKey, List[Row]],
    pr90: Dict[GroupKey, List[Row]],
    top3: Dict[GroupKey, List[Row]],
    *,
    include_sffp: bool,
) -> Sequence[Tuple[str, Dict[GroupKey, List[Row]]]]:
    sections: List[Tuple[str, Dict[GroupKey, List[Row]]]] = []
    if sum3:
        sections.append(("Sum3", sum3))
    if include_sffp and sffp:
        sections.append(("SFFP", sffp))
    if pr90:
        sections.append(("PR90", pr90))
    if top3:
        sections.append(("Top3", top3))
    return sections


def write_rank_tsv(f, section: str, rank_table: RankTable) -> None:
    f.write(f"# {section}\n")
    f.write("\t".join(("metric", *rank_table.algos)) + "\n")
    for label, row in zip(rank_table.row_labels, rank_table.ranks):
        f.write("\t".join((label, *row)) + "\n")
    f.write("\t".join(("median", *rank_table.medians)) + "\n")
    f.write("\n")


def rank_column_widths(rank_table: RankTable) -> List[int]:
    columns = ("metric", *rank_table.algos)
    widths = [len(c) for c in columns]
    for label, row in zip(rank_table.row_labels, rank_table.ranks):
        widths[0] = max(widths[0], len(label))
        for i, value in enumerate(row, start=1):
            if value:
                widths[i] = max(widths[i], len(value))
    for i, median in enumerate(rank_table.medians, start=1):
        widths[i] = max(widths[i], len(median))
    widths[0] = max(widths[0], len("median"))
    return widths


def format_rank_table(rank_table: RankTable) -> List[str]:
    widths = rank_column_widths(rank_table)
    columns = ("metric", *rank_table.algos)
    lines = [
        "  ".join(c.ljust(widths[i]) for i, c in enumerate(columns)),
    ]
    for label, row in zip(rank_table.row_labels, rank_table.ranks):
        cells = [label.ljust(widths[0])]
        cells.extend(
            value.rjust(widths[i + 1]) if value else " " * widths[i + 1]
            for i, value in enumerate(row)
        )
        lines.append("  ".join(cells))
    median_cells = ["median".ljust(widths[0])]
    median_cells.extend(m.rjust(widths[i + 1]) for i, m in enumerate(rank_table.medians))
    lines.append("  ".join(median_cells))
    return lines


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
    pr90: Dict[GroupKey, List[Row]],
    top3: Dict[GroupKey, List[Row]],
) -> None:
    rank_table = build_rank_table(
        rank_sections(sum3, sffp, pr90, top3, include_sffp=True)
    )
    rank_no_sffp = build_rank_table(
        rank_sections(sum3, sffp, pr90, top3, include_sffp=False)
    )
    rank_pvalues = build_rank_pvalue_table(
        rank_sections(sum3, sffp, pr90, top3, include_sffp=True)
    )
    rank_no_sffp_pvalues = build_rank_pvalue_table(
        rank_sections(sum3, sffp, pr90, top3, include_sffp=False)
    )
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        if sum3:
            write_section_tsv(f, "Sum3", SUM3_COLUMNS, sum3)
        if sffp:
            write_section_tsv(f, "SFFP", SFFP_COLUMNS, sffp)
        if pr90:
            write_section_tsv(f, "PR90", PR90_COLUMNS, pr90)
        if top3:
            write_section_tsv(f, "Top3", TOP3_COLUMNS, top3)
        if rank_table is not None:
            write_rank_tsv(f, "Rank", rank_table)
        if rank_no_sffp is not None:
            write_rank_tsv(f, "RankNoSFFP", rank_no_sffp)
        if rank_pvalues is not None:
            write_rank_pvalue_tsv(f, "RankPValues", rank_pvalues)
        if rank_no_sffp_pvalues is not None:
            write_rank_pvalue_tsv(f, "RankNoSFFPPValues", rank_no_sffp_pvalues)


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
    pr90: Dict[GroupKey, List[Row]],
    top3: Dict[GroupKey, List[Row]],
) -> None:
    lines: List[str] = []
    if sum3:
        write_section_txt(lines, "Sum3", SUM3_COLUMNS, sum3)
    if sffp:
        write_section_txt(lines, "SFFP", SFFP_COLUMNS, sffp)
    if pr90:
        write_section_txt(lines, "PR90", PR90_COLUMNS, pr90)
    if top3:
        write_section_txt(lines, "Top3", TOP3_COLUMNS, top3)
    rank_sections_with_sffp = rank_sections(
        sum3, sffp, pr90, top3, include_sffp=True
    )
    rank_sections_no_sffp = rank_sections(
        sum3, sffp, pr90, top3, include_sffp=False
    )
    rank_table = build_rank_table(rank_sections_with_sffp)
    rank_no_sffp = build_rank_table(rank_sections_no_sffp)
    rank_pvalues = build_rank_pvalue_table(rank_sections_with_sffp)
    rank_no_sffp_pvalues = build_rank_pvalue_table(rank_sections_no_sffp)
    _, rank_rows_with_sffp = collect_rank_rows(rank_sections_with_sffp)
    _, rank_rows_no_sffp = collect_rank_rows(rank_sections_no_sffp)
    if rank_table is not None:
        lines.append("Rank")
        lines.extend(format_rank_table(rank_table))
        lines.append("")
    if rank_no_sffp is not None:
        lines.append("RankNoSFFP")
        lines.extend(format_rank_table(rank_no_sffp))
        lines.append("")
    if rank_pvalues is not None and rank_table is not None:
        lines.append("RankPValues")
        lines.extend(format_rank_pvalue_table(
            rank_pvalues, rank_table.algos, rank_rows_with_sffp
        ))
        lines.append("")
    if rank_no_sffp_pvalues is not None and rank_no_sffp is not None:
        lines.append("RankNoSFFPPValues")
        lines.extend(format_rank_pvalue_table(
            rank_no_sffp_pvalues, rank_no_sffp.algos, rank_rows_no_sffp
        ))
        lines.append("")
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

    sum3, sffp, pr90, top3 = group_rows(args.inputs)
    if not sum3 and not sffp and not pr90 and not top3:
        print("error: no summary statistics found in any input file", file=sys.stderr)
        return 1

    base = Path(args.name)
    tsv_path = base.with_suffix(".tsv")
    txt_path = base.with_suffix(".txt")
    write_tsv(tsv_path, sum3, sffp, pr90, top3)
    write_txt(txt_path, sum3, sffp, pr90, top3)
    print(f"wrote {tsv_path} and {txt_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
