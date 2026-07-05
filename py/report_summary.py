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

MetricAgreement quantifies rank concordance between metrics: Spearman correlation
(cell-pooled and per-algo median ranks) plus SFFP discordance vs Sum3/PR90.

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

METRIC_LABELS = ("Sum3", "SFFP", "PR90", "Top3")


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


def build_rank_table(
    sections: Sequence[Tuple[str, Dict[GroupKey, List[Row]]]],
) -> Optional[RankTable]:
    row_labels: List[str] = []
    rank_rows: List[Dict[str, int]] = []
    for section, groups in sections:
        for truth, reference in sorted_group_keys(groups):
            rows = groups[(truth, reference)]
            row_labels.append(section_subhead(section, truth, reference))
            rank_rows.append(competition_ranks(rows))

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


def _rank_values(values: Sequence[float]) -> List[float]:
    """Average ranks for ties (1-based)."""
    ranks = [0.0] * len(values)
    indexed = sorted(enumerate(values), key=lambda t: t[1])
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = avg_rank
        i = j
    return ranks


def spearman_rho(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    rx = _rank_values(xs)
    ry = _rank_values(ys)
    mx = statistics.mean(rx)
    my = statistics.mean(ry)
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(len(xs)))
    den_x = math.sqrt(sum((r - mx) ** 2 for r in rx))
    den_y = math.sqrt(sum((r - my) ** 2 for r in ry))
    if den_x == 0.0 or den_y == 0.0:
        return None
    return num / (den_x * den_y)


def collect_cell_rank_pairs(
    groups_a: Dict[GroupKey, List[Row]],
    groups_b: Dict[GroupKey, List[Row]],
) -> Tuple[List[float], List[float], int]:
    """Pool (rank, rank) pairs across benchmark cells where both metrics exist."""
    xs: List[float] = []
    ys: List[float] = []
    n_cells = 0
    for key in sorted(set(groups_a) & set(groups_b)):
        ranks_a = competition_ranks(groups_a[key])
        ranks_b = competition_ranks(groups_b[key])
        n_cells += 1
        for algo in sorted(set(ranks_a) & set(ranks_b)):
            xs.append(float(ranks_a[algo]))
            ys.append(float(ranks_b[algo]))
    return xs, ys, n_cells


def per_algo_median_ranks(groups: Dict[GroupKey, List[Row]]) -> Dict[str, float]:
    """Median competition rank per algo across all benchmark cells for one metric."""
    by_algo: Dict[str, List[int]] = {}
    for rows in groups.values():
        for algo, rank in competition_ranks(rows).items():
            by_algo.setdefault(algo, []).append(rank)
    return {algo: statistics.median(ranks) for algo, ranks in by_algo.items()}


@dataclass(frozen=True)
class MetricPairStats:
    metric_a: str
    metric_b: str
    n_cells: int
    n_pairs: int
    spearman: Optional[float]
    mean_abs_rank_diff: Optional[float]


@dataclass(frozen=True)
class SffpDiscordanceRow:
    algo: str
    median_delta: float
    pct_sffp_worse: float
    n_cells: int


@dataclass(frozen=True)
class MetricAgreement:
    cell_pairs: Tuple[MetricPairStats, ...]
    algo_pairs: Tuple[MetricPairStats, ...]
    sffp_discordance: Tuple[SffpDiscordanceRow, ...]


def build_metric_agreement(
    sum3: Dict[GroupKey, List[Row]],
    sffp: Dict[GroupKey, List[Row]],
    pr90: Dict[GroupKey, List[Row]],
    top3: Dict[GroupKey, List[Row]],
) -> Optional[MetricAgreement]:
    groups = {
        "Sum3": sum3,
        "SFFP": sffp,
        "PR90": pr90,
        "Top3": top3,
    }
    active = [name for name in METRIC_LABELS if groups[name]]
    if len(active) < 2:
        return None

    cell_pairs: List[MetricPairStats] = []
    algo_pairs: List[MetricPairStats] = []
    for i, name_a in enumerate(active):
        for name_b in active[i + 1:]:
            xs, ys, n_cells = collect_cell_rank_pairs(
                groups[name_a], groups[name_b]
            )
            rho = spearman_rho(xs, ys) if xs else None
            mad = (
                statistics.mean(abs(xs[j] - ys[j]) for j in range(len(xs)))
                if xs
                else None
            )
            cell_pairs.append(MetricPairStats(
                metric_a=name_a,
                metric_b=name_b,
                n_cells=n_cells,
                n_pairs=len(xs),
                spearman=rho,
                mean_abs_rank_diff=mad,
            ))

            med_a = per_algo_median_ranks(groups[name_a])
            med_b = per_algo_median_ranks(groups[name_b])
            algos = sorted(set(med_a) & set(med_b))
            ax = [med_a[a] for a in algos]
            ay = [med_b[a] for a in algos]
            algo_pairs.append(MetricPairStats(
                metric_a=name_a,
                metric_b=name_b,
                n_cells=len(algos),
                n_pairs=len(algos),
                spearman=spearman_rho(ax, ay),
                mean_abs_rank_diff=(
                    statistics.mean(abs(ax[j] - ay[j]) for j in range(len(ax)))
                    if ax
                    else None
                ),
            ))

    sffp_rows: List[SffpDiscordanceRow] = []
    if sffp:
        deltas_by_algo: Dict[str, List[float]] = {}
        worse_by_algo: Dict[str, List[bool]] = {}
        for key in sorted(sffp):
            ranks_sffp = competition_ranks(sffp[key])
            ranks_sum3 = (
                competition_ranks(sum3[key]) if key in sum3 else {}
            )
            ranks_pr90 = (
                competition_ranks(pr90[key]) if key in pr90 else {}
            )
            for algo, rank_sffp in ranks_sffp.items():
                others: List[int] = []
                if algo in ranks_sum3:
                    others.append(ranks_sum3[algo])
                if algo in ranks_pr90:
                    others.append(ranks_pr90[algo])
                if not others:
                    continue
                consensus = statistics.median(others)
                delta = float(rank_sffp) - consensus
                deltas_by_algo.setdefault(algo, []).append(delta)
                worse_by_algo.setdefault(algo, []).append(delta > 0)

        for algo in sorted(deltas_by_algo):
            deltas = deltas_by_algo[algo]
            worse = worse_by_algo[algo]
            sffp_rows.append(SffpDiscordanceRow(
                algo=algo,
                median_delta=statistics.median(deltas),
                pct_sffp_worse=sum(worse) / len(worse),
                n_cells=len(deltas),
            ))

    return MetricAgreement(
        cell_pairs=tuple(cell_pairs),
        algo_pairs=tuple(algo_pairs),
        sffp_discordance=tuple(sffp_rows),
    )


def _fmt_rho(rho: Optional[float]) -> str:
    if rho is None:
        return "."
    return f"{rho:.3f}"


def _fmt_float(value: Optional[float]) -> str:
    if value is None:
        return "."
    return f"{value:.3f}"


def write_metric_agreement_tsv(f, agreement: MetricAgreement) -> None:
    f.write("# RankCorrCells\n")
    f.write("\t".join((
        "metric_a", "metric_b", "n_cells", "n_pairs",
        "spearman", "mean_abs_rank_diff",
    )) + "\n")
    for row in agreement.cell_pairs:
        f.write("\t".join((
            row.metric_a,
            row.metric_b,
            str(row.n_cells),
            str(row.n_pairs),
            _fmt_rho(row.spearman),
            _fmt_float(row.mean_abs_rank_diff),
        )) + "\n")
    f.write("\n")

    f.write("# RankCorrAlgos\n")
    f.write("\t".join((
        "metric_a", "metric_b", "n_algos", "spearman", "mean_abs_rank_diff",
    )) + "\n")
    for row in agreement.algo_pairs:
        f.write("\t".join((
            row.metric_a,
            row.metric_b,
            str(row.n_pairs),
            _fmt_rho(row.spearman),
            _fmt_float(row.mean_abs_rank_diff),
        )) + "\n")
    f.write("\n")

    if agreement.sffp_discordance:
        f.write("# SFFPDiscordance\n")
        f.write("\t".join((
            "algo", "median_delta", "pct_sffp_worse", "n_cells",
        )) + "\n")
        for row in agreement.sffp_discordance:
            f.write("\t".join((
                row.algo,
                _fmt_float(row.median_delta),
                f"{row.pct_sffp_worse:.3f}",
                str(row.n_cells),
            )) + "\n")
        f.write("\n")


def format_metric_agreement_txt(agreement: MetricAgreement) -> List[str]:
    lines: List[str] = []

    lines.append("RankCorrCells")
    cols = ("metric_a", "metric_b", "n_cells", "n_pairs", "spearman", "mad")
    widths = [len(c) for c in cols]
    rows = [
        (
            p.metric_a,
            p.metric_b,
            str(p.n_cells),
            str(p.n_pairs),
            _fmt_rho(p.spearman),
            _fmt_float(p.mean_abs_rank_diff),
        )
        for p in agreement.cell_pairs
    ]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    lines.append("  ".join(c.ljust(widths[i]) for i, c in enumerate(cols)))
    for row in rows:
        lines.append("  ".join(
            (row[0].ljust(widths[0]), row[1].ljust(widths[1]),
             row[2].rjust(widths[2]), row[3].rjust(widths[3]),
             row[4].rjust(widths[4]), row[5].rjust(widths[5]))
        ))
    lines.append("")

    lines.append("RankCorrAlgos")
    cols = ("metric_a", "metric_b", "n_algos", "spearman", "mad")
    widths = [len(c) for c in cols]
    rows = [
        (
            p.metric_a,
            p.metric_b,
            str(p.n_pairs),
            _fmt_rho(p.spearman),
            _fmt_float(p.mean_abs_rank_diff),
        )
        for p in agreement.algo_pairs
    ]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))
    lines.append("  ".join(c.ljust(widths[i]) for i, c in enumerate(cols)))
    for row in rows:
        lines.append("  ".join(
            (row[0].ljust(widths[0]), row[1].ljust(widths[1]),
             row[2].rjust(widths[2]), row[3].rjust(widths[3]),
             row[4].rjust(widths[4]))
        ))
    lines.append("")

    if agreement.sffp_discordance:
        lines.append("SFFPDiscordance")
        cols = ("algo", "median_delta", "pct_worse", "n_cells")
        widths = [len(c) for c in cols]
        rows = [
            (
                r.algo,
                _fmt_float(r.median_delta),
                f"{r.pct_sffp_worse:.3f}",
                str(r.n_cells),
            )
            for r in agreement.sffp_discordance
        ]
        for row in rows:
            for i, cell in enumerate(row):
                widths[i] = max(widths[i], len(cell))
        lines.append("  ".join(c.ljust(widths[i]) for i, c in enumerate(cols)))
        for row in rows:
            lines.append("  ".join(
                (row[0].ljust(widths[0]), row[1].rjust(widths[1]),
                 row[2].rjust(widths[2]), row[3].rjust(widths[3]))
            ))
        lines.append("")

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
    agreement = build_metric_agreement(sum3, sffp, pr90, top3)
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
        if agreement is not None:
            write_metric_agreement_tsv(f, agreement)


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
    rank_table = build_rank_table(
        rank_sections(sum3, sffp, pr90, top3, include_sffp=True)
    )
    rank_no_sffp = build_rank_table(
        rank_sections(sum3, sffp, pr90, top3, include_sffp=False)
    )
    agreement = build_metric_agreement(sum3, sffp, pr90, top3)
    if rank_table is not None:
        lines.append("Rank")
        lines.extend(format_rank_table(rank_table))
        lines.append("")
    if rank_no_sffp is not None:
        lines.append("RankNoSFFP")
        lines.extend(format_rank_table(rank_no_sffp))
        lines.append("")
    if agreement is not None:
        lines.extend(format_metric_agreement_txt(agreement))
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
