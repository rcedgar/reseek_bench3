#!/usr/bin/env python3
"""Classify queries as viral RdRp, decoy, RT, or ambiguous from strumm a3m.

Reads profile-search a3m (match-state columns = uppercase + gaps) and
consensus motif column ranges from rdrp_map_motifs.py output. Motif
extraction allows internal deletions; catalytic checks use fixed MSA
columns (palmfinder A/B/C positions). Gate is ungapped motif-A index 10.
RdRp requires catalytic A/B/C, a known gate, and xDx in xDD or GDN.

Usage:
  python py/rdrp_classify.py rdrp_self_test/search.a3m
  python py/rdrp_classify.py search.a3m -o out.tsv -rdrpfasta rdrp.fa
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

DefaultMotifCols = os.path.join(
	PalmDir, "strumm", "rdrp_abc_strumm_seed_motif_cols.tsv"
)

MOTIF_NAMES = ["F1", "F2", "A", "B", "C", "D", "E"]

# 1-based catalytic residue positions inside motifs A, B, C
CATALYTIC = {
	"A": (6, "D"),
	"B": (2, "G"),
	"C": (6, "D"),
}

# Reverse-transcriptase signatures
RT_XDX = set(["ADD", "VDD", "IDD", "MDD", "TDD"])
RT_XDX1 = set(["L", "V", "I"])
RT_GATE = set(["Y", "F", "W"])

# RdRp catalytic triad motif: any xDD, plus GDN
RDRP_GDN = "GDN"


def trunc_label(label: str) -> str:
	n = len(label)
	for i, c in enumerate(label):
		if c in " \t\r\n|/":
			n = i
			break
	return label[:n]


def is_gap(c: str) -> bool:
	return c in "-."


def read_a3m(path: str) -> Dict[str, str]:
	"""Read a3m; return label -> raw sequence (upper+lower+gap)."""
	seqs: Dict[str, str] = {}
	label: Optional[str] = None
	parts: List[str] = []

	def flush() -> None:
		nonlocal label, parts
		if label is None:
			return
		key = trunc_label(label)
		seq = "".join(parts)
		if key in seqs:
			raise SystemExit("%s: duplicate label %s" % (path, label))
		seqs[key] = seq
		parts = []

	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if line.startswith(">"):
			flush()
			label = line[1:]
			parts = []
			continue
		if label is None:
			continue
		parts.append(line.strip())
	flush()
	if not seqs:
		raise SystemExit("%s: empty a3m" % path)
	return seqs


def match_row(raw: str) -> str:
	"""Keep match-state columns: uppercase letters and gaps."""
	return "".join(c for c in raw if c == "-" or c.isupper())


def load_motif_cols(path: str) -> Dict[str, Tuple[int, int]]:
	"""Load motif -> (start, end) 1-based inclusive column range."""
	cols: Dict[str, Tuple[int, int]] = {}
	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if not line or line.startswith("motif"):
			continue
		fields = line.split("\t")
		if len(fields) != 3:
			raise SystemExit("%s: bad line: %s" % (path, line[:80]))
		name, lo_s, hi_s = fields
		if lo_s == "." or hi_s == ".":
			raise SystemExit("%s: missing columns for motif %s" % (path, name))
		lo, hi = int(lo_s), int(hi_s)
		if hi < lo:
			raise SystemExit("%s: end < start for motif %s" % (path, name))
		cols[name] = (lo, hi)
	for name in MOTIF_NAMES:
		if name not in cols:
			raise SystemExit("%s: missing motif %s" % (path, name))
	return cols


def slice_window(row: str, lo1: int, hi1: int) -> str:
	"""1-based inclusive column slice."""
	return row[lo1 - 1 : hi1]


def ungap(seq: str) -> str:
	return seq.replace("-", "")


def query_start_in_window(row: str, lo1: int, hi1: int) -> Optional[int]:
	"""1-based ungapped query residue at first non-gap in motif window."""
	before = sum(1 for c in row[: lo1 - 1] if not is_gap(c))
	for col in range(lo1 - 1, hi1):
		ch = row[col]
		if not is_gap(ch):
			return before + 1
	return None


def catalytic_at_col(row: str, lo1: int, pos1: int, letter: str) -> bool:
	"""True if fixed MSA column lo1+pos1-1 equals letter."""
	col = lo1 - 1 + (pos1 - 1)
	if col < 0 or col >= len(row):
		return False
	ch = row[col]
	if is_gap(ch):
		return False
	return ch == letter


def catalytic_present(row: str, lo1: int, pos1: int) -> bool:
	"""True if catalytic column is non-gap."""
	col = lo1 - 1 + (pos1 - 1)
	if col < 0 or col >= len(row):
		return False
	return not is_gap(row[col])


def extract_motifs(
	row: str, cols: Dict[str, Tuple[int, int]]
) -> Tuple[Dict[str, Optional[str]], Dict[str, Optional[int]]]:
	"""Extract ungapped motif strings and query starts from match row."""
	sequences: Dict[str, Optional[str]] = {}
	starts: Dict[str, Optional[int]] = {}
	for name in MOTIF_NAMES:
		lo, hi = cols[name]
		window = slice_window(row, lo, hi)
		letters = ungap(window)
		if letters:
			sequences[name] = letters
			starts[name] = query_start_in_window(row, lo, hi)
		else:
			sequences[name] = None
			starts[name] = None
	return sequences, starts


def motif_gate(sequences: Dict[str, Optional[str]]) -> Optional[str]:
	"""Gate = 0-based index 10 in ungapped motif A (may differ from MSA col)."""
	seq = sequences.get("A")
	if seq is None or len(seq) <= 10:
		return None
	return seq[10]


def fixed_xdx(row: str, cols: Dict[str, Tuple[int, int]]) -> Optional[str]:
	"""xDx = letters at C_start+4..+6 (1-based pos 5-7), or None if any gap."""
	lo, _ = cols["C"]
	chars = []
	for off in (4, 5, 6):
		col = lo - 1 + off
		if col < 0 or col >= len(row):
			return None
		ch = row[col]
		if is_gap(ch):
			return None
		chars.append(ch)
	return "".join(chars)


def is_rdrp_xdx(xDx: Optional[str]) -> bool:
	"""True if xDx is xDD (any residue + DD) or GDN."""
	if xDx is None or len(xDx) != 3:
		return False
	if xDx == RDRP_GDN:
		return True
	return xDx[1:] == "DD"


def fixed_catalytic_ok(row: str, cols: Dict[str, Tuple[int, int]]) -> bool:
	for name, (pos1, letter) in CATALYTIC.items():
		lo, _ = cols[name]
		if not catalytic_at_col(row, lo, pos1, letter):
			return False
	return True


def fixed_abc_catalytic_present(row: str, cols: Dict[str, Tuple[int, int]]) -> bool:
	for name, (pos1, _) in CATALYTIC.items():
		lo, _ = cols[name]
		if not catalytic_present(row, lo, pos1):
			return False
	return True


def format_motifs_field(
	sequences: Dict[str, Optional[str]], starts: Dict[str, Optional[int]]
) -> str:
	found_indiv = []
	for name in MOTIF_NAMES:
		if sequences[name] is not None and starts[name] is not None:
			found_indiv.append((starts[name], name))
	found_indiv.sort()
	found_names = [name for _, name in found_indiv]

	f1_found = sequences["F1"] is not None and starts["F1"] is not None
	f2_found = sequences["F2"] is not None and starts["F2"] is not None
	if not f1_found and not f2_found:
		f_char = "."
		f_entry = (1, 0, 0, ".")
	else:
		f_start = min(
			starts[n]
			for n, ok in (("F1", f1_found), ("F2", f2_found))
			if ok
		)
		if f1_found and f2_found:
			f_pos = [j for j, n in enumerate(found_names) if n in ("F1", "F2")]
			f_char = "F" if f_pos[1] == f_pos[0] + 1 else "@"
		else:
			f_char = "F"
		f_entry = (0, f_start, 0, f_char)

	entries = [f_entry]
	for i, name in enumerate(["A", "B", "C", "D", "E"], start=1):
		if sequences[name] is not None and starts[name] is not None:
			entries.append((0, starts[name], i, name))
		else:
			entries.append((1, 0, i, "."))
	entries.sort()
	return "".join(ch for _, _, _, ch in entries)


def is_cab_order(starts: Dict[str, Optional[int]]) -> bool:
	for name in ("A", "B", "C"):
		if starts.get(name) is None:
			return False
	return starts["C"] < starts["A"] < starts["B"]


def apply_cab_suppress_de(
	sequences: Dict[str, Optional[str]], starts: Dict[str, Optional[int]]
) -> bool:
	if not is_cab_order(starts):
		return False
	sequences["D"] = None
	sequences["E"] = None
	starts["D"] = None
	starts["E"] = None
	return True


def classify_row(
	row: str, cols: Dict[str, Tuple[int, int]]
) -> Tuple[str, Dict[str, Optional[str]], Dict[str, Optional[int]], Optional[str], Optional[str]]:
	sequences, starts = extract_motifs(row, cols)
	apply_cab_suppress_de(sequences, starts)

	gate = motif_gate(sequences)
	xDx = fixed_xdx(row, cols)

	if gate in RT_GATE:
		family = "RT"
	elif xDx in RT_XDX:
		family = "RT"
	elif xDx is not None and xDx[0] in RT_XDX1:
		family = "RT"
	elif (
		gate is not None
		and fixed_catalytic_ok(row, cols)
		and is_rdrp_xdx(xDx)
	):
		family = "rdrp"
	elif fixed_abc_catalytic_present(row, cols):
		family = "decoy"
	else:
		family = "ambiguous"

	return family, sequences, starts, xDx, gate


def motif_or_dot(seq: Optional[str]) -> str:
	if seq is None:
		return "."
	return seq


def true_family(label: str) -> Optional[str]:
	if label.endswith("-rdrp"):
		return "rdrp"
	if label.endswith("-decoy"):
		return "decoy"
	return None


def ungapped_upper(raw: str) -> str:
	"""Full query sequence: strip gaps, uppercase letters."""
	return "".join(c.upper() for c in raw if not is_gap(c))


def write_rdrp_fasta(out_fh, rows: List[dict]) -> None:
	"""Write sequences classified as rdrp as ungapped uppercase FASTA."""
	for row in rows:
		if row["family"] != "rdrp":
			continue
		motifs = row["motifs"]
		annots = [
			"A=%s" % motif_or_dot(motifs["A"]),
			"B=%s" % motif_or_dot(motifs["B"]),
			"C=%s" % motif_or_dot(motifs["C"]),
			"GDD=%s" % motif_or_dot(row["xDx"]),
			"gate=%s" % motif_or_dot(row["gate"]),
		]
		seq = ungapped_upper(row["raw"])
		out_fh.write(">%s %s\n" % (row["query"], " ".join(annots)))
		for i in range(0, len(seq), 80):
			out_fh.write(seq[i : i + 80] + "\n")


def write_tsv(out_fh, rows: List[dict]) -> None:
	out_fh.write("query\tfamily\tmotifs\tF1\tF2\tA\tB\tC\tD\tE\tGDD\tgate\n")
	for row in rows:
		motifs = row["motifs"]
		parts = [
			row["query"],
			row["family"],
			row["motifs_str"],
			motif_or_dot(motifs["F1"]),
			motif_or_dot(motifs["F2"]),
			motif_or_dot(motifs["A"]),
			motif_or_dot(motifs["B"]),
			motif_or_dot(motifs["C"]),
			motif_or_dot(motifs["D"]),
			motif_or_dot(motifs["E"]),
			motif_or_dot(row["xDx"]),
			motif_or_dot(row["gate"]),
		]
		out_fh.write("\t".join(parts) + "\n")


def smoke_summary(fn, rows: List[dict]) -> None:
	frep = open(fn, "w")
	n_rdrp = tp = fn = 0
	n_decoy = fp = 0
	for row in rows:
		true = true_family(row["query"])
		pred = row["family"]
		if true == "rdrp":
			n_rdrp += 1
			if pred == "rdrp":
				tp += 1
			else:
				fn += 1
		elif true == "decoy":
			n_decoy += 1
			if pred == "rdrp":
				fp += 1
	if n_rdrp == 0 and n_decoy == 0:
		return
	print("Self-test report:", file=frep)
	if n_rdrp:
		sens = 100.0 * tp / n_rdrp
		print("RdRp queries: %d" % n_rdrp, file=frep)
		print("True positives:  %d" % tp, file=frep)
		print("False negatives: %d" % fn, file=frep)
		print("Sensitivity:     %.1f%%" % sens, file=frep)
	if n_decoy:
		print("Decoy queries: %d" % n_decoy, file=frep)
		print("Decoys classified as RdRp (FP): %d" % fp, file=frep)
	frep.close()


def main(argv: Optional[List[str]] = None) -> int:
	ap = argparse.ArgumentParser(
		description="Classify queries from strumm profile-search a3m"
	)
	ap.add_argument("a3m", help="input a3m from strumm search")
	ap.add_argument(
		"--motif_cols",
		default=DefaultMotifCols,
		help="consensus motif column TSV",
	)
	ap.add_argument(
		"-o",
		"--output",
		default="-",
		help="output TSV (default stdout)",
	)
	ap.add_argument(
		"--selfreport",
		help="self test report",
	)
	ap.add_argument(
		"-rdrpfasta",
		"--rdrpfasta",
		dest="rdrpfasta",
		help="FASTA of sequences classified as rdrp (uppercase, gaps stripped)",
	)

	args = ap.parse_args(argv)

	if not os.path.isfile(args.a3m):
		raise SystemExit("a3m not found: %s" % args.a3m)
	if not os.path.isfile(args.motif_cols):
		raise SystemExit("motif cols not found: %s" % args.motif_cols)

	cols = load_motif_cols(args.motif_cols)

	seqs = read_a3m(args.a3m)
	rows: List[dict] = []
	match_width: Optional[int] = None
	for query in sorted(seqs.keys()):
		raw = seqs[query]
		mrow = match_row(raw)
		if match_width is None:
			match_width = len(mrow)
			for name, (_, hi) in cols.items():
				if hi > match_width:
					raise SystemExit(
						"motif %s end %d > match width %d" %
						(name, hi, match_width)
					)
		elif len(mrow) != match_width:
			raise SystemExit(
				"%s: match row width %d != expected %d" %
				(query, len(mrow), match_width)
			)
		family, motifs, starts, xDx, gate = classify_row(mrow, cols)
		rows.append(
			{
				"query": query,
				"family": family,
				"motifs": motifs,
				"motifs_str": format_motifs_field(motifs, starts),
				"xDx": xDx,
				"gate": gate,
				"raw": raw,
			}
		)

	if args.output == "-":
		write_tsv(sys.stdout, rows)
	else:
		out_dir = os.path.dirname(os.path.abspath(args.output))
		if out_dir and not os.path.isdir(out_dir):
			os.makedirs(out_dir)
		with open(args.output, "w", encoding="utf-8", newline="\n") as fout:
			write_tsv(fout, rows)

	if args.rdrpfasta:
		out_dir = os.path.dirname(os.path.abspath(args.rdrpfasta))
		if out_dir and not os.path.isdir(out_dir):
			os.makedirs(out_dir)
		with open(args.rdrpfasta, "w", encoding="utf-8", newline="\n") as fout:
			write_rdrp_fasta(fout, rows)
		n_rdrp = sum(1 for r in rows if r["family"] == "rdrp")
		print(
			"wrote %d rdrp sequences to %s" % (n_rdrp, args.rdrpfasta),
			file=sys.stderr,
		)

	if args.selfreport:
		smoke_summary(args.selfreport, rows)
	return 0


if __name__ == "__main__":
	sys.exit(main())
