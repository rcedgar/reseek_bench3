#!/usr/bin/env python3
"""Classify queries as cdn, decoy, or ambiguous from strumm a3m.

Reads profile-search a3m (match-state columns = uppercase + gaps) and
consensus motif column ranges from cdn_map_motifs.py output. Motif
extraction allows internal deletions; catalytic checks use fixed MSA
columns (I positions 7–8 = GS; II 3,5 = D,D; III 3 = D).

Usage:
  python py/cdn_classify.py cdn_self_test/search.a3m
  python py/cdn_classify.py search.a3m -o out.tsv --cdnfasta cdn.fa
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

DefaultMotifCols = os.path.join(
	PalmDir, "strumm", "cdn_strumm_seed_motif_cols.tsv"
)

MOTIF_NAMES = ["I", "II", "III"]

# 1-based catalytic residue positions inside motifs I, II, III
# Format GS/D-D/D from I[7,8], II[3,5], III[3]
CATALYTIC_CHECKS: List[Tuple[str, int, str]] = [
	("I", 7, "G"),
	("I", 8, "S"),
	("II", 3, "D"),
	("II", 5, "D"),
	("III", 3, "D"),
]

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


def letter_at_col(row: str, lo1: int, pos1: int) -> Optional[str]:
	"""Letter at fixed MSA column lo1+pos1-1, or None if gap/OOB."""
	col = lo1 - 1 + (pos1 - 1)
	if col < 0 or col >= len(row):
		return None
	ch = row[col]
	if is_gap(ch):
		return None
	return ch


def catalytic_at_col(row: str, lo1: int, pos1: int, letter: str) -> bool:
	"""True if fixed MSA column lo1+pos1-1 equals letter."""
	ch = letter_at_col(row, lo1, pos1)
	return ch == letter


def catalytic_present(row: str, lo1: int, pos1: int) -> bool:
	"""True if catalytic column is non-gap."""
	return letter_at_col(row, lo1, pos1) is not None


def extract_motifs(
	row: str, cols: Dict[str, Tuple[int, int]]
) -> Dict[str, Optional[str]]:
	"""Extract ungapped motif strings from match row."""
	sequences: Dict[str, Optional[str]] = {}
	for name in MOTIF_NAMES:
		lo, hi = cols[name]
		window = slice_window(row, lo, hi)
		letters = ungap(window)
		sequences[name] = letters if letters else None
	return sequences


def fixed_catalytic_string(
	row: str, cols: Dict[str, Tuple[int, int]]
) -> Optional[str]:
	"""Build GS/D-D/D-style string from fixed columns, or None if any gap."""
	chars: List[Optional[str]] = []
	for name, pos1, _ in CATALYTIC_CHECKS:
		lo, _ = cols[name]
		chars.append(letter_at_col(row, lo, pos1))
	if any(c is None for c in chars):
		return None
	# I7 I8 / II3 - II5 / III3
	return "%s%s/%s-%s/%s" % (chars[0], chars[1], chars[2], chars[3], chars[4])


def fixed_catalytic_ok(row: str, cols: Dict[str, Tuple[int, int]]) -> bool:
	for name, pos1, letter in CATALYTIC_CHECKS:
		lo, _ = cols[name]
		if not catalytic_at_col(row, lo, pos1, letter):
			return False
	return True


def fixed_catalytic_present(row: str, cols: Dict[str, Tuple[int, int]]) -> bool:
	for name, pos1, _ in CATALYTIC_CHECKS:
		lo, _ = cols[name]
		if not catalytic_present(row, lo, pos1):
			return False
	return True


def classify_row(
	row: str, cols: Dict[str, Tuple[int, int]]
) -> Tuple[str, Dict[str, Optional[str]], Optional[str]]:
	sequences = extract_motifs(row, cols)
	cat = fixed_catalytic_string(row, cols)

	if fixed_catalytic_ok(row, cols):
		family = "cdn"
	elif fixed_catalytic_present(row, cols):
		family = "decoy"
	else:
		family = "ambiguous"

	return family, sequences, cat


def motif_or_dot(seq: Optional[str]) -> str:
	if seq is None:
		return "."
	return seq


def true_family(label: str) -> Optional[str]:
	"""Gold label: 'Cdn' in name => cdn; else decoy for self-test set."""
	if "Cdn" in label:
		return "cdn"
	return "decoy"


def ungapped_upper(raw: str) -> str:
	"""Full query sequence: strip gaps, uppercase letters."""
	return "".join(c.upper() for c in raw if not is_gap(c))


def write_cdn_fasta(out_fh, rows: List[dict]) -> None:
	"""Write sequences classified as cdn as ungapped uppercase FASTA."""
	for row in rows:
		if row["family"] != "cdn":
			continue
		motifs = row["motifs"]
		annots = [
			"I=%s" % motif_or_dot(motifs["I"]),
			"II=%s" % motif_or_dot(motifs["II"]),
			"III=%s" % motif_or_dot(motifs["III"]),
			"cat=%s" % motif_or_dot(row["catalytic"]),
		]
		seq = ungapped_upper(row["raw"])
		out_fh.write(">%s %s\n" % (row["query"], " ".join(annots)))
		for i in range(0, len(seq), 80):
			out_fh.write(seq[i : i + 80] + "\n")


def write_tsv(out_fh, rows: List[dict]) -> None:
	out_fh.write("query\tfamily\tI\tII\tIII\tcatalytic\n")
	for row in rows:
		motifs = row["motifs"]
		parts = [
			row["query"],
			row["family"],
			motif_or_dot(motifs["I"]),
			motif_or_dot(motifs["II"]),
			motif_or_dot(motifs["III"]),
			motif_or_dot(row["catalytic"]),
		]
		out_fh.write("\t".join(parts) + "\n")


def smoke_summary(fn, rows: List[dict]) -> None:
	frep = open(fn, "w")
	n_cdn = tp = fn_count = 0
	n_decoy = fp = 0
	for row in rows:
		true = true_family(row["query"])
		pred = row["family"]
		if true == "cdn":
			n_cdn += 1
			if pred == "cdn":
				tp += 1
			else:
				fn_count += 1
		elif true == "decoy":
			n_decoy += 1
			if pred == "cdn":
				fp += 1
	if n_cdn == 0 and n_decoy == 0:
		frep.close()
		return
	print("Self-test report:", file=frep)
	if n_cdn:
		sens = 100.0 * tp / n_cdn
		print("CDN queries: %d" % n_cdn, file=frep)
		print("True positives:  %d" % tp, file=frep)
		print("False negatives: %d" % fn_count, file=frep)
		print("Sensitivity:     %.1f%%" % sens, file=frep)
	if n_decoy:
		print("Decoy queries: %d" % n_decoy, file=frep)
		print("Decoys classified as CDN (FP): %d" % fp, file=frep)
	frep.close()


def main(argv: Optional[List[str]] = None) -> int:
	ap = argparse.ArgumentParser(
		description="Classify queries from CDN strumm profile-search a3m"
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
		"-cdnfasta",
		"--cdnfasta",
		dest="cdnfasta",
		help="FASTA of sequences classified as cdn (uppercase, gaps stripped)",
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
		family, motifs, cat = classify_row(mrow, cols)
		rows.append(
			{
				"query": query,
				"family": family,
				"motifs": motifs,
				"catalytic": cat,
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

	if args.cdnfasta:
		out_dir = os.path.dirname(os.path.abspath(args.cdnfasta))
		if out_dir and not os.path.isdir(out_dir):
			os.makedirs(out_dir)
		with open(args.cdnfasta, "w", encoding="utf-8", newline="\n") as fout:
			write_cdn_fasta(fout, rows)
		n_cdn = sum(1 for r in rows if r["family"] == "cdn")
		print(
			"wrote %d cdn sequences to %s" % (n_cdn, args.cdnfasta),
			file=sys.stderr,
		)

	if args.selfreport:
		smoke_summary(args.selfreport, rows)
	return 0


if __name__ == "__main__":
	sys.exit(main())
