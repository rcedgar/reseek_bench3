#!/usr/bin/env python3
"""Infer consensus motif I/II/III columns in the CDN strumm seed MSA.

Locates seqI/seqII/seqIII peptides from data/cdn_motifs.tsv in the
full-length (orig) MSA, maps retained residues onto seed MSA columns
via the msa_prep map, and writes one consensus start/end column range
per motif. Labels are joined by PDB id (token before '-' or '_').

Usage:
  python py/cdn_map_motifs.py
  python py/cdn_map_motifs.py -o strumm/cdn_strumm_seed_motif_cols.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

DefaultMotifs = os.path.join(PalmDir, "data", "cdn_motifs.tsv")
DefaultOrig = os.path.join(PalmDir, "strumm", "cdn.muscle3d.afa")
DefaultMap = os.path.join(PalmDir, "strumm", "cdn.map")
DefaultSeed = os.path.join(PalmDir, "strumm", "cdn_strumm_seed.afa")
DefaultOut = os.path.join(PalmDir, "strumm", "cdn_strumm_seed_motif_cols.tsv")

MOTIF_NAMES = ("I", "II", "III")
MOTIF_TSV_COLS = {
	"I": 1,
	"II": 2,
	"III": 3,
}


def trunc_label(label: str) -> str:
	"""Match muscle default trunc_label: cut at first whitespace, '|', or '/'."""
	n = len(label)
	for i, c in enumerate(label):
		if c in " \t\r\n|/":
			n = i
			break
	return label[:n]


def pdb_id(label: str) -> str:
	"""PDB id = first token before '-' or '_' (after stripping -cdn/-decoy)."""
	s = trunc_label(label)
	for suf in ("-cdn", "-decoy"):
		if s.endswith(suf):
			s = s[: -len(suf)]
			break
	if "_" in s:
		s = s.split("_")[0]
	elif "-" in s:
		s = s.split("-")[0]
	return s.upper()


def is_gap(c: str) -> bool:
	return c in "-."


def motif_or_none(s: str) -> Optional[str]:
	t = s.strip().upper()
	if t == "" or t == ".":
		return None
	return t


def read_fasta(path: str) -> Dict[str, str]:
	"""Read FASTA keyed by PDB id (uppercase)."""
	seqs: Dict[str, str] = {}
	label: Optional[str] = None
	parts: List[str] = []

	def flush() -> None:
		nonlocal label, parts
		if label is None:
			return
		key = pdb_id(label)
		seq = "".join(parts)
		if key in seqs:
			raise SystemExit("%s: duplicate PDB id %s (%s)" % (path, key, label))
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
		parts.append(line)
	flush()
	if len(seqs) == 0:
		raise SystemExit("%s: empty FASTA" % path)
	return seqs


def parse_region_fields(fields: List[str]) -> Tuple[str, Dict[str, str]]:
	kind = fields[0]
	i = 1
	if i < len(fields) and fields[i].isdigit():
		nxt = fields[i + 1] if i + 1 < len(fields) else ""
		if nxt in ("orig_lo", "orig_hi", "seed_lo", "seed_hi", "orig_cols"):
			i += 1
	d: Dict[str, str] = {}
	while i + 1 < len(fields):
		d[fields[i]] = fields[i + 1]
		i += 2
	return kind, d


def parse_msa_prep_map(path: str) -> Tuple[int, List[int]]:
	kept: List[int] = []
	orig_cols: Optional[int] = None
	seed_cols: Optional[int] = None
	saw_header = False
	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if not line:
			continue
		fields = line.split("\t")
		if fields[0] == "msa_prep_map":
			if len(fields) < 2 or fields[1] != "1":
				raise SystemExit("%s: unsupported map version %r" % (path, line))
			saw_header = True
			continue
		if not saw_header:
			raise SystemExit("%s: expected msa_prep_map header, got %r" % (path, line))
		if fields[0] == "orig_cols":
			orig_cols = int(fields[1])
			continue
		if fields[0] == "seed_cols":
			seed_cols = int(fields[1])
			continue
		kind, d = parse_region_fields(fields)
		if kind == "discard":
			continue
		if kind == "block":
			lo = int(d["orig_lo"])
			hi = int(d["orig_hi"])
			if hi < lo:
				raise SystemExit("%s: block orig_hi < orig_lo" % path)
			kept.extend(range(lo, hi + 1))
			continue
		if kind == "spacer":
			cols = d.get("orig_cols", "")
			if cols == "":
				raise SystemExit("%s: spacer missing orig_cols" % path)
			kept.extend(int(x) for x in cols.split(",") if x != "")
			continue
		raise SystemExit("%s: unknown map record %r" % (path, fields[0]))

	if not saw_header:
		raise SystemExit("%s: missing msa_prep_map header" % path)
	if orig_cols is None or seed_cols is None:
		raise SystemExit("%s: missing orig_cols/seed_cols" % path)
	if len(kept) != seed_cols:
		raise SystemExit(
			"%s: listed %d orig columns, seed_cols=%d" %
			(path, len(kept), seed_cols)
		)
	return orig_cols, kept


def find_unique_start(seq: str, motif: str) -> Tuple[Optional[int], int]:
	hits: List[int] = []
	start = 0
	while True:
		i = seq.find(motif, start)
		if i < 0:
			break
		hits.append(i)
		start = i + 1
	if len(hits) == 0:
		return None, 0
	if len(hits) > 1:
		return None, len(hits)
	return hits[0], 1


def ungapped(row: str) -> str:
	return "".join(c for c in row if not is_gap(c)).upper()


def residue_to_orig_col(row: str) -> List[int]:
	"""Map ungapped residue index -> orig MSA column (0-based)."""
	out: List[int] = []
	for col, c in enumerate(row):
		if is_gap(c):
			continue
		out.append(col)
	return out


def median_int(values: Sequence[int]) -> int:
	vals = sorted(values)
	n = len(vals)
	if n == 0:
		raise ValueError("median of empty sequence")
	mid = n // 2
	if n % 2 == 1:
		return vals[mid]
	return (vals[mid - 1] + vals[mid]) // 2


def read_motif_tsv(path: str) -> List[Tuple[str, Dict[str, Optional[str]]]]:
	rows: List[Tuple[str, Dict[str, Optional[str]]]] = []
	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if not line:
			continue
		fields = line.split("\t")
		if fields[0] in ("Label", "label"):
			continue
		if len(fields) < 4:
			raise SystemExit("%s: bad motif row: %s" % (path, line[:80]))
		label = fields[0]
		motifs: Dict[str, Optional[str]] = {}
		for name in MOTIF_NAMES:
			motifs[name] = motif_or_none(fields[MOTIF_TSV_COLS[name]])
		rows.append((label, motifs))
	if not rows:
		raise SystemExit("%s: no motif rows" % path)
	return rows


def map_motif_to_seed_cols(
	row: str,
	motif: str,
	orig_to_seed: Dict[int, int],
) -> Optional[List[int]]:
	u = ungapped(row)
	start, n_hits = find_unique_start(u, motif)
	if n_hits == 0:
		return None
	if n_hits > 1:
		return None
	assert start is not None
	res_cols = residue_to_orig_col(row)
	seed_cols: List[int] = []
	for offset in range(len(motif)):
		orig_col = res_cols[start + offset]
		seed_col = orig_to_seed.get(orig_col)
		if seed_col is not None:
			seed_cols.append(seed_col)
	return seed_cols


def infer_consensus(
	motif_rows: List[Tuple[str, Dict[str, Optional[str]]]],
	orig_seqs: Dict[str, str],
	orig_to_seed: Dict[int, int],
) -> Tuple[Dict[str, Optional[Tuple[int, int]]], Dict[str, int]]:
	"""Return 1-based inclusive (start, end) per motif, plus hit counts."""
	per_motif_spans: Dict[str, List[Tuple[int, int]]] = {n: [] for n in MOTIF_NAMES}
	hit_counts: Dict[str, int] = {n: 0 for n in MOTIF_NAMES}

	for label, motifs in motif_rows:
		key = pdb_id(label)
		if key not in orig_seqs:
			print(
				"warning: PDB id %r (%s) not in orig MSA; skipping" %
				(key, label),
				file=sys.stderr,
			)
			continue
		row = orig_seqs[key]
		for name in MOTIF_NAMES:
			motif = motifs[name]
			if motif is None:
				continue
			u = ungapped(row)
			start, n_hits = find_unique_start(u, motif)
			if n_hits == 0:
				print(
					"warning: motif %s not found in %s: %s" %
					(name, label, motif),
					file=sys.stderr,
				)
				continue
			if n_hits > 1:
				print(
					"warning: motif %s found %d times in %s; skipping" %
					(name, n_hits, label),
					file=sys.stderr,
				)
				continue
			cols = map_motif_to_seed_cols(row, motif, orig_to_seed)
			if cols is None or len(cols) == 0:
				print(
					"warning: motif %s in %s maps to no seed columns" %
					(name, label),
					file=sys.stderr,
				)
				continue
			lo = min(cols)
			hi = max(cols)
			per_motif_spans[name].append((lo, hi))
			hit_counts[name] += 1

	consensus: Dict[str, Optional[Tuple[int, int]]] = {}
	for name in MOTIF_NAMES:
		spans = per_motif_spans[name]
		if not spans:
			consensus[name] = None
			print(
				"warning: no mapped hits for motif %s" % name,
				file=sys.stderr,
			)
			continue
		start0 = median_int([lo for lo, _ in spans])
		end0 = median_int([hi for _, hi in spans])
		if end0 < start0:
			start0, end0 = end0, start0
		# 1-based inclusive
		consensus[name] = (start0 + 1, end0 + 1)
	return consensus, hit_counts


def main(argv: Optional[List[str]] = None) -> int:
	ap = argparse.ArgumentParser(
		description="Infer consensus I/II/III motif columns in the CDN strumm seed MSA"
	)
	ap.add_argument("--motifs", default=DefaultMotifs, help="CDN motif TSV")
	ap.add_argument("--orig", default=DefaultOrig, help="full-length orig MSA")
	ap.add_argument("--map", dest="map_path", default=DefaultMap, help="msa_prep map")
	ap.add_argument("--seed", default=DefaultSeed, help="seed MSA (width check)")
	ap.add_argument("-o", "--output", default=DefaultOut, help="output TSV")
	args = ap.parse_args(argv)

	for path, label in (
		(args.motifs, "motifs"),
		(args.orig, "orig MSA"),
		(args.map_path, "map"),
		(args.seed, "seed MSA"),
	):
		if not os.path.isfile(path):
			raise SystemExit("%s not found: %s" % (label, path))

	orig_cols, kept = parse_msa_prep_map(args.map_path)
	orig_to_seed = {orig_col: seed_col for seed_col, orig_col in enumerate(kept)}

	orig_seqs = read_fasta(args.orig)
	widths = {len(s) for s in orig_seqs.values()}
	if len(widths) != 1:
		raise SystemExit("%s: uneven sequence lengths" % args.orig)
	width = widths.pop()
	if width != orig_cols:
		raise SystemExit(
			"%s: MSA width %d != map orig_cols %d" %
			(args.orig, width, orig_cols)
		)

	seed_seqs = read_fasta(args.seed)
	seed_widths = {len(s) for s in seed_seqs.values()}
	if len(seed_widths) != 1:
		raise SystemExit("%s: uneven sequence lengths" % args.seed)
	seed_width = seed_widths.pop()
	if seed_width != len(kept):
		raise SystemExit(
			"%s: seed width %d != map seed_cols %d" %
			(args.seed, seed_width, len(kept))
		)

	motif_rows = read_motif_tsv(args.motifs)
	consensus, hit_counts = infer_consensus(motif_rows, orig_seqs, orig_to_seed)

	out_dir = os.path.dirname(os.path.abspath(args.output))
	if out_dir and not os.path.isdir(out_dir):
		os.makedirs(out_dir)

	with open(args.output, "w", encoding="utf-8", newline="\n") as fout:
		fout.write("motif\tstart\tend\n")
		for name in MOTIF_NAMES:
			span = consensus[name]
			if span is None:
				fout.write("%s\t.\t.\n" % name)
			else:
				fout.write("%s\t%d\t%d\n" % (name, span[0], span[1]))

	for name in MOTIF_NAMES:
		span = consensus[name]
		if span is None:
			print(
				"%s: hits=%d consensus=." % (name, hit_counts[name]),
				file=sys.stderr,
			)
		else:
			print(
				"%s: hits=%d consensus=%d-%d" %
				(name, hit_counts[name], span[0], span[1]),
				file=sys.stderr,
			)
	print("wrote %s" % args.output, file=sys.stderr)
	return 0


if __name__ == "__main__":
	sys.exit(main())
