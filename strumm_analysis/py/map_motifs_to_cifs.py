#!/usr/bin/python3

import argparse
import os
import sys
from pathlib import Path

from ref_label_map import BaseOf
from cif_seq import (
	MISSING,
	MakeQuery,
	ParseChainSequences,
	ResolveChain,
	Die,
)

Usage = \
(
"Map RCE motif peptides from a motifs TSV onto mmCIF chain coordinates. "
"Output TSV matches cif_make_sifts_annotations (db=rce_motif)."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("mode", choices=("rdrp", "cdn"),
	help="Motif set: rdrp (F1–E) or cdn (seqI–III)")
AP.add_argument("--motifs", default=None, help="Motifs TSV (mode default)")
AP.add_argument("--cifs", default=None, help="mmCIF directory (mode default)")
AP.add_argument("--output", default="-", help="Output TSV (default stdout)")
AP.add_argument("--chainchar", default="_",
	help="Separator between stem and chain in query id (default '_')")
AP.add_argument("--report", default=None, metavar="TEXTFILE",
	help="Write motif corrections (original → inferred) to TEXTFILE")
Args = AP.parse_args()

DB_NAME = "rce_motif"
CDN_CHAIN = "A"
MIN_HALF = 4

RDRP_SKIP_COLS = {"Label", "label"}
CDN_SEQ_COLS = ("seqI", "seqII", "seqIII")
CDN_SKIP_COLS = {"label", "catalytic"}


def DefaultsFor(mode):
	if mode == "rdrp":
		return (
			os.path.join(PalmDir, "data", "rdrp_motifs.tsv"),
			os.path.join(PalmDir, "data", "rdrp_cif"),
			None,
		)
	return (
		os.path.join(PalmDir, "data", "cdn_motifs.tsv"),
		os.path.join(PalmDir, "data", "cdn_cif"),
		CDN_SEQ_COLS,
	)


def MotifOrNone(s):
	if s is None:
		return None
	t = s.strip().upper()
	if t in MISSING:
		return None
	return t


def PdbId(label):
	s = label.strip()
	for suf in ("-cdn", "-decoy"):
		if s.endswith(suf):
			s = s[: -len(suf)]
			break
	if "_" in s:
		s = s.split("_")[0]
	elif "-" in s:
		s = s.split("-")[0]
	return s


def ParseRdrpLabel(label):
	base = BaseOf(label.strip())
	if "_" not in base:
		Die("rdrp label lacks '_': %s" % label)
	pdb, chain = base.rsplit("_", 1)
	if pdb == "" or chain == "":
		Die("bad rdrp label: %s" % label)
	return pdb.lower(), chain


def ParseCdnLabel(label):
	return PdbId(label).lower()


def FindCifFile(cifs_dir, pdb_id):
	pdb_id = pdb_id.lower()
	direct = Path(cifs_dir) / (pdb_id + ".cif")
	if direct.is_file():
		return direct
	for path in Path(cifs_dir).glob("*.cif"):
		if path.stem.lower() == pdb_id:
			return path
	return None


def _find_all(seq, peptide):
	hits = []
	start = 0
	while True:
		i = seq.find(peptide, start)
		if i < 0:
			break
		hits.append(i)
		start = i + 1
	return hits


def _search_views(chain_residues):
	"""Full sequence and gap-removed (nonstandard omitted) views with index maps."""
	full_seq = "".join(aa for _sid, aa, _gap in chain_residues)
	full_map = list(range(len(chain_residues)))
	gap_seq_chars = []
	gap_map = []
	for i, (_sid, aa, is_gap) in enumerate(chain_residues):
		if is_gap:
			continue
		gap_seq_chars.append(aa)
		gap_map.append(i)
	return (full_seq, full_map), ("".join(gap_seq_chars), gap_map)


def _span_from_hit(chain_residues, index_map, lo, length):
	"""Map a hit on a search string back to full-chain label_seq_id span."""
	full_lo = index_map[lo]
	full_hi = index_map[lo + length - 1]
	start_id = chain_residues[full_lo][0]
	end_id = chain_residues[full_hi][0]
	matched = "".join(
		aa for _sid, aa, _gap in chain_residues[full_lo : full_hi + 1])
	return start_id, end_id, matched


def FindUniqueMotif(chain_residues, peptide):
	"""Exact unique match on the full chain sequence.

	Returns None, ("multi", hits), or ("ok", (start_id, end_id, matched)).
	"""
	if peptide is None or len(peptide) == 0:
		return None
	(full_seq, full_map), _gap = _search_views(chain_residues)
	hits = _find_all(full_seq, peptide)
	if len(hits) == 0:
		return None
	if len(hits) > 1:
		return "multi", hits
	start_id, end_id, matched = _span_from_hit(
		chain_residues, full_map, hits[0], len(peptide))
	return "ok", (start_id, end_id, matched)


def InferByHalves(chain_residues, peptide):
	"""Infer motif by uniquely locating first and second halves.

	Tries the full sequence first, then a gap-removed view where nonstandard
	residues (e.g. CAS) are omitted — matching PyMOL 'remove gaps' behavior.
	Returns ((start_id, end_id, corrected), None) or (None, reason).
	"""
	n = len(peptide)
	if n < 2 * MIN_HALF:
		return None, "too_short"
	left = peptide[: n // 2]
	right = peptide[n // 2 :]
	if len(left) < MIN_HALF or len(right) < MIN_HALF:
		return None, "too_short"

	(full_seq, full_map), (gap_seq, gap_map) = _search_views(chain_residues)
	for seq, index_map in ((full_seq, full_map), (gap_seq, gap_map)):
		left_hits = _find_all(seq, left)
		right_hits = _find_all(seq, right)
		if len(left_hits) != 1 or len(right_hits) != 1:
			continue
		li = left_hits[0]
		ri = right_hits[0]
		# Halves must be contiguous on this search view (allow small
		# insertions between them, e.g. a single modified residue kept
		# in the full view).
		if ri < li + len(left) or ri > li + len(left) + 3:
			continue
		full_lo = index_map[li]
		full_hi = index_map[ri + len(right) - 1]
		if full_hi < full_lo:
			continue
		start_id = chain_residues[full_lo][0]
		end_id = chain_residues[full_hi][0]
		corrected = "".join(
			aa for _sid, aa, _gap in chain_residues[full_lo : full_hi + 1])
		# Reject runaway spans (halves matched far apart via insertions).
		if len(corrected) > len(peptide) + 5:
			continue
		return (start_id, end_id, corrected), None
	return None, "halves_not_unique"


def LocateMotif(chain_residues, peptide):
	"""Exact match, else half-inference.

	Returns ("ok"|"corrected", (start_id, end_id, seq)),
	("multi", hits), or ("missing", reason).
	"""
	result = FindUniqueMotif(chain_residues, peptide)
	if result is not None:
		status, payload = result
		if status == "multi":
			return "multi", payload
		return "ok", payload
	inferred, reason = InferByHalves(chain_residues, peptide)
	if inferred is not None:
		return "corrected", inferred
	return "missing", reason


def ReadMotifTsvFile(path, seq_cols):
	rows = []
	header_cols = None
	for line in open(path, encoding="utf-8", errors="replace"):
		line = line.rstrip("\n\r")
		if line == "":
			continue
		fields = line.split("\t")
		if fields[0] in RDRP_SKIP_COLS or fields[0] == "label":
			if header_cols is None:
				if seq_cols is not None:
					header_cols = list(seq_cols)
				elif fields[0] in RDRP_SKIP_COLS:
					header_cols = [
						f for f in fields[1:] if f not in RDRP_SKIP_COLS]
				else:
					header_cols = [
						c for c in fields[1:] if c not in CDN_SKIP_COLS]
			continue
		if header_cols is None:
			Die("no header row in %s" % path)
		if len(fields) < len(header_cols) + 1:
			Die("bad motif row in %s: %s" % (path, line[:80]))
		label = fields[0]
		motifs = {}
		for j, col in enumerate(header_cols):
			motifs[col] = MotifOrNone(fields[j + 1])
		rows.append((label, motifs))
	if len(rows) == 0:
		Die("no motif rows in %s" % path)
	return rows, header_cols


def ProcessRow(mode, label, motifs, seq_cols, cifs_dir, chainchar, chain_cache):
	if mode == "rdrp":
		pdb_id, chain = ParseRdrpLabel(label)
	else:
		pdb_id = ParseCdnLabel(label)
		chain = CDN_CHAIN

	cache_key = (str(cifs_dir), pdb_id)
	if cache_key not in chain_cache:
		cif_path = FindCifFile(cifs_dir, pdb_id)
		if cif_path is None:
			sys.stderr.write(
				"warning: CIF not found for %s (%s)\n" % (label, pdb_id))
			chain_cache[cache_key] = None
		else:
			chains, label_to_auth = ParseChainSequences(
				cif_path, chain_filter=chain)
			chain_cache[cache_key] = (cif_path, chains, label_to_auth)
	cached = chain_cache[cache_key]
	if cached is None:
		return [], []
	cif_path, all_chains, label_to_auth = cached
	auth_chain = ResolveChain(chain, all_chains, label_to_auth)
	if auth_chain is None:
		sys.stderr.write(
			"warning: chain %s not in %s (%s)\n"
			% (chain, cif_path, label))
		return [], []
	chain_residues = all_chains[auth_chain]
	query = MakeQuery(cif_path.stem, auth_chain, chainchar)
	out = []
	corrections = []
	for col in seq_cols:
		peptide = motifs.get(col)
		if peptide is None:
			continue
		status, payload = LocateMotif(chain_residues, peptide)
		if status == "missing":
			sys.stderr.write(
				"warning: motif %s not found in %s (%s): %s\n"
				% (col, query, label, peptide))
			continue
		if status == "multi":
			sys.stderr.write(
				"warning: motif %s found %d times in %s (%s); skipping\n"
				% (col, len(payload), query, label))
			continue
		start_id, end_id, matched = payload
		out.append((query, DB_NAME, col, start_id, end_id, matched))
		if status == "corrected" and matched != peptide:
			corrections.append(
				(query, col, peptide, matched, start_id, end_id))
			sys.stderr.write(
				"warning: motif %s corrected in %s (%s): %s -> %s\n"
				% (col, query, label, peptide, matched))
	return out, corrections


def main():
	default_motifs, default_cifs, default_cols = DefaultsFor(Args.mode)
	motifs_path = Args.motifs or default_motifs
	cifs_dir = Args.cifs or default_cifs
	if not os.path.isfile(motifs_path):
		Die("motifs file not found: %s" % motifs_path)
	if not os.path.isdir(cifs_dir):
		Die("cifs directory not found: %s" % cifs_dir)

	motif_rows, seq_cols = ReadMotifTsvFile(motifs_path, default_cols)
	chain_cache = {}
	all_rows = []
	all_corrections = []
	for label, motifs in motif_rows:
		rows, corrections = ProcessRow(
			Args.mode, label, motifs, seq_cols,
			cifs_dir, Args.chainchar, chain_cache)
		all_rows.extend(rows)
		all_corrections.extend(corrections)

	if Args.output == "-":
		fOut = sys.stdout
	else:
		fOut = open(Args.output, "w")
	fOut.write("query\tdb\tdbinfo\tstart\tend\tseq\n")
	for query, db, dbinfo, start, end, seq in all_rows:
		fOut.write("%s\t%s\t%s\t%d\t%d\t%s\n" % (query, db, dbinfo, start, end, seq))
	if Args.output != "-":
		fOut.close()

	if Args.report is not None:
		with open(Args.report, "w") as fRep:
			fRep.write("query\tdbinfo\toriginal\tcorrected\tstart\tend\n")
			for query, col, original, corrected, start, end in all_corrections:
				fRep.write(
					"%s\t%s\t%s\t%s\t%d\t%d\n"
					% (query, col, original, corrected, start, end))


if __name__ == "__main__":
	main()
