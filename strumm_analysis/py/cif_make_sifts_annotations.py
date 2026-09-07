#!/usr/bin/python3

import argparse
import sys
from pathlib import Path

from cif_seq import (
	Die,
	IterCifLoops,
	MakeQuery,
	ParseChainSequences,
	Present,
	ReadPathList,
	ResolveChain,
	SequenceSlice,
)

Usage = \
(
"Extract SIFTS domain annotations from mmCIF files. Input is a text "
"file with one CIF pathname per line, optionally followed by TAB and "
"chain id. Output is TSV with query, db, dbinfo, start, end, seq. "
"Segment sequences use the shared cif_seq polymer construction."
)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--input", required=True,
	help="Text file: one CIF path per line; optional TAB + chain id")
AP.add_argument("--output", default="-", help="Output TSV (default stdout)")
AP.add_argument("--chainchar", default="_",
	help="Separator between stem and chain in query id (default '_')")
AP.add_argument("--report", default=None, metavar="TEXTFILE",
	help="Write per-file annotation counts to TEXTFILE")
Args = AP.parse_args()

SEGMENTS_CAT = "pdbx_sifts_xref_db_segments"


def ParseSiftsSegments(cif_path):
	segments = []
	for cat, rows in IterCifLoops(cif_path):
		if cat == SEGMENTS_CAT:
			segments.extend(rows)
	return segments


def FormatDbinfo(xref_db, acc, domain_name):
	acc = (acc or "").strip()
	domain = (domain_name or "").strip()
	db = (xref_db or "").strip().upper()
	parts = []
	if Present(acc):
		parts.append(acc)
	if db == "PFAM":
		return parts[0] if parts else ""
	if Present(domain):
		if parts:
			return parts[0] + "/" + domain
		return domain
	return parts[0] if parts else ""


def ProcessFile(path, chain_filter, chainchar):
	cif_path = Path(path)
	if not cif_path.is_file():
		Die("not a file: %s" % path)
	segments = ParseSiftsSegments(cif_path)
	chains, label_to_auth = ParseChainSequences(
		cif_path, chain_filter=chain_filter)
	stem = cif_path.stem
	rows = []
	db_counts = {}
	if len(segments) == 0:
		sys.stderr.write("warning: no SIFTS domain segments: %s\n" % path)
		return rows, db_counts

	auth_filter = None
	if chain_filter is not None:
		auth_filter = ResolveChain(chain_filter, chains, label_to_auth)
		if auth_filter is None:
			sys.stderr.write(
				"warning: chain %s not in %s\n" % (chain_filter, path))
			return rows, db_counts

	for seg in segments:
		asym = (seg.get("asym_id") or "").strip()
		auth_chain = ResolveChain(asym, chains, label_to_auth)
		if auth_chain is None:
			continue
		if auth_filter is not None and auth_chain != auth_filter:
			continue
		xref_db = (seg.get("xref_db") or "").strip()
		if not Present(xref_db):
			continue
		try:
			start = int((seg.get("seq_id_start") or "").strip())
			end = int((seg.get("seq_id_end") or "").strip())
		except ValueError:
			sys.stderr.write(
				"warning: bad coords in %s chain %s db %s\n"
				% (path, auth_chain, xref_db))
			continue
		if end < start:
			sys.stderr.write(
				"warning: inverted coords in %s chain %s db %s: %d-%d\n"
				% (path, auth_chain, xref_db, start, end))
			continue
		dbinfo = FormatDbinfo(
			xref_db,
			seg.get("xref_db_acc"),
			seg.get("domain_name"))
		seq = SequenceSlice(chains[auth_chain], start, end)
		if "X" in seq:
			sys.stderr.write(
				"warning: missing residue letters in %s chain %s "
				"db %s coords %d-%d\n"
				% (path, auth_chain, xref_db, start, end))
		query = MakeQuery(stem, auth_chain, chainchar)
		rows.append((query, auth_chain, xref_db, dbinfo, start, end, seq))
		db_key = xref_db.lower()
		db_counts[db_key] = db_counts.get(db_key, 0) + 1
	rows.sort(key=lambda r: (r[1], r[2], r[4], r[5]))
	return rows, db_counts


def FormatReport(db_counts):
	if len(db_counts) == 0:
		return "(none)"
	parts = []
	for db in sorted(db_counts.keys()):
		parts.append("F%s=%d" % (db, db_counts[db]))
	return " ".join(parts)


def main():
	entries = ReadPathList(Args.input)
	if len(entries) == 0:
		Die("no input paths in %s" % Args.input)
	all_rows = []
	fReport = open(Args.report, "w") if Args.report is not None else None
	try:
		for path, chain_filter, _lineno in entries:
			rows, db_counts = ProcessFile(path, chain_filter, Args.chainchar)
			all_rows.extend(rows)
			if fReport is not None:
				fReport.write("%s\t%s\n" % (path, FormatReport(db_counts)))
	finally:
		if fReport is not None:
			fReport.close()
	if Args.output == "-":
		fOut = sys.stdout
	else:
		fOut = open(Args.output, "w")
	fOut.write("query\tdb\tdbinfo\tstart\tend\tseq\n")
	for query, _asym, xref_db, dbinfo, start, end, seq in all_rows:
		fOut.write(
			"%s\t%s\t%s\t%d\t%d\t%s\n"
			% (query, xref_db, dbinfo, start, end, seq))
	if Args.output != "-":
		fOut.close()


if __name__ == "__main__":
	main()
