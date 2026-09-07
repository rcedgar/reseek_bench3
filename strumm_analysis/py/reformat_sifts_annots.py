#!/usr/bin/python3
"""Resolve SIFTS domain annotations to superfamily/family + description for figures.

Reads annots/{family}_sifts_annots.tsv and writes annots/{family}_sifts_fig_annots.tsv
with dbinfo rewritten:

  CATH  → "<SF> <description>"          e.g. 1.10.490.60 Phage p2 RNA dependent...
  SCOP  → "<SF_node> <description>"     e.g. 3000129 Nucleotidyltransferase-like
  Pfam  → "<family> <clan> <description>" e.g. PF00680 CL0027 Viral RNA-dependent...

SCOP2/SCOP2B are normalized to SCOP. Overlapping or adjacent segments that
share the same CATH homologous superfamily (4-level), SCOP SF node, or Pfam
family id are merged (seq rebuilt from FASTA). Other db types are dropped.
"""

import argparse
import os
import sys

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

SCOP_DBS = frozenset(("SCOP2", "SCOP2B", "SCOP"))


def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)


def ReadCathNames(fn):
	names = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		parts = Line.split(None, 1)
		sf = parts[0]
		names[sf] = parts[1].strip() if len(parts) > 1 else ""
	return names


def ReadScopClaToSf(fn):
	"""Map FA-DOMID and SF-DOMID → SF node id (e.g. 8041162 → 3000129)."""
	id_to_sf = {}
	for Line in open(fn):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		Fields = Line.split()
		if len(Fields) < 6:
			continue
		sf_node = None
		for part in Fields[-1].split(","):
			if part.startswith("SF="):
				sf_node = part[3:]
				break
		if sf_node is None:
			continue
		id_to_sf[Fields[0]] = sf_node
		id_to_sf[Fields[5]] = sf_node
	return id_to_sf


def ReadScopDes(fn):
	names = {}
	for Line in open(fn):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		parts = Line.split(None, 1)
		if len(parts) < 2:
			continue
		names[parts[0]] = parts[1].strip()
	return names


def ReadPfamClans(fn):
	"""acc → (clan_acc_or_dot, description)."""
	out = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0 or Line.startswith("#"):
			continue
		Fields = Line.split("\t")
		if len(Fields) < 4:
			continue
		acc = Fields[0].split(".", 1)[0]
		clan = Fields[1].strip() if len(Fields) > 1 else ""
		if clan == "" or clan.upper() in ("NO_CLAN", "NONE", "NA"):
			clan = "."
		desc = Fields[4].strip() if len(Fields) > 4 else ""
		if desc == "" and len(Fields) > 3:
			desc = Fields[3].strip()
		if acc not in out:
			out[acc] = (clan, desc)
	return out


def CathSf(dbinfo):
	"""Homologous superfamily id: C.A.T.H (four levels)."""
	raw = (dbinfo or "").strip().split("/", 1)[0]
	parts = raw.split(".")
	if len(parts) >= 4:
		return ".".join(parts[:4])
	return raw


def ScopDomId(dbinfo):
	return (dbinfo or "").strip().split("/", 1)[0]


def PfamAcc(dbinfo):
	return (dbinfo or "").strip().split(".", 1)[0]


def MergeAdjacentRows(rows, seqs_by_query):
	"""Merge overlapping or adjacent intervals with the same id (first dbinfo token).

	Adjacent means the next start is at most cur_end+1 (touching or overlapping).
	rows: list of (query, db, dbinfo, start, end, seq)
	Merged seq is taken from seqs_by_query[query][lo-1:hi] when available.
	"""
	by_key = {}
	for query, db, dbinfo, start, end, seq in rows:
		dom_id = dbinfo.split(None, 1)[0]
		by_key.setdefault((query, db, dom_id, dbinfo), []).append(
			(start, end, seq))
	out = []
	for (query, db, dom_id, dbinfo), ivals in by_key.items():
		ivals = sorted(ivals)
		cur_lo, cur_hi, _cur_seq = ivals[0]
		for lo, hi, _seq in ivals[1:]:
			if lo <= cur_hi + 1:
				cur_hi = max(cur_hi, hi)
			else:
				full = seqs_by_query.get(query)
				if full is not None and cur_lo >= 1 and cur_hi <= len(full):
					cur_seq = full[cur_lo - 1:cur_hi]
				else:
					cur_seq = _cur_seq
				out.append((query, db, dbinfo, cur_lo, cur_hi, cur_seq))
				cur_lo, cur_hi, _cur_seq = lo, hi, _seq
		full = seqs_by_query.get(query)
		if full is not None and cur_lo >= 1 and cur_hi <= len(full):
			cur_seq = full[cur_lo - 1:cur_hi]
		else:
			cur_seq = _cur_seq
		out.append((query, db, dbinfo, cur_lo, cur_hi, cur_seq))
	out.sort(key=lambda r: (r[0], r[1], r[3], r[4], r[2]))
	return out


def ResolveRow(db, dbinfo, cath_names, scop_id_to_sf, scop_des, pfam_clans):
	"""Return (out_db, out_dbinfo) or None to drop."""
	db = (db or "").strip()
	info = (dbinfo or "").strip()
	if db == "CATH":
		sf = CathSf(info)
		if sf not in cath_names:
			desc = "(missing)"
		elif cath_names[sf] == "":
			desc = "(unnamed)"
		else:
			desc = cath_names[sf]
		return "CATH", "%s %s" % (sf, desc)
	if db in SCOP_DBS:
		dom = ScopDomId(info)
		sf_node = scop_id_to_sf.get(dom)
		if sf_node is None:
			return "SCOP", "%s (missing)" % dom
		if sf_node not in scop_des:
			desc = "(missing)"
		elif scop_des[sf_node] == "":
			desc = "(unnamed)"
		else:
			desc = scop_des[sf_node]
		return "SCOP", "%s %s" % (sf_node, desc)
	if db == "Pfam":
		acc = PfamAcc(info)
		if acc not in pfam_clans:
			return "Pfam", "%s . (missing)" % acc
		clan, desc = pfam_clans[acc]
		if desc == "":
			desc = "(unnamed)"
		return "Pfam", "%s %s %s" % (acc, clan, desc)
	return None


def ReadFastaLens(fn):
	"""query_label → uppercase sequence (for merge seq rebuild)."""
	seqs = {}
	label = None
	parts = []

	def Flush():
		nonlocal label, parts
		if label is None:
			return
		seqs[label] = "".join(parts).replace(" ", "").upper()
		# also case-insensitive key
		seqs.setdefault(label.upper(), seqs[label])
		parts = []

	if not os.path.isfile(fn):
		return seqs
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		if Line.startswith(">"):
			Flush()
			label = Line[1:].split()[0]
			parts = []
		else:
			parts.append(Line.strip())
	Flush()
	return seqs


def Main():
	AP = argparse.ArgumentParser(description=__doc__)
	AP.add_argument("--family", choices=("rdrp", "cdn"), required=True)
	AP.add_argument("--input", default=None)
	AP.add_argument("--output", default=None)
	Args = AP.parse_args()

	ann = os.path.join(PalmDir, "annots")
	dat = os.path.join(PalmDir, "data")
	in_path = Args.input or os.path.join(ann, "%s_sifts_annots.tsv" % Args.family)
	out_path = Args.output or os.path.join(ann, "%s_sifts_fig_annots.tsv" % Args.family)
	fasta_path = os.path.join(dat, "%s_cif.fa" % Args.family)
	cath_fn = os.path.join(dat, "cath-b-newest-names")
	scop_cla = os.path.join(dat, "scop-cla-latest.txt")
	scop_des_fn = os.path.join(dat, "scop-des-latest.txt")
	pfam_fn = os.path.join(dat, "Pfam-A.clans.tsv")
	for fn in (in_path, cath_fn, scop_cla, scop_des_fn, pfam_fn):
		if not os.path.isfile(fn):
			Die("missing %s" % fn)

	cath_names = ReadCathNames(cath_fn)
	scop_id_to_sf = ReadScopClaToSf(scop_cla)
	scop_des = ReadScopDes(scop_des_fn)
	pfam_clans = ReadPfamClans(pfam_fn)
	seqs = ReadFastaLens(fasta_path)
	# Index by case-insensitive exact label
	seqs_by_query = {}
	for lab, seq in seqs.items():
		seqs_by_query[lab] = seq
		seqs_by_query[lab.upper()] = seq

	mergeable_rows = []
	other_rows = []
	n_in = 0
	n_drop = 0
	for Line in open(in_path):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] == "query":
			continue
		if len(Fields) < 6:
			Die("bad row: %s" % Line[:80])
		n_in += 1
		query, db, dbinfo = Fields[0], Fields[1], Fields[2]
		try:
			start = int(Fields[3])
			end = int(Fields[4])
		except ValueError:
			Die("bad coords: %s" % Line[:80])
		seq = Fields[5]
		resolved = ResolveRow(db, dbinfo, cath_names, scop_id_to_sf, scop_des,
			pfam_clans)
		if resolved is None:
			n_drop += 1
			continue
		out_db, out_info = resolved
		row = (query, out_db, out_info, start, end, seq)
		# Merge adjacent/overlapping CATH / SCOP / Pfam segments.
		if out_db in ("CATH", "SCOP", "Pfam"):
			mergeable_rows.append(row)
		else:
			other_rows.append(row)

	merge_seqs = {}
	for q, _db, _i, _s, _e, _seq in mergeable_rows:
		if q in merge_seqs:
			continue
		if q in seqs_by_query:
			merge_seqs[q] = seqs_by_query[q]
		elif q.upper() in seqs_by_query:
			merge_seqs[q] = seqs_by_query[q.upper()]

	mergeable_rows = MergeAdjacentRows(mergeable_rows, merge_seqs)
	out_rows = mergeable_rows + other_rows
	out_rows.sort(key=lambda r: (r[0], r[1], r[3], r[4], r[2]))

	with open(out_path, "w", encoding="utf-8") as fh:
		fh.write("query\tdb\tdbinfo\tstart\tend\tseq\n")
		for query, db, dbinfo, start, end, seq in out_rows:
			fh.write("%s\t%s\t%s\t%d\t%d\t%s\n" %
				(query, db, dbinfo, start, end, seq))

	sys.stderr.write("wrote %s (%d in, %d out, %d dropped)\n" %
		(out_path, n_in, len(out_rows), n_drop))


if __name__ == "__main__":
	Main()
