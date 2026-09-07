#!/usr/bin/python3
"""Draw SIFTS domain lanes + rce motif trapezoids for rdrp/cdn figure tips."""

import argparse
import os
import re
import sys

try:
	import matplotlib
	matplotlib.use("Agg")
	matplotlib.rcParams["pdf.fonttype"] = 42
	matplotlib.rcParams["ps.fonttype"] = 42
	import matplotlib.pyplot as plt
	from matplotlib.colors import to_rgba
	from matplotlib.patches import Polygon, Rectangle
except ImportError:
	sys.stderr.write("matplotlib is required\n")
	sys.exit(1)

Usage = (
"Draw amino-acid-scale protein bars with SIFTS domain fills "
"(CATH / SCOP / Pfam lanes) and rce motif trapezoid connectors. "
"Tip order comes from data/{family}_figure_order.txt."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

# ---------------------------------------------------------------------------
# Heatmap color scales per SIFTS db type (edit here)
# ---------------------------------------------------------------------------

# MATLAB parula (sampled); matplotlib has no built-in "parula".
_PARULA_RGB = [
	(0.2081, 0.1663, 0.5292),
	(0.2116, 0.1898, 0.5777),
	(0.2123, 0.2138, 0.6268),
	(0.2081, 0.2386, 0.6771),
	(0.1959, 0.2645, 0.7279),
	(0.1707, 0.2919, 0.7792),
	(0.1253, 0.3242, 0.8303),
	(0.0591, 0.3598, 0.8683),
	(0.0117, 0.3875, 0.8820),
	(0.0060, 0.4086, 0.8828),
	(0.0165, 0.4266, 0.8786),
	(0.0329, 0.4430, 0.8720),
	(0.0498, 0.4586, 0.8641),
	(0.0629, 0.4737, 0.8554),
	(0.0723, 0.4887, 0.8467),
	(0.0779, 0.5040, 0.8384),
	(0.0793, 0.5200, 0.8312),
	(0.0749, 0.5375, 0.8263),
	(0.0641, 0.5570, 0.8240),
	(0.0488, 0.5772, 0.8228),
	(0.0343, 0.5966, 0.8199),
	(0.0265, 0.6137, 0.8135),
	(0.0239, 0.6287, 0.8038),
	(0.0358, 0.6418, 0.7912),
	(0.0631, 0.6540, 0.7765),
	(0.1048, 0.6658, 0.7606),
	(0.1606, 0.6773, 0.7441),
	(0.2396, 0.6886, 0.7274),
	(0.3351, 0.6997, 0.7104),
	(0.4362, 0.7099, 0.6931),
	(0.5301, 0.7189, 0.6760),
	(0.6104, 0.7268, 0.6598),
	(0.6777, 0.7337, 0.6438),
	(0.7359, 0.7400, 0.6272),
	(0.7883, 0.7459, 0.6096),
	(0.8353, 0.7520, 0.5900),
	(0.8772, 0.7588, 0.5678),
	(0.9139, 0.7670, 0.5418),
	(0.9455, 0.7773, 0.5111),
	(0.9721, 0.7901, 0.4751),
	(0.9931, 0.8056, 0.4329),
	(0.9979, 0.8263, 0.3859),
	(0.9841, 0.8494, 0.3402),
	(0.9565, 0.8739, 0.3027),
	(0.9231, 0.8981, 0.2743),
	(0.8882, 0.9213, 0.2519),
	(0.8550, 0.9427, 0.2325),
	(0.8265, 0.9621, 0.2139),
	(0.8039, 0.9793, 0.1955),
	(0.7882, 0.9931, 0.1769),
	(0.7764, 0.9973, 0.1580),
	(0.7671, 0.9981, 0.1402),
	(0.7591, 0.9967, 0.1236),
	(0.7517, 0.9938, 0.1090),
	(0.7446, 0.9898, 0.0965),
	(0.7372, 0.9850, 0.0860),
	(0.7287, 0.9796, 0.0773),
	(0.7181, 0.9739, 0.0700),
	(0.7046, 0.9680, 0.0639),
	(0.6871, 0.9621, 0.0590),
	(0.6650, 0.9560, 0.0554),
	(0.6379, 0.9495, 0.0532),
	(0.6062, 0.9425, 0.0526),
	(0.5716, 0.9349, 0.0535),
]
_ParulaRegistered = False


def EnsureParula():
	global _ParulaRegistered
	if _ParulaRegistered:
		return
	try:
		plt.get_cmap("parula")
		_ParulaRegistered = True
		return
	except ValueError:
		pass
	from matplotlib.colors import LinearSegmentedColormap
	cmap = LinearSegmentedColormap.from_list("parula", _PARULA_RGB, N=256)
	try:
		plt.colormaps.register(cmap)
	except AttributeError:
		plt.register_cmap(cmap=cmap)
	_ParulaRegistered = True


def CmapHex(cmap_name, i, n, t0=0.30, t1=0.85):
	if cmap_name == "parula":
		EnsureParula()
	cmap = plt.get_cmap(cmap_name)
	if n <= 1:
		t = 0.5 * (t0 + t1)
	else:
		t = t0 + (t1 - t0) * (float(i) / float(max(1, n - 1)))
	r, g, b, _a = cmap(t)
	return "#%02X%02X%02X" % (int(round(r * 255)), int(round(g * 255)),
		int(round(b * 255)))

def PfamColor(i, n):
	return CmapHex("parula", i, n, t0=0.15, t1=0.90)

def ScopColor(i, n):
	return CmapHex("RdPu", i, n, t0=0.35, t1=0.85)

def CathColor(i, n):
	return CmapHex("GnBu", i, n, t0=0.35, t1=0.85)

DbColorFn = {
	"Pfam": PfamColor,
	"SCOP": ScopColor,
	"CATH": CathColor,
}

# ---------------------------------------------------------------------------

LocusColorsAbc = {
	"A": "dodgerblue",
	"B": "palegreen",
	"C": "magenta",
}

LocusColorsIii = {
	"I": "dodgerblue",
	"II": "palegreen",
	"III": "magenta",
}

LabelColors = {
	"rdrp": "#1B4F72",
	"decoy": "#922B21",
	"other": "#1B4F72",
}

LaneUnit = 0.27  # one annotation-lane height; bar height = n units (min 1)
BarGap = LaneUnit * 0.70  # uniform vertical gap between bars (also label↔bar gap)
SfFs = 5.5
LabelFs = 8.0
LegFs = 6.0
FillAlpha = 0.10
LocusBarAlpha = 0.20
HeaderH = 0.32
AbcPad = 30
LeftExtra = 20
SfHiPad = 20
ResViewLo = 1
# Crop to the window currently labeled 400..1200 on the prior 1..1500
# renumbering of original ticks 500..2000 (= original ticks 899..1699),
# then renumber displayed ticks 1..(1200-400)=800.
CropTickLo = 500 + 400 - 1  # 899
CropTickHi = 500 + 1200 - 1  # 1699
TruncTriW = 10  # amino-acid width of truncation triangle
BarWidthScale = 0.70  # physical width of bar axes vs prior layout

RdrpLocusCols = ("A_seq", "B_seq", "C_seq")
RdrpLocusNames = ("A", "B", "C")
CdnLocusCols = ("seqI", "seqII", "seqIII")
CdnLocusNames = ("I", "II", "III")

SCOP_DBS = frozenset(("SCOP2", "SCOP2B"))


def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)


def NormKey(label):
	return label.strip().upper().replace("_", "-")


def PdbOf(label):
	s = label.strip()
	if len(s) >= 4:
		return s[:4].upper()
	return s.upper()


def StripChainKeepCategory(label):
	"""1hhs_A-rdrp -> 1hhs-rdrp; 4TY0-CdnA unchanged; 4KLQ_Polbe -> 4KLQ-Polbe."""
	s = label.strip()
	m = re.match(r"^([0-9A-Za-z]{4})_([A-Za-z0-9])-(.+)$", s)
	if m:
		return "%s-%s" % (m.group(1), m.group(3))
	m = re.match(r"^([0-9A-Za-z]{4})_([A-Za-z0-9])$", s)
	if m:
		return m.group(1)
	m = re.match(r"^([0-9A-Za-z]{4})_(.+)$", s)
	if m and len(m.group(2)) > 1:
		return "%s-%s" % (m.group(1), m.group(2))
	return s


def ClassOf(disp_label):
	s = disp_label.strip()
	if s.endswith("-rdrp") or s.upper().endswith("-RDRP"):
		return "rdrp"
	if s.endswith("-decoy") or s.upper().endswith("-DECOY"):
		return "decoy"
	return "other"


def ParseDbinfo(db, dbinfo):
	"""Return (bar_id, clan_or_None, description).

	CATH/SCOP dbinfo: \"<id> <description...>\"
	Pfam dbinfo: \"<family> <clan> <description...>\"
	"""
	parts = (dbinfo or "").strip().split()
	if len(parts) == 0:
		return "?", None, "(missing)"
	if db == "Pfam":
		fam = parts[0]
		clan = parts[1] if len(parts) > 1 else "."
		desc = " ".join(parts[2:]) if len(parts) > 2 else "(missing)"
		return fam, clan, desc
	dom_id = parts[0]
	desc = " ".join(parts[1:]) if len(parts) > 1 else "(missing)"
	return dom_id, None, desc


def LaneDbOf(tsv_db, wanted):
	"""Map TSV db field to lane name, or None if not selected."""
	db = (tsv_db or "").strip()
	if db == "Pfam" and "Pfam" in wanted:
		return "Pfam"
	if db == "CATH" and "CATH" in wanted:
		return "CATH"
	if db == "SCOP" and "SCOP" in wanted:
		return "SCOP"
	# Accept legacy SCOP2/SCOP2B if present
	if db in SCOP_DBS and "SCOP" in wanted:
		return "SCOP"
	return None


def FamilyPaths(family):
	ann = os.path.join(PalmDir, "annots")
	dat = os.path.join(PalmDir, "data")
	return {
		"sifts": os.path.join(ann, "%s_sifts_fig_annots.tsv" % family),
		"motifs": os.path.join(ann, "%s_motif_annots.tsv" % family),
		"order": os.path.join(dat, "%s_figure_order.txt" % family),
		"fasta": os.path.join(dat, "%s_cif.fa" % family),
		"category": os.path.join(dat, "%s_motifs.tsv" % family),
		"critical": os.path.join(ann, "%s_critical.tsv" % family),
	}


def ReadCriticalTsv(fn):
	"""Return {query: (GDD, gate)} from extract_critical.py output."""
	out = {}
	if not os.path.isfile(fn):
		return out
	with open(fn) as f:
		for line in f:
			line = line.rstrip("\n")
			if len(line) == 0:
				continue
			fields = line.split("\t")
			if fields[0] == "query":
				continue
			if len(fields) < 4:
				continue
			query, _dgd, gdd, gate = fields[0], fields[1], fields[2], fields[3]
			out[query] = (gdd, gate)
	return out


def ParseDbtypes(s):
	wanted = []
	seen = set()
	for tok in s.split(","):
		t = tok.strip()
		if t == "":
			continue
		# Normalize spelling
		low = t.lower()
		if low == "scop":
			name = "SCOP"
		elif low == "cath":
			name = "CATH"
		elif low == "pfam":
			name = "Pfam"
		else:
			Die("unknown dbtype %s (use CATH,SCOP,Pfam)" % t)
		if name not in seen:
			seen.add(name)
			wanted.append(name)
	if len(wanted) == 0:
		Die("empty --dbtypes")
	return wanted


def ReadFasta(fn):
	seqs = {}
	by_key = {}
	label = None
	parts = []

	def Flush():
		nonlocal label, parts
		if label is None:
			return
		seq = "".join(parts).replace(" ", "").upper()
		if len(seq) == 0:
			Die("empty sequence %s" % label)
		if label not in seqs:
			seqs[label] = seq
		key = NormKey(label)
		if key not in by_key:
			by_key[key] = label
		# Also index by PDB for cdn order ids
		pdb = PdbOf(label)
		if pdb not in by_key:
			by_key[pdb] = label
		parts = []

	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		if Line.startswith(">"):
			Flush()
			label = Line[1:].split()[0]
			parts = []
		else:
			if label is None:
				Die("sequence before header in %s" % fn)
			parts.append(Line.strip())
	Flush()
	if len(seqs) == 0:
		Die("no sequences in fasta %s" % fn)
	return seqs, by_key


def HasChainSuffix(label):
	"""True for PDB_chain like 6f5p_B (single-char chain), not 4KLQ_Polbe."""
	return re.match(r"^[0-9A-Za-z]{4}_[A-Za-z0-9]$", label.strip()) is not None


def ResolveFastaLabel(query, seqs, by_key, allow_pdb_fallback=False):
	"""Map query to a FASTA header.

	Chain-qualified queries (e.g. 6f5p_A) must match that chain; PDB-only
	order ids (e.g. 4TY0) may use allow_pdb_fallback=True.
	"""
	if query in seqs:
		return query
	key = NormKey(query)
	# Prefer exact normalized label (not PDB-only index).
	for lab in seqs:
		if NormKey(lab) == key:
			return lab
	for lab in seqs:
		if lab.upper() == query.upper():
			return lab
	if allow_pdb_fallback and not HasChainSuffix(query):
		pdb = PdbOf(query)
		if pdb in by_key:
			return by_key[pdb]
	return None


def CheckSubseq(fasta_lab, seqs, start, end, seq, tag):
	full = seqs[fasta_lab]
	if start < 1 or end < start or end > len(full):
		Die("%s coords %d-%d out of range for %s (len %d)" %
			(tag, start, end, fasta_lab, len(full)))
	got = full[start - 1 : end]
	want = seq.replace(" ", "").upper()
	if got != want:
		Die("%s seq mismatch %s %d-%d\n  fasta=%s\n  tsv =%s" %
			(tag, fasta_lab, start, end, got[:80], want[:80]))


def ReadOrder(fn):
	order = []
	seen = set()
	for Line in open(fn):
		tok = Line.strip()
		if len(tok) == 0 or tok.startswith("#"):
			continue
		key = tok.upper()
		if key in seen:
			Die("duplicate order id %s" % tok)
		seen.add(key)
		order.append(tok)
	if len(order) == 0:
		Die("empty order file %s" % fn)
	return order


def ReadCategoryMap(fn, family):
	"""Map NormKey(chain tip) / PDB -> category-suffixed motif label."""
	by_chain = {}
	by_pdb = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		lab = Fields[0].strip()
		if lab.lower() in ("label",):
			continue
		# rdrp: 1hhs_A-rdrp → base 1hhs_A
		base = lab
		if family == "rdrp":
			if base.endswith("-rdrp"):
				base = base[:-len("-rdrp")]
			elif base.endswith("-decoy"):
				base = base[:-len("-decoy")]
			by_chain[NormKey(base)] = lab
			by_pdb[PdbOf(base)] = lab
		else:
			# cdn labels are PDB_Category or PDB-Category (no chain)
			by_pdb[PdbOf(lab)] = lab
			by_chain[NormKey(lab)] = lab
	return by_chain, by_pdb


def TipCategoryLabel(fasta_lab, by_chain, by_pdb):
	k = NormKey(fasta_lab)
	if k in by_chain:
		return by_chain[k]
	# strip optional chain for lookup
	m = re.match(r"^([0-9A-Za-z]{4})_([A-Za-z0-9])$", fasta_lab)
	if m:
		k2 = NormKey(m.group(1) + "_" + m.group(2))
		if k2 in by_chain:
			return by_chain[k2]
	pdb = PdbOf(fasta_lab)
	if pdb in by_pdb:
		return by_pdb[pdb]
	return fasta_lab


def CathSortKey(sf):
	parts = sf.split(".")
	vals = []
	for p in parts[:4]:
		try:
			vals.append(int(p))
		except ValueError:
			vals.append(10 ** 9)
	while len(vals) < 4:
		vals.append(0)
	return tuple(vals)


def ScopSortKey(sf):
	try:
		return (int(sf),)
	except ValueError:
		return (10 ** 18,)


def ContainedTextMask(doms):
	show = []
	for i, (_sf, s, e) in enumerate(doms):
		has_inner = False
		for j, (_sf2, s2, e2) in enumerate(doms):
			if i == j:
				continue
			if s <= s2 and e2 <= e and (s2 > s or e2 < e):
				has_inner = True
				break
		show.append(not has_inner)
	return show


def MergeSameIdDoms(doms):
	"""Merge consecutive same-id spans into (id, min_start, max_end)."""
	if len(doms) == 0:
		return []
	ordered = sorted(doms, key=lambda d: (d[1], d[2], d[0]))
	out = []
	cur_sf, cur_lo, cur_hi = ordered[0]
	for sf, lo, hi in ordered[1:]:
		if sf == cur_sf:
			if lo < cur_lo:
				cur_lo = lo
			if hi > cur_hi:
				cur_hi = hi
			continue
		out.append((cur_sf, cur_lo, cur_hi))
		cur_sf, cur_lo, cur_hi = sf, lo, hi
	out.append((cur_sf, cur_lo, cur_hi))
	return out


def SfArchitecture(doms):
	accs = []
	seen = set()
	for sf, start, end in sorted(doms, key=lambda d: (d[1], d[2], d[0])):
		if sf in seen:
			continue
		seen.add(sf)
		accs.append(sf)
	return tuple(accs)


def MotifCenter(m_mid, qlen):
	if m_mid is None:
		return (1.0 + float(qlen)) / 2.0
	return (m_mid[0] + m_mid[1]) / 2.0


def ReadAnnotTsv(fn, wanted, seqs, by_key, locus_cols):
	"""Return (domains_by_tip, motifs_by_tip, id_meta).

	domains_by_tip[fasta_lab][lane_db] = [(id, lo, hi), ...]
	motifs_by_tip[fasta_lab] = {locus_col: (lo, hi)}
	id_meta[(lane_db, id)] = (clan_or_None, description)
	"""
	domains = {}
	motifs = {}
	id_meta = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] == "query":
			continue
		if len(Fields) < 6:
			Die("bad annot row: %s" % Line[:80])
		query, db, dbinfo = Fields[0], Fields[1], Fields[2]
		try:
			start = int(Fields[3])
			end = int(Fields[4])
		except ValueError:
			Die("bad coords: %s" % Line[:80])
		seq = Fields[5]
		fasta_lab = ResolveFastaLabel(query, seqs, by_key, allow_pdb_fallback=False)
		if fasta_lab is None:
			sys.stderr.write("skip no fasta for annot query %s\n" % query)
			continue
		CheckSubseq(fasta_lab, seqs, start, end, seq, "%s:%s" % (fn, db))

		if db == "rce_motif":
			if dbinfo in locus_cols:
				motifs.setdefault(fasta_lab, {})[dbinfo] = (start, end)
			continue

		lane = LaneDbOf(db, wanted)
		if lane is None:
			continue
		did, clan, desc = ParseDbinfo(lane, dbinfo)
		domains.setdefault(fasta_lab, {}).setdefault(lane, []).append(
			(did, start, end))
		if (lane, did) not in id_meta:
			id_meta[(lane, did)] = (clan, desc)

	for lab in domains:
		for lane in domains[lab]:
			domains[lab][lane] = MergeSameIdDoms(domains[lab][lane])
	return domains, motifs, id_meta


def Main():
	AP = argparse.ArgumentParser(description=Usage)
	AP.add_argument("--family", choices=("rdrp", "cdn"), required=True)
	AP.add_argument("--dbtypes", default="CATH,SCOP,Pfam",
		help="Comma-separated SIFTS db types (SCOP from fig annots)")
	AP.add_argument("--output", default=None,
		help="Comma-separated output files (default FAMILY.pdf,FAMILY.png)")
	Args = AP.parse_args()

	family = Args.family
	wanted = ParseDbtypes(Args.dbtypes)
	paths = FamilyPaths(family)
	for key in ("sifts", "motifs", "order", "fasta", "category"):
		if not os.path.isfile(paths[key]):
			Die("missing %s (run reformat_sifts_annots.py first for sifts)" %
				paths[key])

	if Args.output is None:
		out_paths = ["%s.pdf" % family, "%s.png" % family]
	else:
		out_paths = [p.strip() for p in Args.output.split(",") if p.strip()]
		if len(out_paths) == 0:
			Die("empty --output")

	if family == "rdrp":
		locus_cols = RdrpLocusCols
		locus_names = RdrpLocusNames
		locus_colors = LocusColorsAbc
	else:
		locus_cols = CdnLocusCols
		locus_names = CdnLocusNames
		locus_colors = LocusColorsIii

	seqs, by_key = ReadFasta(paths["fasta"])
	order_ids = ReadOrder(paths["order"])
	cat_by_chain, cat_by_pdb = ReadCategoryMap(paths["category"], family)
	critical = ReadCriticalTsv(paths["critical"])

	# Load + validate both annot tables against FASTA.
	sifts_doms, _, id_meta = ReadAnnotTsv(
		paths["sifts"], wanted, seqs, by_key, locus_cols)
	_, motif_spans, _ = ReadAnnotTsv(
		paths["motifs"], wanted, seqs, by_key, locus_cols)
	domains = sifts_doms

	tip_order = []
	tip_motifs = {}
	tip_domains = {}
	disp_labels = {}
	for oid in order_ids:
		fasta_lab = ResolveFastaLabel(oid, seqs, by_key, allow_pdb_fallback=True)
		if fasta_lab is None:
			sys.stderr.write("skip no fasta for order id %s\n" % oid)
			continue
		qlen = len(seqs[fasta_lab])
		spans = []
		ms = motif_spans.get(fasta_lab, {})
		for col in locus_cols:
			spans.append(ms.get(col))
		tip_order.append(fasta_lab)
		tip_motifs[fasta_lab] = (spans[0], spans[1], spans[2], qlen)
		tip_domains[fasta_lab] = domains.get(fasta_lab, {})
		cat = TipCategoryLabel(fasta_lab, cat_by_chain, cat_by_pdb)
		disp_labels[fasta_lab] = StripChainKeepCategory(cat)

	if len(tip_order) == 0:
		Die("no tips to draw")

	def TipAnnotKey(label):
		"""Sort: rdrp before decoy, then similar domain architectures together."""
		cls = ClassOf(disp_labels[label])
		cls_rank = {"rdrp": 0, "decoy": 1, "other": 2}.get(cls, 3)
		doms = tip_domains.get(label, {})
		arch = []
		for db in wanted:
			arch.append(SfArchitecture(doms.get(db, [])))
		gdd, gate = critical.get(label, (".", "."))
		return (cls_rank, arch, gdd or ".", gate or ".", disp_labels[label])

	tip_order.sort(key=TipAnnotKey)

	def TipLanes(label):
		lanes = []
		doms = tip_domains.get(label, {})
		for db in wanted:
			if len(doms.get(db, [])) > 0:
				lanes.append(db)
		return lanes

	def TipBarH(label):
		# Each DB is one LaneUnit tall; ≤1 DB still uses one unit.
		n = len(TipLanes(label))
		return float(max(1, n)) * LaneUnit

	def TipLaneH(label):
		return LaneUnit

	tip_y = {}
	y_cursor = 0.0
	for i, label in enumerate(tip_order):
		bh = TipBarH(label)
		if i > 0:
			y_cursor -= BarGap
		tip_top = y_cursor
		tip_y[label] = tip_top - bh / 2.0
		y_cursor = tip_top - bh
	y_bottom = y_cursor

	sf_sets = {db: set() for db in wanted}
	for label in tip_order:
		doms = tip_domains.get(label, {})
		for db in wanted:
			for sf in SfArchitecture(doms.get(db, [])):
				sf_sets[db].add(sf)

	sf_order = {}
	for db in wanted:
		sfs = list(sf_sets[db])
		if db == "CATH":
			sfs.sort(key=CathSortKey)
		elif db == "SCOP":
			sfs.sort(key=ScopSortKey)
		else:
			sfs.sort()
		sf_order[db] = sfs

	SfColors = {}
	for db in wanted:
		sfs = sf_order[db]
		fn = DbColorFn[db]
		for i, sf in enumerate(sfs):
			SfColors[(db, sf)] = fn(i, max(1, len(sfs)))

	def LegText(db, sf):
		clan, desc = id_meta.get((db, sf), (None, "(missing)"))
		if db == "Pfam":
			return "%s %s %s" % (sf, clan if clan else ".", desc)
		return "%s %s" % (sf, desc)

	def AnchorOf(label):
		mA, mB, mC, qlen = tip_motifs[label]
		return MotifCenter(mB, qlen)

	x0 = None
	for label in tip_order:
		center = AnchorOf(label)
		left = 1.0 - center
		if x0 is None or left < x0:
			x0 = left
	shift = -x0 + 1.0

	def XOf(res, center):
		return float(res) - center + shift

	def SegX(lo, hi, center):
		return XOf(lo, center) - 0.5, XOf(hi, center) + 0.5

	def ClipSeg(x1, x2, lo, hi):
		a = max(x1, lo)
		b = min(x2, hi)
		if b <= a:
			return None
		return a, b

	# Full-length bars first, then crop to former tick window 500..2000.
	XViewLo = None
	XViewHi = None
	for label in tip_order:
		mA, mB, mC, qlen = tip_motifs[label]
		center = AnchorOf(label)
		x1, x2 = SegX(1, qlen, center)
		if XViewLo is None or x1 < XViewLo:
			XViewLo = x1
		if XViewHi is None or x2 > XViewHi:
			XViewHi = x2
		if mA is not None:
			lo = max(1, mA[0] - AbcPad)
			ax1 = XOf(lo, center) - 0.5
			if ax1 < XViewLo:
				XViewLo = ax1
	XViewLo -= LeftExtra
	XViewHi += SfHiPad
	if XViewHi is None or XViewHi <= XViewLo:
		Die("empty bar x range")
	# Tick label s was at XViewLo + (s - 1); crop to labels CropTickLo..CropTickHi.
	XViewLo = XViewLo + float(CropTickLo - 1)
	XViewHi = XViewLo + float(CropTickHi - CropTickLo)
	if XViewHi <= XViewLo:
		Die("empty cropped bar x range")

	fig_w = 12.0
	y_top_extent = HeaderH + 1.05
	y_bot_extent = y_bottom - 0.20
	right_in = 0.18
	n_leg_items = sum(len(sf_order[db]) for db in wanted)
	n_leg_hdrs = sum(1 for db in wanted if len(sf_order[db]) > 0)
	leg_row_in = 0.155
	leg_hdr_in = 0.20
	leg_est_in = (0.35 + n_leg_hdrs * leg_hdr_in + n_leg_items * leg_row_in
		+ 0.10 * max(0, n_leg_hdrs - 1))
	bottom_in = max(2.2, leg_est_in + 0.70)
	top_in = 0.50
	data_span = y_top_extent - y_bot_extent
	fig_h = max(8.5, data_span * 0.42 + bottom_in + top_in + 1.0)
	ax_h_frac = 1.0 - (bottom_in + top_in) / fig_h
	ax_y_frac = bottom_in / fig_h

	fig = plt.figure(figsize=(fig_w, fig_h))
	ax_bar = fig.add_axes([0.20, ax_y_frac, 0.74, ax_h_frac])
	ax_bar.set_ylim(y_bot_extent, y_top_extent)

	def ProbeIn(s, fs=LabelFs, family_font="monospace"):
		t = fig.text(0, 0, s, fontsize=fs, fontfamily=family_font)
		fig.canvas.draw()
		w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
		t.remove()
		return w

	bar_ax_in = fig_w * 0.74
	bar_left_in = fig_w * 0.20
	plot_w = bar_left_in + bar_ax_in + right_in

	leg_groups = [(db, sf_order[db]) for db in wanted]
	longest_leg = "Pfam"
	for db, sfs in leg_groups:
		for sf in sfs:
			txt = LegText(db, sf)
			if len(txt) > len(longest_leg):
				longest_leg = txt
	leg_txt_in = ProbeIn(longest_leg, fs=LegFs, family_font="sans-serif") if longest_leg else 0.0
	needed_leg_in = 0.28 + 0.14 + leg_txt_in + 0.50
	fig_w = max(plot_w, needed_leg_in + 0.40)
	xoff = 0.5 * (fig_w - plot_w)
	fig.set_size_inches(fig_w, fig_h)
	bar_ax_in = fig_w - (bar_left_in + right_in)
	# Shrink bar panel to BarWidthScale; give freed width to left margin.
	freed_in = bar_ax_in * (1.0 - BarWidthScale)
	bar_ax_in *= BarWidthScale
	bar_left_in += freed_in
	ax_bar.set_position([(xoff + bar_left_in) / fig_w, ax_y_frac,
		bar_ax_in / fig_w, ax_h_frac])

	# Match a small pad between tip-text block and leftmost possible bar.
	ax_h_in = fig_h * ax_h_frac
	gap_in = BarGap / data_span * ax_h_in
	x_span = XViewHi - XViewLo
	label_gap_x = gap_in * x_span / bar_ax_in

	# Tip text: right-aligned label, then GDD gate (GDD column vertically aligned).
	max_lab_len = 1
	tip_gdd = {}
	tip_gate = {}
	for label in tip_order:
		gdd, gate = critical.get(label, (None, None))
		if gdd is None:
			sys.stderr.write("warning: no critical row for %s\n" % label)
			gdd, gate = ".", "."
		tip_gdd[label] = gdd
		tip_gate[label] = gate
		max_lab_len = max(max_lab_len, len(disp_labels[label]))
	lab_field = "X" * max_lab_len
	gdd_gate_sample = "WWW X"
	full_sample = "%s %s" % (lab_field, gdd_gate_sample)
	tip_txt_in = ProbeIn(full_sample, fs=LabelFs, family_font="monospace")
	space_in = ProbeIn(" ", fs=LabelFs, family_font="monospace")
	need_left_in = tip_txt_in + 0.30
	if bar_left_in < need_left_in:
		extra = need_left_in - bar_left_in
		fig_w += extra
		bar_left_in += extra
		plot_w = bar_left_in + bar_ax_in + right_in
		xoff = 0.5 * (fig_w - plot_w)
		fig.set_size_inches(fig_w, fig_h)
		ax_bar.set_position([(xoff + bar_left_in) / fig_w, ax_y_frac,
			bar_ax_in / fig_w, ax_h_frac])
		label_gap_x = gap_in * x_span / bar_ax_in
	space_x = space_in * x_span / bar_ax_in
	gdd_gate_w_x = ProbeIn(gdd_gate_sample, fs=LabelFs, family_font="monospace") * x_span / bar_ax_in
	# Right edge of tip block sits label_gap_x left of XViewLo.
	block_right = XViewLo - label_gap_x
	gdd_x = block_right - gdd_gate_w_x
	lab_right_x = gdd_x - space_x

	for label in tip_order:
		cls = ClassOf(disp_labels[label])
		col = LabelColors[cls]
		ax_bar.text(lab_right_x, tip_y[label], disp_labels[label],
			va="center", ha="right", fontsize=LabelFs, fontfamily="monospace",
			color=col, clip_on=False, zorder=2)
		ax_bar.text(gdd_x, tip_y[label],
			"%s %s" % (tip_gdd[label], tip_gate[label]),
			va="center", ha="left", fontsize=LabelFs, fontfamily="monospace",
			color=col, clip_on=False, zorder=2)

	used_sf = {}
	by_y = list(tip_order)

	# Header badges
	for mi, name in enumerate(locus_names):
		top_lab = None
		for lab in by_y:
			if tip_motifs[lab][mi] is not None:
				top_lab = lab
				break
		if top_lab is None:
			continue
		mot = tip_motifs[top_lab][mi]
		cen_top = AnchorOf(top_lab)
		y_top_bar = tip_y[top_lab] + TipBarH(top_lab) / 2.0
		clip = ClipSeg(*SegX(mot[0], mot[1], cen_top), XViewLo, XViewHi)
		if clip is None:
			continue
		vx1, vx2 = clip
		col = locus_colors[name]
		ax_bar.add_patch(Polygon(
			((vx1, y_top_bar), (vx2, y_top_bar),
				(vx2, y_top_bar + HeaderH), (vx1, y_top_bar + HeaderH)),
			closed=True, facecolor=col, edgecolor="none", zorder=0))
		ax_bar.text((vx1 + vx2) / 2.0, y_top_bar + HeaderH + 0.04, name,
			ha="center", va="bottom", fontsize=7.0, color="black",
			zorder=8, clip_on=False)

	# Trapezoids
	for mi, name in enumerate(locus_names):
		have = [lab for lab in by_y if tip_motifs[lab][mi] is not None]
		for i in range(len(have) - 1):
			lab_a = have[i]
			lab_b = have[i + 1]
			mot_a = tip_motifs[lab_a][mi]
			mot_b = tip_motifs[lab_b][mi]
			cen_a = AnchorOf(lab_a)
			cen_b = AnchorOf(lab_b)
			y_top = tip_y[lab_a] - TipBarH(lab_a) / 2.0
			y_bot = tip_y[lab_b] + TipBarH(lab_b) / 2.0
			clip_a = ClipSeg(*SegX(mot_a[0], mot_a[1], cen_a), XViewLo, XViewHi)
			clip_b = ClipSeg(*SegX(mot_b[0], mot_b[1], cen_b), XViewLo, XViewHi)
			if clip_a is None or clip_b is None:
				continue
			x1_a, x2_a = clip_a
			x1_b, x2_b = clip_b
			col = locus_colors[name]
			ax_bar.add_patch(Polygon(
				((x1_a, y_top), (x2_a, y_top), (x2_b, y_bot), (x1_b, y_bot)),
				closed=True, facecolor=col, edgecolor="none", zorder=0))

	def DrawLaneDomains(center, doms, db, y_lo, y_hi):
		text_ok = ContainedTextMask(doms)
		for (sf, start, end), draw_text in zip(doms, text_ok):
			clip_d = ClipSeg(*SegX(start, end, center), XViewLo, XViewHi)
			if clip_d is None:
				continue
			dvx1, dvx2 = clip_d
			used_sf[(db, sf)] = True
			ax_bar.add_patch(Rectangle((dvx1, y_lo), dvx2 - dvx1, y_hi - y_lo,
				facecolor=SfColors[(db, sf)], edgecolor="none", zorder=2))
			if draw_text:
				ax_bar.text((dvx1 + dvx2) / 2.0, (y_lo + y_hi) / 2.0, sf,
					ha="center", va="center", fontsize=SfFs, fontweight="bold",
					color="white", zorder=3, clip_on=True)

	for label in tip_order:
		y = tip_y[label]
		cls = ClassOf(disp_labels[label])
		col = LabelColors[cls]
		mA, mB, mC, qlen = tip_motifs[label]
		center = AnchorOf(label)
		lanes = TipLanes(label)
		bh = TipBarH(label)
		lh = TipLaneH(label)
		x1, x2 = SegX(1, qlen, center)
		trunc_left = x1 < XViewLo - 1e-9
		trunc_right = x2 > XViewHi + 1e-9
		clip = ClipSeg(x1, x2, XViewLo, XViewHi)
		if clip is None:
			continue
		vx1, vx2 = clip
		y0 = y - bh / 2.0
		y1 = y + bh / 2.0
		fill_x1 = vx1
		fill_x2 = vx2
		if trunc_left:
			fill_x1 = min(vx1 + TruncTriW, vx2)
		if trunc_right:
			fill_x2 = max(vx2 - TruncTriW, fill_x1)
		edge = to_rgba(col, 0.55)
		bar_fill = to_rgba(col, FillAlpha)
		if fill_x2 > fill_x1:
			ax_bar.add_patch(Rectangle((fill_x1, y0), fill_x2 - fill_x1, bh,
				facecolor=bar_fill, edgecolor="none", zorder=1))
		# Top/bottom edges stop at triangle base (not apex).
		ax_bar.plot([fill_x1, fill_x2], [y0, y0], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([fill_x1, fill_x2], [y1, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		if trunc_left:
			apex = (vx1, 0.5 * (y0 + y1))
			base_lo = (vx1 + TruncTriW, y0)
			base_hi = (vx1 + TruncTriW, y1)
			ax_bar.add_patch(Polygon(
				(base_lo, base_hi, apex),
				closed=True, facecolor=bar_fill, edgecolor="none", zorder=1))
			ax_bar.plot([base_lo[0], apex[0]], [base_lo[1], apex[1]],
				color=edge, lw=0.6, zorder=4, solid_capstyle="butt")
			ax_bar.plot([base_hi[0], apex[0]], [base_hi[1], apex[1]],
				color=edge, lw=0.6, zorder=4, solid_capstyle="butt")
		else:
			ax_bar.plot([vx1, vx1], [y0, y1], color=edge, lw=0.6, zorder=4,
				solid_capstyle="butt")
		if trunc_right:
			apex = (vx2, 0.5 * (y0 + y1))
			base_lo = (vx2 - TruncTriW, y0)
			base_hi = (vx2 - TruncTriW, y1)
			ax_bar.add_patch(Polygon(
				(base_lo, base_hi, apex),
				closed=True, facecolor=bar_fill, edgecolor="none", zorder=1))
			ax_bar.plot([base_lo[0], apex[0]], [base_lo[1], apex[1]],
				color=edge, lw=0.6, zorder=4, solid_capstyle="butt")
			ax_bar.plot([base_hi[0], apex[0]], [base_hi[1], apex[1]],
				color=edge, lw=0.6, zorder=4, solid_capstyle="butt")
		else:
			ax_bar.plot([vx2, vx2], [y0, y1], color=edge, lw=0.6, zorder=4,
				solid_capstyle="butt")
		for li in range(1, len(lanes)):
			yy = y1 - li * lh
			ax_bar.plot([fill_x1, fill_x2], [yy, yy], color=to_rgba(col, 0.25),
				lw=0.5, zorder=4, solid_capstyle="butt", linestyle=":")
		for li, db in enumerate(lanes):
			lane_hi = y1 - li * lh
			lane_lo = lane_hi - lh
			DrawLaneDomains(center, tip_domains[label][db], db, lane_lo, lane_hi)

	for label in tip_order:
		y = tip_y[label]
		bh = TipBarH(label)
		mA, mB, mC, qlen = tip_motifs[label]
		center = AnchorOf(label)
		for name, mot in zip(locus_names, (mA, mB, mC)):
			if mot is None:
				continue
			clip = ClipSeg(*SegX(mot[0], mot[1], center), XViewLo, XViewHi)
			if clip is None:
				continue
			vx1, vx2 = clip
			ax_bar.add_patch(Rectangle((vx1, y - bh / 2.0), vx2 - vx1, bh,
				facecolor=to_rgba(locus_colors[name], LocusBarAlpha),
				edgecolor="none", zorder=5))

	ax_bar.set_xlim(XViewLo, XViewHi)
	# Renumbered ticks: cropped window labeled 1..(CropTickHi-CropTickLo).
	shown = [1]
	k = 100
	tick_hi = CropTickHi - CropTickLo  # 800
	while k <= tick_hi + 0.01:
		shown.append(k)
		k += 100
	ax_bar.set_xticks([XViewLo + float(s) - 1.0 for s in shown])
	ax_bar.set_xticklabels(["%d" % int(s) for s in shown])
	ax_bar.tick_params(axis="y", left=False, labelleft=False)
	ax_bar.spines["left"].set_visible(False)
	ax_bar.spines["right"].set_visible(False)
	ax_bar.spines["top"].set_visible(False)
	ax_bar.tick_params(axis="x", labelsize=8, pad=2)

	fig.canvas.draw()
	renderer = fig.canvas.get_renderer()
	tick_bbs = [t.get_window_extent(renderer=renderer) for t in ax_bar.get_xticklabels()
		if t.get_text() != ""]
	inv = fig.transFigure.inverted()
	if len(tick_bbs) > 0:
		ax_bb = ax_bar.get_window_extent(renderer=renderer)
		cx = 0.5 * (ax_bb.x0 + ax_bb.x1)
		tick_bottom = min([b.y0 for b in tick_bbs])
		fx, fy = inv.transform((cx, tick_bottom - 1.0))
		xlab = fig.text(fx, fy, "amino acids", fontsize=9, ha="center", va="top",
			transform=fig.transFigure)
		fig.canvas.draw()
		renderer = fig.canvas.get_renderer()
		lab_bottom = xlab.get_window_extent(renderer=renderer).y0
		_, leg_top = inv.transform((0.0, lab_bottom - 6.0))
	else:
		leg_top = ax_y_frac - 0.02

	visible_groups = []
	for db, sfs in leg_groups:
		items = []
		for sf in sfs:
			if (db, sf) not in used_sf:
				continue
			items.append((db, sf, LegText(db, sf)))
		if len(items) > 0:
			visible_groups.append((db, items))

	if len(visible_groups) > 0:
		# Legend spacing in inches so tall figures (many rdrp tips) do not
		# stretch swatches / row gaps in figure-fraction coordinates.
		n_vis_items = sum(len(items) for _t, items in visible_groups)
		n_vis_hdrs = len(visible_groups)
		row_dy_in = 0.12
		hdr_dy_in = 0.15
		gap_grp_in = 0.06
		swatch_h_in = 0.07
		swatch_w = 0.010
		leg_need_in = (0.25 + n_vis_hdrs * hdr_dy_in + n_vis_items * row_dy_in
			+ gap_grp_in * max(0, n_vis_hdrs - 1))
		# Space from xlab bottom to figure bottom, in inches.
		leg_top_in = leg_top * fig_h
		avail_in = max(0.5, leg_top_in - 0.15)
		if leg_need_in > avail_in:
			extra_in = leg_need_in - avail_in + 0.35
			fig_h += extra_in
			bottom_in += extra_in
			fig.set_size_inches(fig_w, fig_h)
			ax_h_frac = 1.0 - (bottom_in + top_in) / fig_h
			ax_y_frac = bottom_in / fig_h
			ax_bar.set_position([(xoff + bar_left_in) / fig_w, ax_y_frac,
				bar_ax_in / fig_w, ax_h_frac])
			fig.canvas.draw()
			renderer = fig.canvas.get_renderer()
			tick_bbs = [t.get_window_extent(renderer=renderer)
				for t in ax_bar.get_xticklabels() if t.get_text() != ""]
			inv = fig.transFigure.inverted()
			if len(tick_bbs) > 0:
				ax_bb = ax_bar.get_window_extent(renderer=renderer)
				cx = 0.5 * (ax_bb.x0 + ax_bb.x1)
				tick_bottom = min([b.y0 for b in tick_bbs])
				fx, fy = inv.transform((cx, tick_bottom - 1.0))
				for txt in list(fig.texts):
					if getattr(txt, "get_text", lambda: "")() == "amino acids":
						txt.remove()
				xlab = fig.text(fx, fy, "amino acids", fontsize=9, ha="center",
					va="top", transform=fig.transFigure)
				fig.canvas.draw()
				renderer = fig.canvas.get_renderer()
				lab_bottom = xlab.get_window_extent(renderer=renderer).y0
				_, leg_top = inv.transform((0.0, lab_bottom - 6.0))
		row_dy = row_dy_in / fig_h
		hdr_dy = hdr_dy_in / fig_h
		swatch_h = swatch_h_in / fig_h
		gap = 0.008
		# Align legend left edge with the GDD column start.
		ax_left_fig = (xoff + bar_left_in) / fig_w
		gdd_fig = ax_left_fig + (gdd_x - XViewLo) / x_span * (bar_ax_in / fig_w)
		left_x = gdd_fig
		y = leg_top
		for gi, (title, items) in enumerate(visible_groups):
			if gi > 0:
				y -= gap_grp_in / fig_h
			fig.text(left_x, y, title, fontsize=LegFs + 0.5, fontweight="bold",
				ha="left", va="top", transform=fig.transFigure, color="#333333")
			y -= hdr_dy
			for db, sf, txt in items:
				# Align swatch to text midline (va=top text).
				fig.patches.append(Rectangle(
					(left_x, y - 0.70 * swatch_h), swatch_w, swatch_h,
					transform=fig.transFigure, facecolor=SfColors[(db, sf)],
					edgecolor="none", clip_on=False))
				fig.text(left_x + swatch_w + gap, y, txt, fontsize=LegFs,
					ha="left", va="top", transform=fig.transFigure)
				y -= row_dy

	for OutPath in out_paths:
		outdir = os.path.dirname(OutPath)
		if outdir and not os.path.isdir(outdir):
			os.makedirs(outdir)
		fig.savefig(OutPath, dpi=200)
		sys.stderr.write("wrote %s\n" % OutPath)
	plt.close(fig)


if __name__ == "__main__":
	Main()
