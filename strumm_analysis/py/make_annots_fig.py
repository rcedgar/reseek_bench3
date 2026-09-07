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

def CmapHex(cmap_name, i, n, t0=0.30, t1=0.85):
	cmap = plt.get_cmap(cmap_name)
	if n <= 1:
		t = 0.5 * (t0 + t1)
	else:
		t = t0 + (t1 - t0) * (float(i) / float(max(1, n - 1)))
	r, g, b, _a = cmap(t)
	return "#%02X%02X%02X" % (int(round(r * 255)), int(round(g * 255)),
		int(round(b * 255)))

def PfamColor(i, n):
	return CmapHex("Greys", i, n, t0=0.45, t1=0.80)

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
	}


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
			domains[lab][lane].sort(key=lambda d: (d[1], d[2], d[0]))
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

	# Full-length bars: viewport covers every tip end-to-end (no truncation).
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
	ax_bar.set_position([(xoff + bar_left_in) / fig_w, ax_y_frac,
		bar_ax_in / fig_w, ax_h_frac])

	# Match label↔bar gap to vertical BarGap in physical (inch) units.
	ax_h_in = fig_h * ax_h_frac
	gap_in = BarGap / data_span * ax_h_in
	label_gap_x = gap_in * (XViewHi - XViewLo) / bar_ax_in

	for label in tip_order:
		cls = ClassOf(disp_labels[label])
		col = LabelColors[cls]
		center = AnchorOf(label)
		mA, mB, mC, qlen = tip_motifs[label]
		bar_x1, _ = SegX(1, qlen, center)
		ax_bar.text(bar_x1 - label_gap_x, tip_y[label], disp_labels[label],
			va="center", ha="right", fontsize=LabelFs, fontfamily="monospace",
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
		vx1, vx2 = x1, x2
		y0 = y - bh / 2.0
		y1 = y + bh / 2.0
		ax_bar.add_patch(Rectangle((vx1, y0), vx2 - vx1, bh,
			facecolor=to_rgba(col, FillAlpha), edgecolor="none", zorder=1))
		edge = to_rgba(col, 0.55)
		ax_bar.plot([vx1, vx2], [y0, y0], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([vx1, vx2], [y1, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([vx1, vx1], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([vx2, vx2], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		for li in range(1, len(lanes)):
			yy = y1 - li * lh
			ax_bar.plot([vx1, vx2], [yy, yy], color=to_rgba(col, 0.25),
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
	shown = [1]
	k = 100
	while XViewLo + float(k) - 1.0 <= XViewHi + 0.01:
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
		left_x = 0.06
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
