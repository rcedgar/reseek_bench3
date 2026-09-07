#!/usr/bin/python3

import argparse
import os
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

Usage = \
(
"Draw amino-acid-scale protein bars with Pfam, SCOP, and CATH domain fills "
"and motif trapezoid connectors. Each present database is one equal-height "
"lane (bar height 1x/2x/3x). Tip order comes from a figure-order file."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--scop-domains", required=True,
	help="SCOP domain annot TSV from hits_to_domain_annots.py")
AP.add_argument("--cath-domains", required=True,
	help="CATH domain annot TSV from hits_to_domain_annots.py")
AP.add_argument("--pfam-domains", default=None,
	help="Pfam domain TSV (label pfamid start end length); optional")
AP.add_argument("--pfam-names", default=None,
	help="Pfam accession + name file (e.g. pfams.txt)")
AP.add_argument("--motifs", required=True,
	help="Motif TSV (rdrp_motifs.tsv or cdn_motifs.tsv)")
AP.add_argument("--fasta", required=True, action="append",
	help="FASTA matching domain coordinates (repeatable; first wins on key clash)")
AP.add_argument("--motif-style", choices=("abc", "iii"), required=True,
	help="abc: A/B/C from rdrp motifs; iii: I/II/III from cdn motifs")
AP.add_argument("--order", required=True,
	help="Tip order file (rdrp_figure_order.txt or cdn_figure_order.txt)")
AP.add_argument("--scop-names", required=True,
	help="SCOP superfamily name table")
AP.add_argument("--cath-names", required=True,
	help="CATH superfamily name table")
AP.add_argument("--output", required=True, help="Output PDF or PNG path")
AP.add_argument("--view-lo", type=int, default=1,
	help="Left residue coordinate for viewport (default 1)")
Args = AP.parse_args()

ResViewLo = Args.view_lo

def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)

def MotifOrNone(s):
	if s is None:
		return None
	t = s.strip().upper()
	if t == "" or t == ".":
		return None
	return t

def StripChainize(label):
	s = label.strip()
	for suf in ("-rdrp", "-decoy"):
		i = s.rfind(suf)
		if i >= 0:
			rest = s[i + len(suf):]
			if len(rest) == 2 and rest[0] == "_" and rest[1].isalnum():
				s = s[: i + len(suf)]
			break
	return s

def NormKey(label):
	return StripChainize(label).upper().replace("_", "-")

def DisplayLabel(label):
	s = StripChainize(label)
	for suf in ("-rdrp", "-decoy"):
		i = s.rfind(suf)
		if i >= 0:
			s = s[:i]
			break
	return s

def ClassOf(label):
	s = StripChainize(label)
	if s.endswith("-rdrp"):
		return "rdrp"
	if s.endswith("-decoy"):
		return "decoy"
	return "other"

def OrderIdOf(label, motif_style):
	disp = DisplayLabel(label)
	if motif_style == "iii":
		for sep in ("-", "_"):
			if sep in disp:
				return disp.split(sep, 1)[0].upper()
		return disp.upper()
	return disp

def FormatSf(sf, cath):
	if cath:
		return sf.replace("_", ".")
	return sf

# Shared hues for paired SCOP/CATH; SCOP = cooler/higher sat, CATH = warmer/lower sat.
PairHues = [
	0.58, 0.48, 0.68, 0.40, 0.62, 0.52, 0.72, 0.35,
	0.55, 0.45, 0.65, 0.38,
]
UnpairedScopHues = [0.60, 0.55, 0.70, 0.50, 0.65, 0.45]
UnpairedCathHues = [0.08, 0.02, 0.12, 0.95, 0.05, 0.15]

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

# One equal-height lane per present DB (Pfam / SCOP / CATH) → bar 1x, 2x, or 3x.
LaneH = 0.54
BarGap = LaneH * 0.20
SfFs = 5.5
LabelFs = 8.0
LegFs = 6.5
FillAlpha = 0.10
LocusBarAlpha = 0.20
HeaderH = 0.32
AbcPad = 30
LeftExtra = 50
SfHiPad = 20
OverlapFracMin = 0.50

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
	# Darker greys so white labels stay readable.
	return CmapHex("Greys", i, n, t0=0.45, t1=0.80)

def ScopColor(i, n):
	return CmapHex("RdPu", i, n, t0=0.35, t1=0.85)

def CathColor(i, n):
	return CmapHex("GnBu", i, n, t0=0.35, t1=0.85)

def ReadFasta(fns):
	# Keep every distinct label (needed so Pfam full-length coords can map
	# onto palmcore tips). by_key still first-wins for tip selection.
	seqs = {}
	by_key = {}
	for fn in fns:
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
		Die("no sequences in fasta")
	return seqs, by_key

def ReadMotifs(fn, style):
	annot = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] in ("Label", "label"):
			continue
		label = Fields[0]
		key = NormKey(label)
		if style == "abc":
			if len(Fields) < 6:
				Die("bad rdrp motif row: %s" % Line[:80])
			vals = (MotifOrNone(Fields[3]), MotifOrNone(Fields[4]),
				MotifOrNone(Fields[5]))
			names = ("A", "B", "C")
		else:
			if len(Fields) < 4:
				Die("bad cdn motif row: %s" % Line[:80])
			vals = (MotifOrNone(Fields[1]), MotifOrNone(Fields[2]),
				MotifOrNone(Fields[3]))
			names = ("I", "II", "III")
		if key in annot:
			Die("duplicate motif key %s (%s)" % (key, label))
		annot[key] = (names, vals)
	return annot

def ReadDomains(fn, cath):
	by_query = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] == "query":
			continue
		if len(Fields) < 5:
			Die("bad domain row: %s" % Line[:80])
		query = Fields[0]
		sf = FormatSf(Fields[1], cath)
		try:
			env_lo = int(Fields[2])
			env_hi = int(Fields[3])
		except ValueError:
			Die("bad domain coords: %s" % Line[:80])
		if env_lo < 1 or env_hi < env_lo:
			Die("bad domain coords %s %s %d-%d" % (query, sf, env_lo, env_hi))
		by_query.setdefault(query, []).append((sf, env_lo, env_hi))
	return by_query

def ReadPfamDomains(fn):
	"""label pfamid start end length → by NormKey: list of (acc, lo, hi, length)."""
	by_key = {}
	if fn is None or not os.path.isfile(fn):
		return by_key
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] in ("label", "Label"):
			continue
		if len(Fields) < 5:
			Die("bad pfam domain row: %s" % Line[:80])
		label = Fields[0]
		acc = Fields[1]
		try:
			start = int(Fields[2])
			end = int(Fields[3])
			length = int(Fields[4])
		except ValueError:
			Die("bad pfam coords: %s" % Line[:80])
		key = NormKey(label)
		by_key.setdefault(key, []).append((acc, start, end, length))
	return by_key

def ReadPfamNames(fn):
	"""pfams.txt (acc name) or hmmscan domain tbl (accession + trailing description)."""
	names = {}
	if fn is None or not os.path.isfile(fn):
		return names
	# Detect hmmscan domain tabular by comment header.
	is_tbl = False
	with open(fn) as fh:
		for Line in fh:
			if Line.startswith("#"):
				if "target name" in Line and "accession" in Line:
					is_tbl = True
					break
			elif len(Line.strip()) > 0:
				break
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0 or Line.startswith("#"):
			continue
		if is_tbl:
			parts = Line.split()
			if len(parts) < 23:
				continue
			acc = parts[1].split(".", 1)[0]
			desc = " ".join(parts[22:])
			if acc not in names:
				names[acc] = desc
		else:
			parts = Line.split(None, 1)
			acc = parts[0].split(".", 1)[0]
			names[acc] = parts[1].strip() if len(parts) > 1 else ""
	return names

def FindLenMatchedSeq(key, length, seqs):
	for lab, seq in seqs.items():
		if NormKey(lab) == key and len(seq) == length:
			return seq
	return None

def ProjectPfamToTip(pfam_rows, tip_seq, seqs, key):
	"""Map full-length Pfam coords onto tip sequence coordinates."""
	if len(pfam_rows) == 0:
		return []
	length = pfam_rows[0][3]
	full = FindLenMatchedSeq(key, length, seqs)
	if full is None:
		# Tip itself may be the annotated length.
		if len(tip_seq) == length:
			full = tip_seq
		else:
			return []
	if tip_seq == full:
		return [(acc, lo, hi) for acc, lo, hi, _L in pfam_rows]
	idx = full.find(tip_seq)
	if idx < 0:
		return []
	off = idx
	tip_len = len(tip_seq)
	out = []
	for acc, lo, hi, _L in pfam_rows:
		lo2 = lo - off
		hi2 = hi - off
		if hi2 < 1 or lo2 > tip_len:
			continue
		lo2 = max(1, lo2)
		hi2 = min(tip_len, hi2)
		if lo2 <= hi2:
			out.append((acc, lo2, hi2))
	return out

def ReadOrder(fn):
	order = []
	seen = set()
	for Line in open(fn):
		tok = Line.strip()
		if len(tok) == 0 or tok.startswith("#"):
			continue
		key = tok.upper() if Args.motif_style == "iii" else tok
		if key in seen:
			Die("duplicate order id %s" % tok)
		seen.add(key)
		order.append(key)
	if len(order) == 0:
		Die("empty order file %s" % fn)
	return order

def ReadCathNames(fn):
	names = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		parts = Line.split(None, 1)
		sf = parts[0]
		if len(parts) == 1:
			names[sf] = ""
		else:
			names[sf] = parts[1].strip()
	return names

def ReadScopNames(fn):
	names = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if len(Fields) < 5:
			continue
		if Fields[1] != "sf":
			continue
		sf = Fields[2]
		name = Fields[4].strip()
		names[sf] = name
	return names

def SfDisplayName(sf, name_map):
	if sf not in name_map:
		return "(missing)"
	nm = name_map[sf]
	if nm == "":
		return "(unnamed)"
	return nm

def FindMotif(seq, motif):
	if motif is None:
		return None
	hits = []
	start = 0
	while True:
		i = seq.find(motif, start)
		if i < 0:
			break
		hits.append(i)
		start = i + 1
	if len(hits) != 1:
		return None
	lo = hits[0] + 1
	hi = lo + len(motif) - 1
	return lo, hi

def MotifCenter(m_mid, qlen):
	if m_mid is None:
		return (1.0 + float(qlen)) / 2.0
	return (m_mid[0] + m_mid[1]) / 2.0

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

def IntervalOverlapFrac(a_lo, a_hi, b_lo, b_hi):
	"""Intersection / min(length); 0 if either empty."""
	inter = max(0, min(a_hi, b_hi) - max(a_lo, b_lo) + 1)
	if inter == 0:
		return 0.0
	la = a_hi - a_lo + 1
	lb = b_hi - b_lo + 1
	m = min(la, lb)
	if m <= 0:
		return 0.0
	return float(inter) / float(m)

def BestPairOverlap(sdoms, cdoms, scop_sf, cath_sf):
	best = 0.0
	for sf, slo, shi in sdoms:
		if sf != scop_sf:
			continue
		for cf, clo, chi in cdoms:
			if cf != cath_sf:
				continue
			best = max(best, IntervalOverlapFrac(slo, shi, clo, chi))
	return best

def PairScopCath(tip_order, scop_domains, cath_domains, scop_sfs, cath_sfs):
	"""
	Pair SCOP/CATH superfamilies when, in every tip where both occur, some
	domain interval pair overlaps by >50% of the shorter span.
	"""
	co_ok = {}
	co_n = {}
	for ssf in scop_sfs:
		for csf in cath_sfs:
			oks = []
			for label in tip_order:
				sdoms = [d for d in scop_domains[label] if d[0] == ssf]
				cdoms = [d for d in cath_domains[label] if d[0] == csf]
				if len(sdoms) == 0 or len(cdoms) == 0:
					continue
				oks.append(BestPairOverlap(scop_domains[label],
					cath_domains[label], ssf, csf) > OverlapFracMin)
			if len(oks) == 0:
				continue
			co_n[(ssf, csf)] = len(oks)
			co_ok[(ssf, csf)] = all(oks)

	candidates = []
	for (ssf, csf), ok in co_ok.items():
		if not ok:
			continue
		candidates.append((-co_n[(ssf, csf)], ssf, csf))
	candidates.sort()

	paired = []
	used_s = set()
	used_c = set()
	for _n, ssf, csf in candidates:
		if ssf in used_s or csf in used_c:
			continue
		used_s.add(ssf)
		used_c.add(csf)
		paired.append((ssf, csf))

	unpaired_s = [s for s in scop_sfs if s not in used_s]
	unpaired_c = [c for c in cath_sfs if c not in used_c]
	return paired, unpaired_s, unpaired_c

def IndexByOrderId(by_query, motif_style):
	out = {}
	for query in by_query:
		oid = OrderIdOf(query, motif_style)
		if motif_style == "iii":
			oid = oid.upper()
		out.setdefault(oid, []).append(query)
	return out

def PickQuery(dom_cands):
	for q in dom_cands:
		if NormKey(q) in fasta_by_key:
			return q
	return None

def ResolveDoms(query, domains_raw, qlen, tag):
	if query is None:
		return []
	doms = sorted(domains_raw[query], key=lambda d: (d[1], d[2], d[0]))
	for sf, lo, hi in doms:
		if hi > qlen:
			Die("domain past end %s %s %s %d>%d" %
				(tag, query, sf, hi, qlen))
	return doms

def DrawTruncArrow(ax, x_edge, y, inward, col, arrow_w, bar_h):
	y0 = y - bar_h / 2.0
	y1 = y + bar_h / 2.0
	tip_x = x_edge
	base_x = x_edge + inward * arrow_w
	ax.add_patch(Polygon(
		((tip_x, y), (base_x, y0), (base_x, y1)),
		closed=True, facecolor=to_rgba(col, FillAlpha), edgecolor="none",
		zorder=6))

seqs, fasta_by_key = ReadFasta(Args.fasta)
motif_annot = ReadMotifs(Args.motifs, Args.motif_style)
scop_raw = ReadDomains(Args.scop_domains, cath=False)
cath_raw = ReadDomains(Args.cath_domains, cath=True)
pfam_raw = ReadPfamDomains(Args.pfam_domains)
order_ids = ReadOrder(Args.order)
scop_name_map = ReadScopNames(Args.scop_names)
cath_name_map = ReadCathNames(Args.cath_names)
pfam_name_map = ReadPfamNames(Args.pfam_names)

if Args.motif_style == "abc":
	LocusColors = LocusColorsAbc
	LocusNames = ("A", "B", "C")
else:
	LocusColors = LocusColorsIii
	LocusNames = ("I", "II", "III")

scop_by_oid = IndexByOrderId(scop_raw, Args.motif_style)
cath_by_oid = IndexByOrderId(cath_raw, Args.motif_style)

fa_by_oid = {}
for key, fa_label in fasta_by_key.items():
	oid = OrderIdOf(fa_label, Args.motif_style)
	if Args.motif_style == "iii":
		oid = oid.upper()
	fa_by_oid.setdefault(oid, []).append(fa_label)

motifs = {}
scop_domains = {}
cath_domains = {}
pfam_domains = {}
tip_order = []
for oid in order_ids:
	fa_cands = fa_by_oid.get(oid, [])
	scop_q = PickQuery(scop_by_oid.get(oid, []))
	cath_q = PickQuery(cath_by_oid.get(oid, []))
	chosen = None
	if scop_q is not None:
		chosen = fasta_by_key[NormKey(scop_q)]
	elif cath_q is not None:
		chosen = fasta_by_key[NormKey(cath_q)]
	elif len(fa_cands) > 0:
		chosen = fa_cands[0]
	else:
		sys.stderr.write("skip no fasta for order id %s\n" % oid)
		continue

	key = NormKey(chosen)
	seq = seqs[chosen]
	qlen = len(seq)
	sdoms = ResolveDoms(scop_q, scop_raw, qlen, "scop")
	cdoms = ResolveDoms(cath_q, cath_raw, qlen, "cath")
	pdoms = ProjectPfamToTip(pfam_raw.get(key, []), seq, seqs, key)

	spans = [None, None, None]
	if key in motif_annot:
		names, vals = motif_annot[key]
		spans = []
		for name, pep in zip(names, vals):
			span = FindMotif(seq, pep)
			if pep is not None and span is None:
				sys.stderr.write("warn motif %s not unique in %s\n" %
					(name, chosen))
			spans.append(span)
		while len(spans) < 3:
			spans.append(None)

	tip_order.append(chosen)
	motifs[chosen] = (spans[0], spans[1], spans[2], qlen)
	scop_domains[chosen] = sdoms
	cath_domains[chosen] = cdoms
	pfam_domains[chosen] = pdoms

if len(tip_order) == 0:
	Die("no tips to draw")
n_leaf = len(tip_order)

def TipLanes(label):
	"""Top→bottom lane DBs present on this tip."""
	lanes = []
	if len(pfam_domains.get(label, [])) > 0:
		lanes.append("pfam")
	if len(scop_domains.get(label, [])) > 0:
		lanes.append("scop")
	if len(cath_domains.get(label, [])) > 0:
		lanes.append("cath")
	return lanes

def TipBarH(label):
	n = len(TipLanes(label))
	if n == 0:
		return LaneH
	return float(n) * LaneH

# Stack tips top→bottom with uniform edge gaps.
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

pfam_sf_order = []
scop_sf_order = []
cath_sf_order = []
sf_seen = set()
for label in tip_order:
	for sf in SfArchitecture(pfam_domains[label]):
		key = ("pfam", sf)
		if key not in sf_seen:
			sf_seen.add(key)
			pfam_sf_order.append(sf)
	for sf in SfArchitecture(scop_domains[label]):
		key = ("scop", sf)
		if key not in sf_seen:
			sf_seen.add(key)
			scop_sf_order.append(sf)
	for sf in SfArchitecture(cath_domains[label]):
		key = ("cath", sf)
		if key not in sf_seen:
			sf_seen.add(key)
			cath_sf_order.append(sf)

SfColors = {}
for i, sf in enumerate(pfam_sf_order):
	SfColors[("pfam", sf)] = PfamColor(i, max(1, len(pfam_sf_order)))
for i, sf in enumerate(scop_sf_order):
	SfColors[("scop", sf)] = ScopColor(i, max(1, len(scop_sf_order)))
for i, sf in enumerate(cath_sf_order):
	SfColors[("cath", sf)] = CathColor(i, max(1, len(cath_sf_order)))

x0 = None
for label in tip_order:
	mA, mB, mC, qlen = motifs[label]
	center = MotifCenter(mB, qlen)
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

ref_label = tip_order[0]
mA_ref, mB_ref, mC_ref, qlen_ref = motifs[ref_label]
ref_center = MotifCenter(mB_ref, qlen_ref)
XViewLo = XOf(ResViewLo, ref_center) - 0.5
XViewHi = None
for label in tip_order:
	mA, mB, mC, qlen = motifs[label]
	center = MotifCenter(mB, qlen)
	if mA is not None:
		lo = max(1, mA[0] - AbcPad)
		ax1 = XOf(lo, center) - 0.5
		if ax1 < XViewLo:
			XViewLo = ax1
	for doms in (pfam_domains[label], scop_domains[label], cath_domains[label]):
		for sf, start, end in doms:
			hx = XOf(end + SfHiPad, center) + 0.5
			if XViewHi is None or hx > XViewHi:
				XViewHi = hx
	hx = XOf(qlen, center) + 0.5
	if XViewHi is None or hx > XViewHi:
		XViewHi = hx
XViewLo -= LeftExtra
if XViewHi <= XViewLo:
	Die("empty bar x range")

disp_labels = {}
for label in tip_order:
	disp_labels[label] = DisplayLabel(label)

fig_w = 12.0
y_top_extent = HeaderH + 1.05
y_bot_extent = y_bottom - 0.20
right_in = 0.18
# Reserve enough figure bottom for single-column Pfam/SCOP/CATH legend.
n_leg_items = (len(pfam_sf_order) + len(scop_sf_order) + len(cath_sf_order))
n_leg_hdrs = sum(1 for s in (pfam_sf_order, scop_sf_order, cath_sf_order) if len(s) > 0)
leg_est_in = 0.55 + n_leg_hdrs * 0.22 + n_leg_items * 0.18 + 0.12 * max(0, n_leg_hdrs - 1)
bottom_in = max(3.6, leg_est_in + 0.90)
top_in = 0.50
data_span = y_top_extent - y_bot_extent
fig_h = max(8.5, data_span * 0.38 + 5.0 + max(0.0, bottom_in - 3.6))
ax_h_frac = 1.0 - (bottom_in + top_in) / fig_h
ax_y_frac = bottom_in / fig_h

fig = plt.figure(figsize=(fig_w, fig_h))
ax_bar = fig.add_axes([0.20, ax_y_frac, 0.74, ax_h_frac])
ax_bar.set_ylim(y_bot_extent, y_top_extent)

def ProbeIn(s, fs=LabelFs, family="monospace"):
	t = fig.text(0, 0, s, fontsize=fs, fontfamily=family)
	fig.canvas.draw()
	w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
	t.remove()
	return w

upper_in = ProbeIn("MM") - ProbeIn("M")

bar_ax_in = fig_w * 0.74
bar_left_in = fig_w * 0.20
plot_w = bar_left_in + bar_ax_in + right_in

# Single-column legend grouped Pfam → SCOP → CATH.
def LegText(db, sf):
	if db == "pfam":
		return "%s %s" % (sf, SfDisplayName(sf, pfam_name_map))
	if db == "scop":
		return "%s %s" % (sf, SfDisplayName(sf, scop_name_map))
	return "%s %s" % (sf, SfDisplayName(sf, cath_name_map))

leg_groups = [
	("Pfam", "pfam", pfam_sf_order),
	("SCOP", "scop", scop_sf_order),
	("CATH", "cath", cath_sf_order),
]
longest_leg = "Pfam"
for _title, db, sfs in leg_groups:
	for sf in sfs:
		txt = LegText(db, sf)
		if len(txt) > len(longest_leg):
			longest_leg = txt
leg_txt_in = ProbeIn(longest_leg, fs=LegFs, family="sans-serif") if longest_leg else 0.0
needed_leg_in = 0.28 + 0.14 + leg_txt_in + 0.50
fig_w = max(plot_w, needed_leg_in + 0.40)
xoff = 0.5 * (fig_w - plot_w)
fig.set_size_inches(fig_w, fig_h)
bar_ax_in = fig_w - (bar_left_in + right_in)
ax_bar.set_position([(xoff + bar_left_in) / fig_w, ax_y_frac,
	bar_ax_in / fig_w, ax_h_frac])
ArrowW = 1.5 * upper_in * (XViewHi - XViewLo) / bar_ax_in

label_dx = 0.015 * (XViewHi - XViewLo)
for label in tip_order:
	cls = ClassOf(label)
	col = LabelColors[cls]
	ax_bar.text(XViewLo - label_dx, tip_y[label], disp_labels[label],
		va="center", ha="right", fontsize=LabelFs, fontfamily="monospace",
		color=col, clip_on=False, zorder=2,
		bbox=dict(boxstyle="square,pad=0.12", facecolor=to_rgba(col, FillAlpha),
			edgecolor="none"))

used_sf = {}
cut_marks = []
# tip_order is already top→bottom; tip_y decreases downward.
by_y = list(tip_order)

# Header badges: first tip from the top that has each motif.
for mi, name in enumerate(LocusNames):
	top_lab = None
	for lab in by_y:
		if motifs[lab][mi] is not None:
			top_lab = lab
			break
	if top_lab is None:
		continue
	mot = motifs[top_lab][mi]
	cen_top = MotifCenter(motifs[top_lab][1], motifs[top_lab][3])
	y_top_bar = tip_y[top_lab] + TipBarH(top_lab) / 2.0
	clip = ClipSeg(*SegX(mot[0], mot[1], cen_top), XViewLo, XViewHi)
	if clip is None:
		continue
	vx1, vx2 = clip
	col = LocusColors[name]
	ax_bar.add_patch(Polygon(
		((vx1, y_top_bar), (vx2, y_top_bar),
			(vx2, y_top_bar + HeaderH), (vx1, y_top_bar + HeaderH)),
		closed=True, facecolor=col, edgecolor="none", zorder=0))
	ax_bar.text((vx1 + vx2) / 2.0, y_top_bar + HeaderH + 0.04, name,
		ha="center", va="bottom", fontsize=7.0, color="black",
		zorder=8, clip_on=False)

# Trapezoids: for each motif, connect consecutive tips that have it
# (skip tips missing that motif so bands stay continuous).
for mi, name in enumerate(LocusNames):
	have = [lab for lab in by_y if motifs[lab][mi] is not None]
	for i in range(len(have) - 1):
		lab_a = have[i]
		lab_b = have[i + 1]
		mot_a = motifs[lab_a][mi]
		mot_b = motifs[lab_b][mi]
		cen_a = MotifCenter(motifs[lab_a][1], motifs[lab_a][3])
		cen_b = MotifCenter(motifs[lab_b][1], motifs[lab_b][3])
		y_top = tip_y[lab_a] - TipBarH(lab_a) / 2.0
		y_bot = tip_y[lab_b] + TipBarH(lab_b) / 2.0
		clip_a = ClipSeg(*SegX(mot_a[0], mot_a[1], cen_a), XViewLo, XViewHi)
		clip_b = ClipSeg(*SegX(mot_b[0], mot_b[1], cen_b), XViewLo, XViewHi)
		if clip_a is None or clip_b is None:
			continue
		x1_a, x2_a = clip_a
		x1_b, x2_b = clip_b
		col = LocusColors[name]
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

DomByDb = {
	"pfam": pfam_domains,
	"scop": scop_domains,
	"cath": cath_domains,
}

for label in tip_order:
	y = tip_y[label]
	cls = ClassOf(label)
	col = LabelColors[cls]
	mA, mB, mC, qlen = motifs[label]
	center = MotifCenter(mB, qlen)
	lanes = TipLanes(label)
	bh = TipBarH(label)
	x1, x2 = SegX(1, qlen, center)
	left_trunc = x1 < XViewLo
	right_trunc = x2 > XViewHi
	clip = ClipSeg(x1, x2, XViewLo, XViewHi)
	if clip is None:
		continue
	vx1, vx2 = clip
	y0 = y - bh / 2.0
	y1 = y + bh / 2.0
	draw1 = vx1 + ArrowW if left_trunc else vx1
	draw2 = vx2 - ArrowW if right_trunc else vx2
	if draw2 <= draw1:
		draw1, draw2 = vx1, vx2
		left_trunc = False
		right_trunc = False
	ax_bar.add_patch(Rectangle((draw1, y0), draw2 - draw1, bh,
		facecolor=to_rgba(col, FillAlpha), edgecolor="none", zorder=1))
	edge = to_rgba(col, 0.55)
	ax_bar.plot([draw1, draw2], [y0, y0], color=edge, lw=0.6, zorder=4,
		solid_capstyle="butt")
	ax_bar.plot([draw1, draw2], [y1, y1], color=edge, lw=0.6, zorder=4,
		solid_capstyle="butt")
	# Lane dividers when multiple DBs present.
	for li in range(1, len(lanes)):
		yy = y1 - li * LaneH
		ax_bar.plot([draw1, draw2], [yy, yy], color=to_rgba(col, 0.25),
			lw=0.5, zorder=4, solid_capstyle="butt", linestyle=":")
	if not left_trunc:
		ax_bar.plot([vx1, vx1], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
	if not right_trunc:
		ax_bar.plot([vx2, vx2], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
	for li, db in enumerate(lanes):
		lane_hi = y1 - li * LaneH
		lane_lo = lane_hi - LaneH
		DrawLaneDomains(center, DomByDb[db][label], db, lane_lo, lane_hi)
	if left_trunc:
		cut_marks.append((XViewLo, y, 1.0, col, bh))
	if right_trunc:
		cut_marks.append((XViewHi, y, -1.0, col, bh))

for label in tip_order:
	y = tip_y[label]
	bh = TipBarH(label)
	mA, mB, mC, qlen = motifs[label]
	center = MotifCenter(mB, qlen)
	for name, mot in zip(LocusNames, (mA, mB, mC)):
		if mot is None:
			continue
		clip = ClipSeg(*SegX(mot[0], mot[1], center), XViewLo, XViewHi)
		if clip is None:
			continue
		vx1, vx2 = clip
		ax_bar.add_patch(Rectangle((vx1, y - bh / 2.0), vx2 - vx1, bh,
			facecolor=to_rgba(LocusColors[name], LocusBarAlpha), edgecolor="none",
			zorder=5))

for x_edge, y, inward, col, bh in cut_marks:
	DrawTruncArrow(ax_bar, x_edge, y, inward, col, ArrowW, bh)

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
	leg_top = 0.14

# Single-column legend grouped Pfam / SCOP / CATH.
visible_groups = []
for title, db, sfs in leg_groups:
	items = []
	for sf in sfs:
		if (db, sf) not in used_sf:
			continue
		items.append((db, sf, LegText(db, sf)))
	if len(items) > 0:
		visible_groups.append((title, items))

if len(visible_groups) > 0:
	swatch_w = 0.012
	swatch_h = 0.012
	gap = 0.010
	row_dy = 0.016
	hdr_dy = 0.020
	left_x = 0.08
	hdr_y = leg_top
	y = hdr_y
	for gi, (title, items) in enumerate(visible_groups):
		if gi > 0:
			y -= 0.008
		fig.text(left_x, y, title, fontsize=LegFs + 0.5, fontweight="bold",
			ha="left", va="top", transform=fig.transFigure, color="#333333")
		y -= hdr_dy
		for db, sf, txt in items:
			fig.patches.append(Rectangle(
				(left_x, y - swatch_h * 0.85), swatch_w, swatch_h,
				transform=fig.transFigure, facecolor=SfColors[(db, sf)],
				edgecolor="none", clip_on=False))
			fig.text(left_x + swatch_w + gap, y, txt, fontsize=LegFs,
				ha="left", va="top", transform=fig.transFigure)
			y -= row_dy

OutPath = Args.output
fig.savefig(OutPath, dpi=200)
root, ext = os.path.splitext(OutPath)
if ext.lower() == ".pdf":
	fig.savefig(root + ".png", dpi=200)
elif ext.lower() == ".png":
	fig.savefig(root + ".pdf", dpi=200)
plt.close(fig)
sys.stderr.write("wrote %s\n" % OutPath)
