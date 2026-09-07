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
	from matplotlib.patches import Patch, Polygon, Rectangle
except ImportError:
	sys.stderr.write("matplotlib is required\n")
	sys.exit(1)

Usage = \
(
"Draw amino-acid-scale protein bars with SCOP/CATH superfamily domain fills "
"and motif trapezoid connectors. Tip order comes from a figure-order file; "
"sequences without domain annotations are still drawn."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--domains", required=True,
	help="Domain annot TSV from hits_to_domain_annots.py")
AP.add_argument("--motifs", required=True,
	help="Motif TSV (rdrp_motifs.tsv or cdn_motifs.tsv)")
AP.add_argument("--fasta", required=True, action="append",
	help="FASTA matching domain coordinates (repeatable; first wins on key clash)")
AP.add_argument("--motif-style", choices=("abc", "iii"), required=True,
	help="abc: A/B/C from rdrp motifs; iii: I/II/III from cdn motifs")
AP.add_argument("--order", required=True,
	help="Tip order file (rdrp_figure_order.txt or cdn_figure_order.txt)")
AP.add_argument("--sf-names", required=True,
	help="Superfamily name table (cath-b-newest-names or scop_superfamily_names.txt)")
AP.add_argument("--cath", action="store_true",
	help="CATH mode: replace '_' with '.' in SF ids; parse cath name file")
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
	"""Join query / motif / fasta labels across chainize and _/- variants."""
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

Palette = [
	"#0077BB", "#33BBEE", "#009988", "#EE7733",
	"#CC3311", "#EE3377", "#BBBBBB", "#AA3377",
	"#44AA99", "#DDAA33", "#332288", "#88CCEE",
	"#117733", "#882255", "#999933", "#661100",
]

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

BarH = 0.72
SfFs = 6.5
LabelFs = 8.0
LegFs = 7.0
FillAlpha = 0.10
LocusBarAlpha = 0.20
HeaderH = 0.32
AbcPad = 30
LeftExtra = 50
SfHiPad = 20

def ReadFasta(fns):
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
			key = NormKey(label)
			# First FASTA wins on key clash (prefer palmcore over full-length).
			if key not in by_key:
				seqs[label] = seq
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

def DrawTruncArrow(ax, x_edge, y, inward, col, arrow_w):
	y0 = y - BarH / 2.0
	y1 = y + BarH / 2.0
	tip_x = x_edge
	base_x = x_edge + inward * arrow_w
	ax.add_patch(Polygon(
		((tip_x, y), (base_x, y0), (base_x, y1)),
		closed=True, facecolor=to_rgba(col, FillAlpha), edgecolor="none",
		zorder=6))

seqs, fasta_by_key = ReadFasta(Args.fasta)
motif_annot = ReadMotifs(Args.motifs, Args.motif_style)
domains_raw = ReadDomains(Args.domains, Args.cath)
order_ids = ReadOrder(Args.order)
if Args.cath:
	sf_name_map = ReadCathNames(Args.sf_names)
else:
	sf_name_map = ReadScopNames(Args.sf_names)

if Args.motif_style == "abc":
	LocusColors = LocusColorsAbc
	LocusNames = ("A", "B", "C")
else:
	LocusColors = LocusColorsIii
	LocusNames = ("I", "II", "III")

# Index domain queries and fasta labels by order id.
dom_by_oid = {}
for query, doms in domains_raw.items():
	oid = OrderIdOf(query, Args.motif_style)
	if Args.motif_style == "iii":
		oid = oid.upper()
	dom_by_oid.setdefault(oid, []).append(query)

fa_by_oid = {}
for key, fa_label in fasta_by_key.items():
	oid = OrderIdOf(fa_label, Args.motif_style)
	if Args.motif_style == "iii":
		oid = oid.upper()
	fa_by_oid.setdefault(oid, []).append(fa_label)

motifs = {}
domains = {}
tip_order = []
for oid in order_ids:
	fa_cands = fa_by_oid.get(oid, [])
	dom_cands = dom_by_oid.get(oid, [])
	# Prefer a domain-query label whose NormKey exists in fasta.
	chosen = None
	for q in dom_cands:
		if NormKey(q) in fasta_by_key:
			chosen = fasta_by_key[NormKey(q)]
			query_for_dom = q
			break
	else:
		query_for_dom = None
		if len(fa_cands) > 0:
			chosen = fa_cands[0]
		else:
			sys.stderr.write("skip no fasta for order id %s\n" % oid)
			continue

	key = NormKey(chosen)
	seq = seqs[chosen]
	qlen = len(seq)
	# Use domain rows from matching query if present.
	doms = []
	if query_for_dom is not None:
		doms = sorted(domains_raw[query_for_dom],
			key=lambda d: (d[1], d[2], d[0]))
		for sf, lo, hi in doms:
			if hi > qlen:
				Die("domain past end %s %s %d>%d" %
					(query_for_dom, sf, hi, qlen))
	elif oid in dom_by_oid:
		sys.stderr.write("warn domains without matching fasta length for %s\n"
			% oid)

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

	tip_label = chosen
	tip_order.append(tip_label)
	motifs[tip_label] = (spans[0], spans[1], spans[2], qlen)
	domains[tip_label] = doms

if len(tip_order) == 0:
	Die("no tips to draw")
n_leaf = len(tip_order)
tip_y = {}
for i, label in enumerate(tip_order):
	tip_y[label] = float(n_leaf - 1 - i)

sf_order = []
sf_seen = set()
for label in tip_order:
	for sf in SfArchitecture(domains[label]):
		if sf not in sf_seen:
			sf_seen.add(sf)
			sf_order.append(sf)
SfColors = {}
for i, sf in enumerate(sf_order):
	SfColors[sf] = Palette[i % len(Palette)]

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
	for sf, start, end in domains.get(label, []):
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
max_pdb_chars = 0
for label in tip_order:
	pdb = DisplayLabel(label)
	disp_labels[label] = pdb
	if len(pdb) > max_pdb_chars:
		max_pdb_chars = len(pdb)

fig_w = 12.0
fig_h = max(5.4, 0.30 * n_leaf + 2.8)
right_in = 0.18
bottom_in = 1.85
top_in = 0.42
ax_h_frac = 1.0 - (bottom_in + top_in) / fig_h
ax_y_frac = bottom_in / fig_h

fig = plt.figure(figsize=(fig_w, fig_h))
ax_bar = fig.add_axes([0.22, ax_y_frac, 0.72, ax_h_frac])
ax_bar.set_ylim(-BarH / 2.0 - 0.12,
	n_leaf - 1.0 + BarH / 2.0 + HeaderH + 1.05)

def ProbeIn(s, fs=LabelFs, family="monospace"):
	t = fig.text(0, 0, s, fontsize=fs, fontfamily=family)
	fig.canvas.draw()
	w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
	t.remove()
	return w

upper_in = ProbeIn("MM") - ProbeIn("M")

bar_ax_in = fig_w * 0.72
bar_left_in = fig_w * 0.22
plot_w = bar_left_in + bar_ax_in + right_in
longest_leg = ""
leg_labels = {}
for sf in sf_order:
	leg_labels[sf] = "%s %s" % (sf, SfDisplayName(sf, sf_name_map))
	if len(leg_labels[sf]) > len(longest_leg):
		longest_leg = leg_labels[sf]
leg_txt_in = ProbeIn(longest_leg, fs=LegFs, family="sans-serif") if longest_leg else 0.0
needed_leg_in = 2.0 * (0.22 + 0.12 + leg_txt_in) + 0.50
fig_w = max(plot_w, needed_leg_in + 0.30)
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
by_y = sorted(tip_order, key=lambda L: -tip_y[L])
top_label = by_y[0]
mA_top, mB_top, mC_top = motifs[top_label][:3]
cen_top = MotifCenter(mB_top, motifs[top_label][3])
y_top_bar = tip_y[top_label] + BarH / 2.0
for name, mot in zip(LocusNames, (mA_top, mB_top, mC_top)):
	if mot is None:
		continue
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

for i in range(len(by_y) - 1):
	lab_a = by_y[i]
	lab_b = by_y[i + 1]
	y_a = tip_y[lab_a]
	y_b = tip_y[lab_b]
	mA_a, mB_a, mC_a = motifs[lab_a][:3]
	mA_b, mB_b, mC_b = motifs[lab_b][:3]
	cen_a = MotifCenter(mB_a, motifs[lab_a][3])
	cen_b = MotifCenter(mB_b, motifs[lab_b][3])
	y_top = y_a - BarH / 2.0
	y_bot = y_b + BarH / 2.0
	for name, mot_a, mot_b in zip(LocusNames,
			(mA_a, mB_a, mC_a), (mA_b, mB_b, mC_b)):
		if mot_a is None or mot_b is None:
			continue
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

for label in tip_order:
	y = tip_y[label]
	cls = ClassOf(label)
	col = LabelColors[cls]
	mA, mB, mC, qlen = motifs[label]
	center = MotifCenter(mB, qlen)
	x1, x2 = SegX(1, qlen, center)
	left_trunc = x1 < XViewLo
	right_trunc = x2 > XViewHi
	clip = ClipSeg(x1, x2, XViewLo, XViewHi)
	if clip is None:
		continue
	vx1, vx2 = clip
	y0 = y - BarH / 2.0
	y1 = y + BarH / 2.0
	draw1 = vx1 + ArrowW if left_trunc else vx1
	draw2 = vx2 - ArrowW if right_trunc else vx2
	if draw2 <= draw1:
		draw1, draw2 = vx1, vx2
		left_trunc = False
		right_trunc = False
	ax_bar.add_patch(Rectangle((draw1, y0), draw2 - draw1, BarH,
		facecolor=to_rgba(col, FillAlpha), edgecolor="none", zorder=1))
	edge = to_rgba(col, 0.55)
	ax_bar.plot([draw1, draw2], [y0, y0], color=edge, lw=0.6, zorder=4,
		solid_capstyle="butt")
	ax_bar.plot([draw1, draw2], [y1, y1], color=edge, lw=0.6, zorder=4,
		solid_capstyle="butt")
	if not left_trunc:
		ax_bar.plot([vx1, vx1], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
	if not right_trunc:
		ax_bar.plot([vx2, vx2], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
	doms = domains.get(label, [])
	text_ok = ContainedTextMask(doms)
	for (sf, start, end), draw_text in zip(doms, text_ok):
		clip_d = ClipSeg(*SegX(start, end, center), XViewLo, XViewHi)
		if clip_d is None:
			continue
		dvx1, dvx2 = clip_d
		used_sf[sf] = True
		ax_bar.add_patch(Rectangle((dvx1, y - BarH / 2.0), dvx2 - dvx1, BarH,
			facecolor=SfColors[sf], edgecolor="none", zorder=2))
		if draw_text:
			ax_bar.text((dvx1 + dvx2) / 2.0, y, sf, ha="center", va="center",
				fontsize=SfFs, fontweight="bold", color="white", zorder=3,
				clip_on=True)
	if left_trunc:
		cut_marks.append((XViewLo, y, 1.0, col))
	if right_trunc:
		cut_marks.append((XViewHi, y, -1.0, col))

for label in tip_order:
	y = tip_y[label]
	mA, mB, mC, qlen = motifs[label]
	center = MotifCenter(mB, qlen)
	for name, mot in zip(LocusNames, (mA, mB, mC)):
		if mot is None:
			continue
		clip = ClipSeg(*SegX(mot[0], mot[1], center), XViewLo, XViewHi)
		if clip is None:
			continue
		vx1, vx2 = clip
		ax_bar.add_patch(Rectangle((vx1, y - BarH / 2.0), vx2 - vx1, BarH,
			facecolor=to_rgba(LocusColors[name], LocusBarAlpha), edgecolor="none",
			zorder=5))

for x_edge, y, inward, col in cut_marks:
	DrawTruncArrow(ax_bar, x_edge, y, inward, col, ArrowW)

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

handles = []
for sf in sf_order:
	if sf not in used_sf:
		continue
	handles.append(Patch(facecolor=SfColors[sf], edgecolor="none",
		label=leg_labels[sf]))

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
	_, leg_top = inv.transform((0.0, lab_bottom - 4.0))
else:
	leg_top = 0.12
if len(handles) > 0:
	fig.legend(handles=handles, loc="upper center", fontsize=LegFs,
		frameon=False, bbox_to_anchor=(0.5, leg_top), bbox_transform=fig.transFigure,
		ncol=1, handletextpad=0.5, columnspacing=1.4, labelspacing=0.35)

OutPath = Args.output
fig.savefig(OutPath, dpi=200)
root, ext = os.path.splitext(OutPath)
if ext.lower() == ".pdf":
	fig.savefig(root + ".png", dpi=200)
elif ext.lower() == ".png":
	fig.savefig(root + ".pdf", dpi=200)
plt.close(fig)
sys.stderr.write("wrote %s\n" % OutPath)
