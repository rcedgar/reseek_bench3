#!/usr/bin/python3
"""Draw cdn CATH-only annotation figure with overlap-aware sub-levels.

Same layout as make_annots_fig.py (motif trapezoids, tip labels, legend),
but CATH domains only. Overlapping domain rectangles are stacked into
vertical sub-levels; repeated SF ids share one exclusive level (overlapping
rects allowed, label shown once). Other domains pack left-to-right.
SF id labels sit in outlined boxes right of the bar, aligned per sub-level.
"""

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

# Import shared helpers / constants from the multi-db annot figure.
import make_annots_fig as maf

Usage = (
"Draw cdn tips with CATH domains only. Overlapping domains/labels are "
"stacked into sub-levels within each bar (left-to-right, first free level)."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)
CathLegFsScale = 1.25
LegRowInBase = 0.155


def IntervalsOverlap(a1, a2, b1, b2):
	"""True if intervals [a1,a2] and [b1,b2] overlap in x."""
	return not (a2 <= b1 or b2 <= a1)


def AssignCathLevels(doms, center, SegX, text_half_x_fn, XViewLo, XViewHi):
	"""Assign each CATH domain to a sub-level.

	Domains whose SF id appears more than once share one exclusive level
	(all rects drawn there, overlap allowed; label shown once). Other
	domains use greedy left-to-right packing with rect/label overlap checks.
	Levels reserved for a duplicate SF cannot host any other SF.
	"""
	from collections import Counter

	ordered = sorted(doms, key=lambda d: (d[1], d[2], d[0]))
	sf_counts = Counter(sf for sf, _s, _e in ordered)
	duplicate_sfs = {sf for sf, c in sf_counts.items() if c > 1}

	def RectSpan(start, end):
		x1, x2 = SegX(start, end, center)
		rx1 = max(x1, XViewLo)
		rx2 = min(x2, XViewHi)
		if rx2 <= rx1:
			return None
		return rx1, rx2

	def TextSpan(rx1, rx2, sf):
		mid = 0.5 * (rx1 + rx2)
		hw = text_half_x_fn(sf)
		return mid - hw, mid + hw

	def LevelOk(rx1, rx2, tx1, tx2, level):
		if level["exclusive"] is not None:
			for irx1, irx2, _, _ in level["items"]:
				if IntervalsOverlap(rx1, rx2, irx1, irx2):
					return False
				if IntervalsOverlap(tx1, tx2, irx1, irx2):
					return False
			ltx1, ltx2 = level["text"]
			if IntervalsOverlap(tx1, tx2, ltx1, ltx2):
				return False
			if IntervalsOverlap(rx1, rx2, ltx1, ltx2):
				return False
			return True
		for irx1, irx2, itx1, itx2 in level["items"]:
			if IntervalsOverlap(rx1, rx2, irx1, irx2):
				return False
			if IntervalsOverlap(tx1, tx2, itx1, itx2):
				return False
			if IntervalsOverlap(rx1, rx2, itx1, itx2):
				return False
			if IntervalsOverlap(tx1, tx2, irx1, irx2):
				return False
		return True

	levels = []  # {exclusive: sf|None, items: [...], text: (tx1,tx2)}
	placed = []

	sf_first = {}
	for sf, start, _end in ordered:
		if sf not in sf_first:
			sf_first[sf] = start
	dup_order = sorted(duplicate_sfs, key=lambda sf: (sf_first[sf], sf))

	for sf in dup_order:
		spans = []
		dom_rows = []
		for ds, start, end in ordered:
			if ds != sf:
				continue
			span = RectSpan(start, end)
			if span is None:
				continue
			spans.append(span)
			dom_rows.append((sf, start, end))
		if len(spans) == 0:
			continue
		union_lo = min(s[0] for s in spans)
		union_hi = max(s[1] for s in spans)
		tx1, tx2 = TextSpan(union_lo, union_hi, sf)
		items = [(rx1, rx2, tx1, tx2) for rx1, rx2 in spans]
		li = len(levels)
		levels.append({"exclusive": sf, "items": items, "text": (tx1, tx2)})
		for row in dom_rows:
			placed.append(row + (li,))

	for sf, start, end in ordered:
		if sf in duplicate_sfs:
			continue
		span = RectSpan(start, end)
		if span is None:
			continue
		rx1, rx2 = span
		tx1, tx2 = TextSpan(rx1, rx2, sf)
		chosen = None
		for li, level in enumerate(levels):
			if level["exclusive"] is not None:
				continue
			if LevelOk(rx1, rx2, tx1, tx2, level):
				chosen = li
				break
		if chosen is None:
			chosen = len(levels)
			levels.append({"exclusive": None, "items": [], "text": (0.0, 0.0)})
		levels[chosen]["items"].append((rx1, rx2, tx1, tx2))
		placed.append((sf, start, end, chosen))

	if len(placed) == 0:
		return [], 1
	return placed, max(1, len(levels))


def LevelStripWidthX(domains, text_w_x_fn, label_gap_x):
	"""Width in data-x of horizontally packed label boxes for one level."""
	if len(domains) == 0:
		return 0.0
	seen = []
	for sf, _s, _e in domains:
		if sf not in seen:
			seen.append(sf)
	w = label_gap_x
	for i, sf in enumerate(seen):
		w += text_w_x_fn(sf)
		if i + 1 < len(seen):
			w += label_gap_x
	return w


def OutlinePadY(ax, lw_pts):
	"""Data-y padding matching one linewidth (for label box inset)."""
	fig = ax.figure
	fig.canvas.draw()
	bbox = ax.get_window_extent()
	yspan = ax.get_ylim()[1] - ax.get_ylim()[0]
	if bbox.height <= 0.0 or yspan <= 0.0:
		return 0.002
	lw_in = lw_pts / 72.0
	return lw_in * yspan / (bbox.height / fig.dpi)


LABEL_OUTLINE_LW = 0.6


def DrawRightLevelLabels(ax, domains, bar_x2, lane_lo, lane_hi, label_gap_x,
		sf_colors, text_w_x_fn, fs):
	"""Outlined label boxes right of bar; each SF id shown once per level."""
	if len(domains) == 0:
		return bar_x2
	seen = []
	for sf, _start, _end in domains:
		if sf not in seen:
			seen.append(sf)
	pad_y = OutlinePadY(ax, LABEL_OUTLINE_LW)
	lane_h = lane_hi - lane_lo
	box_h = max(lane_h - 2.0 * pad_y, lane_h * 0.5)
	box_lo = lane_lo + 0.5 * (lane_h - box_h)
	cy = box_lo + 0.5 * box_h
	x = bar_x2 + label_gap_x
	for sf in seen:
		box_w = text_w_x_fn(sf)
		ax.add_patch(Rectangle((x, box_lo), box_w, box_h,
			facecolor="none", edgecolor=sf_colors[sf],
			linewidth=LABEL_OUTLINE_LW, clip_on=False, zorder=6))
		ax.text(x + 0.5 * box_w, cy, sf, ha="center", va="center",
			fontsize=fs, color="black", clip_on=False, zorder=7)
		x += box_w + label_gap_x
	return x


# Consecutive entries differ in hue family so sorted CATH ids contrast
# (HSV put 3.30.210.10 and 3.30.460.10 both in blue).
CATH_CONTRAST_COLORS = (
	"#E41A1C",  # red
	"#377EB8",  # blue
	"#4DAF4A",  # green
	"#984EA3",  # purple
	"#FF7F00",  # orange
	"#1B9E77",  # teal
	"#F781BF",  # pink
	"#A65628",  # brown
	"#17BECF",  # cyan
	"#6A3D9A",  # deep violet
	"#E6AB02",  # mustard
	"#66A61E",  # olive
)


def CathContrastColor(i, n):
	"""Qualitative colors for CATH fills and label outlines."""
	pal = CATH_CONTRAST_COLORS
	if n <= len(pal):
		return pal[i]
	import colorsys
	h = (float(i) * 0.618033988749895) % 1.0
	s = 0.82
	v = 0.88 if (i % 2 == 0) else 0.70
	r, g, b = colorsys.hsv_to_rgb(h, s, v)
	return "#%02X%02X%02X" % (int(round(r * 255)), int(round(g * 255)),
		int(round(b * 255)))


def Main():
	AP = argparse.ArgumentParser(description=Usage)
	AP.add_argument("--output", default=None,
		help="Comma-separated output files "
			"(default figs/cdn_cath_annots.pdf,figs/cdn_cath_annots.png)")
	Args = AP.parse_args()
	family = "cdn"
	wanted = ["CATH"]
	paths = maf.FamilyPaths(family)
	for key in ("sifts", "motifs", "order", "fasta", "category"):
		if not os.path.isfile(paths[key]):
			maf.Die("missing %s (run reformat_sifts_annots.py first for sifts)" %
				paths[key])

	if Args.output is None:
		out_dir = os.path.join(PalmDir, "figs")
		out_paths = [
			os.path.join(out_dir, "cdn_cath_annots.pdf"),
			os.path.join(out_dir, "cdn_cath_annots.png"),
		]
	else:
		out_paths = [p.strip() for p in Args.output.split(",") if p.strip()]
		if len(out_paths) == 0:
			maf.Die("empty --output")

	locus_cols = maf.CdnLocusCols
	locus_names = maf.CdnLocusNames
	locus_colors = maf.LocusColorsIii

	seqs, by_key = maf.ReadFasta(paths["fasta"])
	order_ids = maf.ReadOrder(paths["order"])
	cat_by_chain, cat_by_pdb = maf.ReadCategoryMap(paths["category"], family)

	sifts_doms, _, id_meta = maf.ReadAnnotTsv(
		paths["sifts"], wanted, seqs, by_key, locus_cols)
	_, motif_spans, _ = maf.ReadAnnotTsv(
		paths["motifs"], wanted, seqs, by_key, locus_cols)

	tip_order = []
	tip_motifs = {}
	tip_domains = {}
	disp_labels = {}
	for oid in order_ids:
		fasta_lab = maf.ResolveFastaLabel(oid, seqs, by_key,
			allow_pdb_fallback=True)
		if fasta_lab is None:
			sys.stderr.write("skip no fasta for order id %s\n" % oid)
			continue
		qlen = len(seqs[fasta_lab])
		ms = motif_spans.get(fasta_lab, {})
		spans = [ms.get(col) for col in locus_cols]
		tip_order.append(fasta_lab)
		tip_motifs[fasta_lab] = (spans[0], spans[1], spans[2], qlen)
		tip_domains[fasta_lab] = sifts_doms.get(fasta_lab, {})
		cat = maf.TipCategoryLabel(fasta_lab, cat_by_chain, cat_by_pdb)
		disp_labels[fasta_lab] = maf.StripChainKeepCategory(cat)

	if len(tip_order) == 0:
		maf.Die("no tips to draw")

	def AnchorOf(label):
		_mA, mB, _mC, qlen = tip_motifs[label]
		return maf.MotifCenter(mB, qlen)

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

	# Full-length bars: viewport covers every tip end-to-end.
	XViewLo = None
	XViewHi = None
	for label in tip_order:
		mA, _mB, _mC, qlen = tip_motifs[label]
		center = AnchorOf(label)
		x1, x2 = SegX(1, qlen, center)
		if XViewLo is None or x1 < XViewLo:
			XViewLo = x1
		if XViewHi is None or x2 > XViewHi:
			XViewHi = x2
		if mA is not None:
			lo = max(1, mA[0] - maf.AbcPad)
			ax1 = XOf(lo, center) - 0.5
			if ax1 < XViewLo:
				XViewLo = ax1
	XViewLo -= maf.LeftExtra
	bar_only_hi = XViewHi
	XViewHi += maf.SfHiPad
	if XViewHi is None or XViewHi <= XViewLo:
		maf.Die("empty bar x range")

	# Provisional figure for text-width probing / level assignment.
	fig_w = 12.0
	right_in = 0.18
	bar_left_in = fig_w * 0.20
	bar_ax_in = fig_w * 0.74
	fig = plt.figure(figsize=(fig_w, 8.0))

	def ProbeIn(s, fs=maf.LabelFs, family_font="monospace"):
		t = fig.text(0, 0, s, fontsize=fs, fontfamily=family_font)
		fig.canvas.draw()
		w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
		t.remove()
		return w

	x_span = max(1.0, bar_only_hi - XViewLo + maf.SfHiPad)

	def TextHalfX(sf):
		return (0.5 * ProbeIn(sf, fs=maf.SfFs, family_font="sans-serif")
			* x_span / bar_ax_in)

	tip_placed = {}
	tip_n_levels = {}
	for label in tip_order:
		doms = tip_domains.get(label, {}).get("CATH", [])
		center = AnchorOf(label)
		placed, n_lev = AssignCathLevels(
			doms, center, SegX, TextHalfX, XViewLo, XViewHi)
		tip_placed[label] = placed
		tip_n_levels[label] = n_lev

	def TipBarH(label):
		return float(tip_n_levels[label]) * maf.LaneUnit

	tip_y = {}
	y_cursor = 0.0
	for i, label in enumerate(tip_order):
		bh = TipBarH(label)
		if i > 0:
			y_cursor -= maf.BarGap
		tip_top = y_cursor
		tip_y[label] = tip_top - bh / 2.0
		y_cursor = tip_top - bh
	y_bottom = y_cursor

	sf_set = set()
	for label in tip_order:
		for sf, _s, _e, _lv in tip_placed[label]:
			sf_set.add(sf)
	sf_order = sorted(sf_set, key=maf.CathSortKey)

	SfColors = {}
	for i, sf in enumerate(sf_order):
		SfColors[sf] = CathContrastColor(i, max(1, len(sf_order)))

	def LegText(sf):
		_clan, desc = id_meta.get(("CATH", sf), (None, "(missing)"))
		return "%s %s" % (sf, desc)

	# Figure size with known bar heights.
	y_top_extent = maf.HeaderH + 1.05
	y_bot_extent = y_bottom - 0.20
	n_leg_items = len(sf_order)
	n_leg_rows = (n_leg_items + 1) // 2 if n_leg_items > 0 else 0
	leg_row_in = LegRowInBase * CathLegFsScale
	leg_est_in = 0.35 + n_leg_rows * leg_row_in
	bottom_in = max(2.2, leg_est_in + 0.70)
	top_in = 0.50
	data_span = y_top_extent - y_bot_extent
	fig_h = max(8.5, data_span * 0.42 + bottom_in + top_in + 1.0)
	ax_h_frac = 1.0 - (bottom_in + top_in) / fig_h
	ax_y_frac = bottom_in / fig_h

	# Estimate right-side label strip and extend viewport.
	x_span_prov = max(1.0, bar_only_hi - XViewLo + maf.SfHiPad)
	gap_in_prov = maf.BarGap / data_span * (fig_h * ax_h_frac)
	label_gap_x_prov = gap_in_prov * x_span_prov / bar_ax_in

	def TextWidthXProv(sf):
		return (ProbeIn(sf, fs=maf.SfFs, family_font="sans-serif")
			* x_span_prov / bar_ax_in)

	max_right = bar_only_hi
	for label in tip_order:
		_mA, _mB, _mC, qlen = tip_motifs[label]
		center = AnchorOf(label)
		_vx1, vx2 = SegX(1, qlen, center)
		placed = tip_placed[label]
		n_lev = tip_n_levels[label]
		tip_right = vx2
		for li in range(n_lev):
			doms_li = sorted(
				[(sf, s, e) for sf, s, e, lv in placed if lv == li],
				key=lambda d: (d[1], d[2], d[0]))
			tip_right = max(tip_right, vx2 + LevelStripWidthX(
				doms_li, TextWidthXProv, label_gap_x_prov))
		max_right = max(max_right, tip_right)
	XViewHi = max(max_right, bar_only_hi + maf.SfHiPad)

	strip_in_prov = (XViewHi - bar_only_hi) * bar_ax_in / x_span_prov
	right_in = max(0.18, strip_in_prov + gap_in_prov + 0.08)
	fig_w = bar_left_in + bar_ax_in + right_in

	plt.close(fig)
	fig = plt.figure(figsize=(fig_w, fig_h))
	ax_bar = fig.add_axes([0.20, ax_y_frac, 0.74, ax_h_frac])
	ax_bar.set_ylim(y_bot_extent, y_top_extent)

	def ProbeIn2(s, fs=maf.LabelFs, family_font="monospace"):
		t = fig.text(0, 0, s, fontsize=fs, fontfamily=family_font)
		fig.canvas.draw()
		w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
		t.remove()
		return w

	plot_w = bar_left_in + bar_ax_in + right_in
	cath_leg_fs = maf.LegFs * CathLegFsScale
	mid_leg = (len(sf_order) + 1) // 2

	def LegEntryW(sf):
		return 0.14 + ProbeIn2(LegText(sf), fs=cath_leg_fs,
			family_font="sans-serif")

	if len(sf_order) > 0:
		w0 = max(LegEntryW(sf) for sf in sf_order[:mid_leg])
		w1 = max(LegEntryW(sf) for sf in sf_order[mid_leg:]) if mid_leg < len(
			sf_order) else 0.0
		needed_leg_in = 0.28 + w0 + w1 + 0.35
	else:
		needed_leg_in = 0.28
	fig_w = max(plot_w, needed_leg_in + 0.40)
	xoff = 0.5 * (fig_w - plot_w)
	fig.set_size_inches(fig_w, fig_h)
	bar_ax_in = fig_w - (bar_left_in + right_in)
	ax_bar.set_position([(xoff + bar_left_in) / fig_w, ax_y_frac,
		bar_ax_in / fig_w, ax_h_frac])

	# Match label↔bar gap to vertical BarGap in physical inches.
	ax_h_in = fig_h * ax_h_frac
	gap_in = maf.BarGap / data_span * ax_h_in
	x_span = XViewHi - XViewLo
	label_gap_x = gap_in * x_span / bar_ax_in

	def TextWidthX(sf):
		return ProbeIn2(sf, fs=maf.SfFs, family_font="sans-serif") * x_span / bar_ax_in

	for label in tip_order:
		cls = maf.ClassOf(disp_labels[label])
		col = maf.LabelColors[cls]
		center = AnchorOf(label)
		_mA, _mB, _mC, qlen = tip_motifs[label]
		bar_x1, _ = SegX(1, qlen, center)
		ax_bar.text(bar_x1 - label_gap_x, tip_y[label], disp_labels[label],
			va="center", ha="right", fontsize=maf.LabelFs, fontfamily="monospace",
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
				(vx2, y_top_bar + maf.HeaderH), (vx1, y_top_bar + maf.HeaderH)),
			closed=True, facecolor=col, edgecolor="none", zorder=0))
		ax_bar.text((vx1 + vx2) / 2.0, y_top_bar + maf.HeaderH + 0.04, name,
			ha="center", va="bottom", fontsize=7.0, color="black",
			zorder=8, clip_on=False)

	# Motif trapezoids
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

	# Bars + CATH domains on assigned sub-levels
	for label in tip_order:
		y = tip_y[label]
		cls = maf.ClassOf(disp_labels[label])
		col = maf.LabelColors[cls]
		_mA, _mB, _mC, qlen = tip_motifs[label]
		center = AnchorOf(label)
		bh = TipBarH(label)
		lh = maf.LaneUnit
		n_lev = tip_n_levels[label]
		x1, x2 = SegX(1, qlen, center)
		vx1, vx2 = x1, x2
		y0 = y - bh / 2.0
		y1 = y + bh / 2.0
		ax_bar.add_patch(Rectangle((vx1, y0), vx2 - vx1, bh,
			facecolor=to_rgba(col, maf.FillAlpha), edgecolor="none", zorder=1))
		edge = to_rgba(col, 0.55)
		ax_bar.plot([vx1, vx2], [y0, y0], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([vx1, vx2], [y1, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([vx1, vx1], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		ax_bar.plot([vx2, vx2], [y0, y1], color=edge, lw=0.6, zorder=4,
			solid_capstyle="butt")
		for li in range(1, n_lev):
			yy = y1 - li * lh
			ax_bar.plot([vx1, vx2], [yy, yy], color=to_rgba(col, 0.25),
				lw=0.5, zorder=4, solid_capstyle="butt", linestyle=":")
		for sf, start, end, li in tip_placed[label]:
			clip_d = ClipSeg(*SegX(start, end, center), XViewLo, XViewHi)
			if clip_d is None:
				continue
			dvx1, dvx2 = clip_d
			lane_hi = y1 - li * lh
			lane_lo = lane_hi - lh
			used_sf[sf] = True
			ax_bar.add_patch(Rectangle((dvx1, lane_lo), dvx2 - dvx1,
				lane_hi - lane_lo, facecolor=SfColors[sf], edgecolor="none",
				zorder=2))
		for li in range(n_lev):
			doms_li = sorted(
				[(sf, s, e) for sf, s, e, lv in tip_placed[label] if lv == li],
				key=lambda d: (d[1], d[2], d[0]))
			lane_hi = y1 - li * lh
			lane_lo = lane_hi - lh
			DrawRightLevelLabels(ax_bar, doms_li, vx2, lane_lo, lane_hi,
				label_gap_x, SfColors, TextWidthX, maf.SfFs)

	# Locus overlays
	for label in tip_order:
		y = tip_y[label]
		bh = TipBarH(label)
		mA, mB, mC, _qlen = tip_motifs[label]
		center = AnchorOf(label)
		for name, mot in zip(locus_names, (mA, mB, mC)):
			if mot is None:
				continue
			clip = ClipSeg(*SegX(mot[0], mot[1], center), XViewLo, XViewHi)
			if clip is None:
				continue
			vx1, vx2 = clip
			ax_bar.add_patch(Rectangle((vx1, y - bh / 2.0), vx2 - vx1, bh,
				facecolor=to_rgba(locus_colors[name], maf.LocusBarAlpha),
				edgecolor="none", zorder=5))

	ax_bar.set_xlim(XViewLo, XViewHi)
	shown = []
	k = 200
	while XViewLo + float(k) - 1.0 <= XViewHi + 0.01:
		shown.append(k)
		k += 200
	if XViewLo <= 1.0 - 1.0 + 0.01:
		shown = [1] + shown
	ax_bar.set_xticks([XViewLo + float(s) - 1.0 for s in shown])
	ax_bar.set_xticklabels(["%d" % int(s) for s in shown])
	ax_bar.tick_params(axis="y", left=False, labelleft=False)
	ax_bar.spines["left"].set_visible(False)
	ax_bar.spines["right"].set_visible(False)
	ax_bar.spines["top"].set_visible(False)
	ax_bar.tick_params(axis="x", labelsize=8, pad=2)

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
		xlab = fig.text(fx, fy, "amino acids", fontsize=9, ha="center", va="top",
			transform=fig.transFigure)
		fig.canvas.draw()
		renderer = fig.canvas.get_renderer()
		lab_bottom = xlab.get_window_extent(renderer=renderer).y0
		_, leg_top = inv.transform((0.0, lab_bottom - 6.0))
	else:
		leg_top = ax_y_frac - 0.02

	visible_items = []
	for sf in sf_order:
		if sf not in used_sf:
			continue
		visible_items.append((sf, LegText(sf)))

	if len(visible_items) > 0:
		cath_leg_fs = maf.LegFs * CathLegFsScale
		n_vis_items = len(visible_items)
		mid_vis = (n_vis_items + 1) // 2
		leg_cols = [visible_items[:mid_vis], visible_items[mid_vis:]]
		n_rows = max(len(c) for c in leg_cols)
		row_dy_in = 0.12 * CathLegFsScale
		swatch_h_in = 0.07 * CathLegFsScale
		swatch_w = 0.010
		gap = 0.008
		col_gap_in = 0.35
		leg_need_in = 0.15 + n_rows * row_dy_in
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
		swatch_h = swatch_h_in / fig_h
		left_x = 0.06

		def ColWidthIn(items):
			if len(items) == 0:
				return 0.0
			return max(0.14 + ProbeIn2(txt, fs=cath_leg_fs,
				family_font="sans-serif") for _sf, txt in items)

		col1_x = left_x + (ColWidthIn(leg_cols[0]) + col_gap_in) / fig_w
		col_xs = [left_x, col1_x]
		for ci, items in enumerate(leg_cols):
			if len(items) == 0:
				continue
			lx = col_xs[ci]
			y = leg_top
			for sf, txt in items:
				fig.patches.append(Rectangle(
					(lx, y - 0.70 * swatch_h), swatch_w, swatch_h,
					transform=fig.transFigure, facecolor=SfColors[sf],
					edgecolor="none", clip_on=False))
				fig.text(lx + swatch_w + gap, y, txt, fontsize=cath_leg_fs,
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
