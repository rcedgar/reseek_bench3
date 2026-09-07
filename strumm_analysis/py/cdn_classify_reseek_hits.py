#!/usr/bin/python3

import argparse
import math
import os
import sys

Usage = \
(
"Transfer motif I/II/III from an annotated target onto each query "
"using the local alignment in a hits TSV."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
CdntaseDir = os.path.dirname(ScriptDir)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--input",
	help="Hits TSV (query+target+pctid+qcovpct+pvalue+qlo+tlo+qrow+trow)")
AP.add_argument("--annot",
	default=os.path.join(CdntaseDir, "rce_motif_analysis", "annot.tsv"),
	help="Reference motif TSV from make_annot_tsv.py")
AP.add_argument("--output", default="-", help="Output TSV (default stdout)")
AP.add_argument("--loo-report", default=None, help="Leave-one-out recovery report file")
AP.add_argument("--max-pctid", type=float, default=None,
	help="Discard hits with pctid above this value (all runs)")
AP.add_argument("--maxp", type=float, default=1e-6,
	help="Pvalue cutoff for cdn classification (default 1e-6)")
AP.add_argument("--maxlog10diff", type=float, default=6,
	help="Min -log10(p_cdn)+log10(p_decoy) for cdn (default 6)")
AP.add_argument("--fields", default="2,1",
	help="1-based query,reference field numbers (default 2,1)")
AP.add_argument("--desc", action="store_true",
	help="Split query label at first space; append description column")
Args = AP.parse_args()

def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)

def Warning(msg):
	sys.stderr.write("*Warning* " + msg + "\n")

def ParseFields(s):
	parts = s.split(",")
	if len(parts) != 2:
		Die("bad --fields %s (want i,j)" % s)
	try:
		i = int(parts[0])
		j = int(parts[1])
	except ValueError:
		Die("bad --fields %s (want i,j)" % s)
	if i not in (1, 2) or j not in (1, 2) or i == j:
		Die("bad --fields %s (want 1,2 or 2,1)" % s)
	return i, j

def FirstToken(s):
	toks = s.split()
	if len(toks) == 0:
		Die("empty label field")
	return toks[0]

def SplitAccDesc(raw):
	toks = raw.split(None, 1)
	if len(toks) == 0:
		Die("empty query field")
	acc = toks[0]
	if len(toks) > 1 and toks[1] != "":
		return acc, toks[1]
	return acc, "(no description)"

QField, RField = ParseFields(Args.fields)

def ReadAnnot(fn):
	order = []
	annot = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] == "label":
			continue
		if len(Fields) < 5:
			Die("bad annot row: %s" % Line[:80])
		label, seqI, seqII, seqIII, catalytic = Fields[:5]
		if label in annot:
			Die("duplicate annot label: %s" % label)
		annot[label] = (seqI, seqII, seqIII)
		order.append(label)
	if len(annot) == 0:
		Die("no annot rows in %s" % fn)
	return order, annot

def AlignmentMaps(qrow, trow):
	if len(qrow) != len(trow):
		Die("qrow/trow length mismatch %d/%d" % (len(qrow), len(trow)))
	pairs = []
	t_to_q = []
	q_ungapped = []
	qi = 0
	for qc, tc in zip(qrow, trow):
		if tc == "-":
			if qc != "-":
				q_ungapped.append(qc)
				qi += 1
			continue
		if qc == "-":
			pairs.append((tc, "-"))
			t_to_q.append(None)
		else:
			pairs.append((tc, qc))
			t_to_q.append(qi)
			q_ungapped.append(qc)
			qi += 1
	return pairs, t_to_q, "".join(q_ungapped)

def FindMotif(pairs, motif):
	ungapped = "".join([t for t, q in pairs])
	idx = ungapped.find(motif)
	if idx < 0:
		return None
	return idx

def Extract(pairs, start, L):
	if start < 0 or start + L > len(pairs):
		return None
	s = []
	for i in range(L):
		s.append(pairs[start + i][1])
	return "".join(s)

def IdentI(qseq, annot):
	return qseq[6] == annot[6] and qseq[7] == annot[7]

def IdentII(qseq, annot):
	def de(a, b):
		if a == b:
			return True
		if a in "DE" and b in "DE":
			return True
		return False
	return de(qseq[2], annot[2]) and de(qseq[4], annot[4])

def IdentIII(qseq, annot):
	a = qseq[2]
	b = annot[2]
	if a == b:
		return True
	if a == "D" or b == "D":
		return True
	return False

def QueryIndexAt(t_to_q, anchor_t):
	if anchor_t < 0 or anchor_t >= len(t_to_q):
		return None
	if t_to_q[anchor_t] is not None:
		return t_to_q[anchor_t]
	prev = None
	for t in range(anchor_t):
		if t_to_q[t] is not None:
			prev = t_to_q[t]
	if prev is not None:
		return prev + 1
	for t in range(anchor_t + 1, len(t_to_q)):
		if t_to_q[t] is not None:
			return t_to_q[t]
	return None

def UngappedExtract(t_to_q, q_ungapped, anchor_t, prefix_len, L):
	if anchor_t < 0 or anchor_t >= len(t_to_q):
		return None
	q_idx = QueryIndexAt(t_to_q, anchor_t)
	if q_idx is None:
		return None
	start = q_idx - prefix_len
	if start < 0 or start + L > len(q_ungapped):
		return None
	return q_ungapped[start:start + L]

def Infer(pairs, t_to_q, q_ungapped, motif, ident_fn, cat_pos):
	start = FindMotif(pairs, motif)
	if start is None:
		return None
	L = len(motif)
	chosen = start
	q0 = Extract(pairs, start, L)
	ident0 = (q0 is not None) and ident_fn(q0, motif)
	if q0 is None or not ident0:
		found = False
		for shift in (-1, 1):
			qs = Extract(pairs, start + shift, L)
			if qs is not None and ident_fn(qs, motif):
				chosen = start + shift
				found = True
				break
		if not found and q0 is None:
			return None
	result = UngappedExtract(t_to_q, q_ungapped, chosen + cat_pos, cat_pos, L)
	q_idx = QueryIndexAt(t_to_q, chosen + cat_pos)
	if q_idx is not None:
		base = q_idx - cat_pos
		cands = []
		for shift in (0, -1, 1):
			st = base + shift
			if st < 0 or st + L > len(q_ungapped):
				continue
			cands.append((shift, q_ungapped[st:st + L]))
		picked = None
		for shift, seq in cands:
			if ident_fn(seq, motif):
				picked = (shift, seq)
				break
		if picked is None and len(cands) > 0:
			def nident(seq):
				n = 0
				for a, b in zip(seq, motif):
					if a == b:
						n += 1
				return n
			picked = cands[0]
			bestn = nident(picked[1])
			for shift, seq in cands[1:]:
				n = nident(seq)
				if n > bestn:
					bestn = n
					picked = (shift, seq)
		if picked is not None:
			result = picked[1]
	return result

def Catalytic(sI, sII, sIII):
	return "%s%s/%s-%s/%s" % (sI[6], sI[7], sII[2], sII[4], sIII[2])

CDN_CATS = ("GS/D-D/D",)

def ClassOf(label):
	if label.endswith("-cdn"):
		return "cdn"
	if label.endswith("-decoy"):
		return "decoy"
	Die("unknown suffix: %s" % label)

def NegLog10(p):
	if p == float("inf"):
		return float("-inf")
	if p <= 0:
		return float("inf")
	return -math.log10(p)

def Classify(catalytic, top, min_cdn, min_decoy, maxp, maxlog10diff, motif_ok):
	topclass = ClassOf(top)
	cdn_sig = (topclass == "cdn" and min_cdn <= maxp)
	log10diff = NegLog10(min_cdn) - NegLog10(min_decoy)
	log10_ok = log10diff >= maxlog10diff
	if not motif_ok:
		if cdn_sig and not log10_ok:
			return "notcdn-maxlog10diff"
		if topclass == "decoy":
			return "notcdn-decoy"
		return "notcdn-motif"
	pvalue_test = cdn_sig and log10_ok
	if pvalue_test and catalytic in CDN_CATS:
		return "cdn"
	if catalytic not in CDN_CATS:
		return "notcdn-catalytic"
	if cdn_sig and not log10_ok:
		return "notcdn-maxlog10diff"
	if topclass == "decoy":
		return "notcdn-decoy"
	return "notcdn-other"

def FailWhy(category, catalytic, min_cdn, min_decoy, maxp, maxlog10diff):
	if category == "notcdn-catalytic":
		return "catalytic %s not GS/D-D/D" % catalytic
	if category == "notcdn-maxlog10diff":
		d = NegLog10(min_cdn) - NegLog10(min_decoy)
		return "log10diff %.2f < %g" % (d, maxlog10diff)
	if category == "notcdn-decoy":
		return "top hit is -decoy"
	if category == "notcdn-other":
		if min_cdn > maxp:
			return "cdn pvalue > maxp"
		return "pvalue test failed"
	if category == "notcdn-motif":
		return "motif not found in target"
	if category == "cdn":
		return "pvalue and catalytic tests passed"
	Die("unexpected category: %s" % category)

def MarkCase(shown, other):
	out = []
	for i in range(len(shown)):
		c = shown[i]
		if c in "/-":
			out.append(c)
			continue
		if i >= len(other):
			out.append(c.upper())
			continue
		o = other[i]
		if o in "/-":
			out.append(c.upper())
		elif c.upper() == o.upper():
			out.append(c.lower())
		else:
			out.append(c.upper())
	return "".join(out)

def SameMotif(aI, aII, aIII, aCat, bI, bII, bIII, bCat):
	return (aI.upper() == bI.upper() and aII.upper() == bII.upper() \
		and aIII.upper() == bIII.upper() and aCat.upper() == bCat.upper())

annot_order, annot = ReadAnnot(Args.annot)
LOO = Args.loo_report is not None

# query -> [(pvalue, target, pvalue_s, pctid_s, classify_row, ref_row), ...]
by_query = {}
query_desc = {}
seen = set()
aln_rows = {1: 7, 2: 8}
for Line in open(Args.input):
	Line = Line.rstrip("\n")
	if len(Line) == 0:
		continue
	Fields = Line.split("\t")
	if len(Fields) < 9:
		Die("unexpected columns (%d): %s" % (len(Fields), Line[:80]))
	query_raw = Fields[QField - 1]
	ref_raw = Fields[RField - 1]
	query = FirstToken(query_raw)
	target = FirstToken(ref_raw)
	if Args.desc and query not in query_desc:
		_, desc = SplitAccDesc(query_raw)
		query_desc[query] = desc
	pctid_s = Fields[2]
	pctid = float(pctid_s)
	pvalue_s = Fields[4]
	classify_row = Fields[aln_rows[QField]]
	ref_row = Fields[aln_rows[RField]]
	if target not in annot:
		Die("target not in annot: %s" % target)
	if LOO and query not in annot:
		Die("query not in annot: %s" % query)
	seen.add(query)
	if query == target:
		continue
	if Args.max_pctid is not None and pctid > Args.max_pctid:
		continue
	by_query.setdefault(query, []).append(
		(float(pvalue_s), target, pvalue_s, pctid_s, classify_row, ref_row))

for query in sorted(seen):
	if query not in by_query:
		Warning("no non-self hit: %s" % query)

if Args.output == "-":
	fOut = sys.stdout
else:
	fOut = open(Args.output, "w")

rows = []
for query, neighbors in by_query.items():
	best_p = min([p for p, _, _, _, _, _ in neighbors])
	best = [(t, ps, pct, qr, tr) for p, t, ps, pct, qr, tr in neighbors if p == best_p]
	best.sort(key=lambda x: x[0])
	if len(best) > 1:
		tied = " ".join([t for t, _, _, _, _ in best])
		sys.stderr.write("warning: tie for %s: %s\n" % (query, tied))
	target, pvalue_s, pctid_s, classify_row, ref_row = best[0]
	rows.append((best_p, query, target, pvalue_s, pctid_s, classify_row, ref_row, neighbors))

rows.sort(key=lambda r: (r[0], r[1]))

n = 0
lifted = {}
for _, query, target, pvalue_s, pctid_s, classify_row, ref_row, neighbors in rows:
	seqI, seqII, seqIII = annot[target]
	pairs, t_to_q, q_ungapped = AlignmentMaps(classify_row, ref_row)
	gotI = Infer(pairs, t_to_q, q_ungapped, seqI, IdentI, 6)
	gotII = Infer(pairs, t_to_q, q_ungapped, seqII, IdentII, 2)
	gotIII = Infer(pairs, t_to_q, q_ungapped, seqIII, IdentIII, 2)
	motif_ok = (gotI is not None and gotII is not None and gotIII is not None)
	infI = gotI if gotI is not None else ("-" * len(seqI))
	infII = gotII if gotII is not None else ("-" * len(seqII))
	infIII = gotIII if gotIII is not None else ("-" * len(seqIII))
	infCat = Catalytic(infI, infII, infIII)
	topclass = ClassOf(target)
	min_cdn = float("inf")
	min_decoy = float("inf")
	opp = []
	for p, t, ps, pct, qr, tr in neighbors:
		c = ClassOf(t)
		if c == "cdn" and p < min_cdn:
			min_cdn = p
		if c == "decoy" and p < min_decoy:
			min_decoy = p
		if c != topclass:
			opp.append((p, t, ps))
	if len(opp) == 0:
		opp_t = "NA"
		opp_p = "NA"
	else:
		opp.sort(key=lambda x: (x[0], x[1]))
		opp_t = opp[0][1]
		opp_p = opp[0][2]
	category = Classify(infCat, target, min_cdn, min_decoy, Args.maxp, Args.maxlog10diff, motif_ok)
	out_fields = [query, category, target, pvalue_s, opp_t, opp_p,
		infI, infII, infIII, infCat]
	if Args.desc:
		out_fields.append(query_desc.get(query, "(no description)"))
	fOut.write("\t".join(out_fields) + "\n")
	lifted[query] = (target, pvalue_s, pctid_s, infI, infII, infIII, infCat,
		category, min_cdn, min_decoy)
	n += 1

if Args.output != "-":
	fOut.close()
sys.stderr.write("%d queries\n" % n)

if LOO:
	fRep = open(Args.loo_report, "w")
	n_ok = 0
	n_bad = 0
	n_tp = 0
	n_fn = 0
	n_fp = 0
	first = True
	for label in annot_order:
		if not first:
			fRep.write("\n")
		if label not in lifted:
			Warning("no non-self hit: %s" % label)
			fRep.write("No non-self hit: %s" % label + "\n")
			continue
		first = False
		trueI, trueII, trueIII = annot[label]
		trueCat = Catalytic(trueI, trueII, trueIII)
		target, pvalue_s, pctid_s, infI, infII, infIII, infCat, \
			category, min_cdn, min_decoy = lifted[label]
		if SameMotif(trueI, trueII, trueIII, trueCat, infI, infII, infIII, infCat):
			fRep.write("\t".join([label, target, pvalue_s, pctid_s + "%", infI, infII, infIII, infCat, " correct"]) + "\n")
			n_ok += 1
		else:
			cI = MarkCase(trueI, infI)
			cII = MarkCase(trueII, infII)
			cIII = MarkCase(trueIII, infIII)
			cCat = MarkCase(trueCat, infCat)
			wI = MarkCase(infI, trueI)
			wII = MarkCase(infII, trueII)
			wIII = MarkCase(infIII, trueIII)
			wCat = MarkCase(infCat, trueCat)
			fRep.write("\t".join([label, target, pvalue_s, pctid_s + "%", cI, cII, cIII, cCat, " *ERROR*"]) + "\n")
			fRep.write("\t".join([label, target, pvalue_s, pctid_s + "%", wI, wII, wIII, wCat]) + "\n")
			n_bad += 1
		gold = ClassOf(label)
		called_cdn = (category == "cdn")
		if gold == "cdn" or called_cdn:
			if gold == "cdn" and called_cdn:
				tf = "TP"
				n_tp += 1
			elif gold == "cdn":
				tf = "FN"
				n_fn += 1
			else:
				tf = "FP"
				n_fp += 1
			fields = [label, tf]
			if tf != "TP":
				why = FailWhy(category, infCat, min_cdn, min_decoy,
					Args.maxp, Args.maxlog10diff)
				fields.extend([category, why])
			fRep.write("\t".join(fields) + "\n")
	summary = "loo-report %s OK=%d BAD=%d TP=%d FN=%d FP=%d\n" % (
		Args.loo_report, n_ok, n_bad, n_tp, n_fn, n_fp)
	fRep.write("\n# " + summary + "\n")
	fRep.close()
	sys.stderr.write(summary)
