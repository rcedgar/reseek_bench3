#!/usr/bin/python3

import argparse
import sys

Usage = \
(
"Convert reseek SCOP/CATH hits to consensus superfamily domain "
"annotations: extend local hits to estimated envelopes, then resolve "
"overlaps with a greedy Pfam-like policy."
)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--input", required=True,
	help="Hits TSV (query target qlo qhi tlo thi qcovpct tcovpct pctid pvalue)")
AP.add_argument("--output", default="-", help="Output TSV (default stdout)")
AP.add_argument("--maxp", type=float, default=None,
	help="Discard hits with pvalue above this (default: no filter)")
AP.add_argument("--max-overlap", type=int, default=8,
	help="Max allowed envelope overlap in aa (default 8)")
Args = AP.parse_args()

def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)

def SuperfamilyOf(target):
	if "/" not in target:
		Die("target lacks '/': %s" % target)
	cls = target.split("/", 1)[1]
	fields = cls.split(".")
	if len(fields) < 4:
		Die("classification has <4 fields: %s" % target)
	return ".".join(fields[:3])

def InferLength(lo, hi, covpct, label):
	if covpct <= 0:
		Die("non-positive coverage for %s" % label)
	aln_len = hi - lo + 1
	if aln_len < 1:
		Die("bad alignment span for %s: %d-%d" % (label, lo, hi))
	return int(round(aln_len * 100.0 / covpct))

def EnvelopeOverlap(a_lo, a_hi, b_lo, b_hi):
	return max(0, min(a_hi, b_hi) - max(a_lo, b_lo) + 1)

def ReadHits(fn, maxp):
	hits_by_query = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0 or Line.startswith("#"):
			continue
		Fields = Line.split("\t")
		if len(Fields) < 10:
			Die("bad hits row (<10 fields): %s" % Line[:80])
		if Fields[0] == "query":
			continue
		query = Fields[0]
		target = Fields[1]
		try:
			qlo = int(Fields[2])
			qhi = int(Fields[3])
			tlo = int(Fields[4])
			thi = int(Fields[5])
			qcovpct = float(Fields[6])
			tcovpct = float(Fields[7])
			pvalue = float(Fields[9])
		except ValueError:
			Die("bad numeric fields: %s" % Line[:80])
		if maxp is not None and pvalue > maxp:
			continue

		sf = SuperfamilyOf(target)
		qlen = InferLength(qlo, qhi, qcovpct, query)
		tlen = InferLength(tlo, thi, tcovpct, target)
		nterm = tlo - 1
		cterm = tlen - thi
		if nterm < 0 or cterm < 0:
			Die("inferred flanks negative for %s vs %s" % (query, target))
		env_lo = max(1, qlo - nterm)
		env_hi = min(qlen, qhi + cterm)
		if env_hi < env_lo:
			Die("bad envelope %d-%d for %s vs %s" %
				(env_lo, env_hi, query, target))

		rec = (pvalue, env_lo, env_hi, qlo, qhi, sf, target)
		hits_by_query.setdefault(query, []).append(rec)
	return hits_by_query

def ResolveOverlaps(hits, max_overlap):
	hits = sorted(hits, key=lambda h: (h[0], h[1], h[2], h[6]))
	accepted = []
	for hit in hits:
		_, env_lo, env_hi = hit[0], hit[1], hit[2]
		ok = True
		for acc in accepted:
			ov = EnvelopeOverlap(env_lo, env_hi, acc[1], acc[2])
			if ov > max_overlap:
				ok = False
				break
		if ok:
			accepted.append(hit)
	accepted.sort(key=lambda h: (h[1], h[2], h[6]))
	return accepted

def main():
	hits_by_query = ReadHits(Args.input, Args.maxp)
	if Args.output == "-":
		Out = sys.stdout
	else:
		Out = open(Args.output, "w")

	Out.write("query\tsf\tenv_lo\tenv_hi\taln_lo\taln_hi\tpvalue\ttarget\n")
	for query in sorted(hits_by_query.keys()):
		for pvalue, env_lo, env_hi, aln_lo, aln_hi, sf, target in \
				ResolveOverlaps(hits_by_query[query], Args.max_overlap):
			Out.write("%s\t%s\t%d\t%d\t%d\t%d\t%g\t%s\n" %
				(query, sf, env_lo, env_hi, aln_lo, aln_hi, pvalue, target))
	if Args.output != "-":
		Out.close()

if __name__ == "__main__":
	main()
