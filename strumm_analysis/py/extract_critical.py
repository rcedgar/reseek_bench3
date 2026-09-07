#!/usr/bin/python3
"""Extract D-G-D, GDD, and gate residues from rdrp motif annots TSV."""

import argparse
import os
import sys

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)

DefaultInput = os.path.join(PalmDir, "annots", "rdrp_motif_annots.tsv")
DefaultOutput = os.path.join(PalmDir, "annots", "rdrp_critical.tsv")

Needed = ("A_seq", "B_seq", "C_seq")


def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)


def ReadMotifAnnots(fn):
	"""Return {query: {A_seq/B_seq/C_seq: seq}} from long motif annots TSV."""
	by_q = {}
	with open(fn) as f:
		for line in f:
			line = line.rstrip("\n")
			if len(line) == 0:
				continue
			fields = line.split("\t")
			if fields[0] == "query":
				continue
			if len(fields) < 6:
				Die("bad annot row: %s" % line[:80])
			query, _db, dbinfo, _start, _end, seq = fields[:6]
			if dbinfo not in Needed:
				continue
			by_q.setdefault(query, {})[dbinfo] = seq
	return by_q


def ExtractRow(query, motifs):
	"""Return (dgd, gdd, gate) or None if motifs incomplete/too short."""
	a = motifs.get("A_seq")
	b = motifs.get("B_seq")
	c = motifs.get("C_seq")
	if a is None or b is None or c is None:
		sys.stderr.write("warning: %s missing A/B/C motif\n" % query)
		return None
	if len(a) <= 10:
		sys.stderr.write("warning: %s A_seq too short (%d)\n" % (query, len(a)))
		return None
	if len(b) <= 1:
		sys.stderr.write("warning: %s B_seq too short (%d)\n" % (query, len(b)))
		return None
	if len(c) < 7:
		sys.stderr.write("warning: %s C_seq too short (%d)\n" % (query, len(c)))
		return None
	d1 = a[5]
	g = b[1]
	d2 = c[5]
	dgd = "%s-%s-%s" % (d1, g, d2)
	if d1 != "D" or g != "G" or d2 != "D":
		sys.stderr.write("warning: %s D-G-D is %s (expected D-G-D)\n" %
			(query, dgd))
	gdd = c[4:7]
	gate = a[10]
	return dgd, gdd, gate


def Main():
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument("--input", default=DefaultInput,
		help="Long motif annots TSV (default annots/rdrp_motif_annots.tsv)")
	ap.add_argument("--output", default=DefaultOutput,
		help="Output TSV (default annots/rdrp_critical.tsv)")
	args = ap.parse_args()

	if not os.path.isfile(args.input):
		Die("missing input %s" % args.input)

	by_q = ReadMotifAnnots(args.input)
	outdir = os.path.dirname(args.output)
	if outdir and not os.path.isdir(outdir):
		os.makedirs(outdir)

	n = 0
	with open(args.output, "w") as out:
		out.write("query\tD-G-D\tGDD\tgate\n")
		for query in sorted(by_q.keys()):
			row = ExtractRow(query, by_q[query])
			if row is None:
				continue
			dgd, gdd, gate = row
			out.write("%s\t%s\t%s\t%s\n" % (query, dgd, gdd, gate))
			n += 1
	sys.stderr.write("wrote %d rows to %s\n" % (n, args.output))


if __name__ == "__main__":
	Main()
