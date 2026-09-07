#!/usr/bin/python3

import argparse
import os
import sys

Usage = \
(
"Write a Jalview sequence-features file from RCE motif annot TSV. "
"seqI/seqII/seqIII are colored blue/green/red. Motif coordinates are "
"taken from --fasta (1-based residue numbers). Annot labels are joined "
"to FASTA IDs by PDB id."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
CdntaseDir = os.path.dirname(ScriptDir)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--annot",
	default=os.path.join(CdntaseDir, "rce_motif_analysis", "annot.tsv"))
AP.add_argument("--fasta",
	default=os.path.join(CdntaseDir, "ref", "ref.fa"))
AP.add_argument("--output",
	default=os.path.join(CdntaseDir, "rce_motif_analysis", "annot.features"))
Args = AP.parse_args()

def Warn(msg):
	sys.stderr.write(msg + "\n")

def OpenLines(fn):
	if not os.path.isfile(fn):
		Warn("missing file %s" % fn)
		return []
	try:
		return open(fn).read().splitlines()
	except IOError as e:
		Warn("cannot read %s: %s" % (fn, e))
		return []

def PdbId(label):
	s = label
	for suf in ("-cdn", "-decoy"):
		if s.endswith(suf):
			s = s[:-len(suf)]
	if "_" in s:
		s = s.split("_")[0]
	elif "-" in s:
		s = s.split("-")[0]
	return s

def ReadFasta(fn):
	seqs = {}
	order = []
	label = None
	parts = []
	def Flush():
		if label is None:
			return
		seq = "".join(parts).replace(" ", "").upper()
		if label in seqs:
			Warn("duplicate fasta label %s, keeping first" % label)
			return
		if len(seq) == 0:
			Warn("empty sequence %s" % label)
			return
		seqs[label] = seq
		order.append(label)
	for Line in OpenLines(fn):
		if len(Line) == 0:
			continue
		if Line.startswith(">"):
			Flush()
			label = Line[1:].split()[0]
			parts = []
		else:
			if label is None:
				Warn("sequence before header in %s" % fn)
				continue
			parts.append(Line.strip())
	Flush()
	if len(seqs) == 0:
		Warn("no sequences in %s" % fn)
	return order, seqs

def ReadAnnot(fn):
	rows = []
	for Line in OpenLines(fn):
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] == "label":
			continue
		if len(Fields) < 4:
			Warn("bad annot row: %s" % Line[:80])
			continue
		src = Fields[0]
		seqI = Fields[1].upper()
		seqII = Fields[2].upper()
		seqIII = Fields[3].upper()
		rows.append((src, seqI, seqII, seqIII))
	if len(rows) == 0:
		Warn("no annot rows in %s" % fn)
	return rows

def FindMotif(seq, motif, seq_id, name):
	if motif is None or motif == "":
		Warn("empty motif %s in %s" % (name, seq_id))
		return None
	hits = []
	start = 0
	while True:
		i = seq.find(motif, start)
		if i < 0:
			break
		hits.append(i)
		start = i + 1
	if len(hits) == 0:
		Warn("motif %s not found in %s" % (name, seq_id))
		return None
	if len(hits) > 1:
		Warn("motif %s not unique in %s (%d hits, using first)" %
			(name, seq_id, len(hits)))
	lo = hits[0] + 1
	hi = lo + len(motif) - 1
	return lo, hi

fasta_order, seqs = ReadFasta(Args.fasta)
fasta_by_pdb = {}
for label in fasta_order:
	pdb = PdbId(label)
	if pdb in fasta_by_pdb and fasta_by_pdb[pdb] != label:
		Warn("duplicate PDB id %s in fasta (%s, %s)" %
			(pdb, fasta_by_pdb[pdb], label))
		continue
	fasta_by_pdb[pdb] = label

annot_rows = ReadAnnot(Args.annot)
matched_pdb = {}
features = []
Loci = (("seqI", "blue"), ("seqII", "green"), ("seqIII", "red"))

for src, seqI, seqII, seqIII in annot_rows:
	pdb = PdbId(src)
	if pdb not in fasta_by_pdb:
		Warn("unmatched annot row: %s" % src)
		continue
	seq_id = fasta_by_pdb[pdb]
	matched_pdb[pdb] = True
	seq = seqs[seq_id]
	motifs = {"seqI": seqI, "seqII": seqII, "seqIII": seqIII}
	for name, _col in Loci:
		span = FindMotif(seq, motifs[name], seq_id, name)
		if span is None:
			continue
		features.append((name, seq_id, span[0], span[1]))

for label in fasta_order:
	pdb = PdbId(label)
	if pdb not in matched_pdb:
		Warn("no annot for %s" % label)

out_dir = os.path.dirname(os.path.abspath(Args.output))
if out_dir != "" and not os.path.isdir(out_dir):
	try:
		os.makedirs(out_dir)
	except OSError as e:
		Warn("cannot create output dir %s: %s" % (out_dir, e))
		sys.exit(1)
try:
	fOut = open(Args.output, "w")
except IOError as e:
	Warn("cannot write %s: %s" % (Args.output, e))
	sys.exit(1)

for name, col in Loci:
	fOut.write("%s\t%s\n" % (name, col))
for name, seq_id, lo, hi in features:
	fOut.write("%s\t%s\t-1\t%d\t%d\t%s\n" % (name, seq_id, lo, hi, name))
fOut.close()
