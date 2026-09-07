#!/usr/bin/python3

import argparse
import os
import sys

Usage = \
(
"Convert a motifseqs TSV (Label F1_seq F2_seq A_seq B_seq C_seq D_seq E_seq) "
	"to a Jalview sequence-features file. Motif peptides are located in "
	"the ungapped FASTA sequence (exactly one hit each). Coordinates are "
	"1-based residue positions so Jalview can map them onto an alignment."
)

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)
DefaultIn = os.path.join(PalmDir, "reference_annotations",
	"ref_full_motifseqs_formatted.tsv")
DefaultFasta = os.path.join(PalmDir, "reference_data", "ref_full_labeled.fa")
DefaultOut = os.path.join(PalmDir, "reference_annotations",
	"ref_full_motifs.features")

MotifCols = (
	("F1", 1),
	("F2", 2),
	("A", 3),
	("B", 4),
	("C", 5),
	("D", 6),
	("E", 7),
)

FeatureColors = (
	("A", "0000ff"),
	("B", "00ff00"),
	("C", "ff0000"),
	("D", "ffa500"),
	("E", "ffff00"),
	("F1", "00ffff"),
	("F2", "00ffff"),
)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("--input", default=DefaultIn,
	help="Motif TSV (ref_full_motifseqs_formatted.tsv)")
AP.add_argument("--fasta", default=DefaultFasta,
	help="FASTA with the same sequence labels")
AP.add_argument("--output", default=DefaultOut,
	help="Jalview features file")
Args = AP.parse_args()

from ref_cab_exclude import MotifOrNone, ReadFasta

def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)

def Warning(msg):
	sys.stderr.write("*Warning* " + msg + "\n")

def NormLabel(s):
	s = s.strip()
	if len(s) == 0:
		Die("empty label")
	cut = len(s)
	for i, c in enumerate(s):
		if c in "-/ \t":
			cut = i
			break
	return s[:cut].lower()

def IndexFasta(seqs):
	by_norm = {}
	for raw in seqs:
		key = NormLabel(raw)
		if key in by_norm:
			Die("duplicate FASTA label %s and %s" % (by_norm[key], raw))
		by_norm[key] = raw
	return by_norm

def LookupSeq(label, seqs, by_norm):
	key = NormLabel(label)
	raw = by_norm.get(key)
	if raw is None:
		return None, None
	return raw, seqs[raw]

def FindUniqueStart(seq, motif):
	hits = []
	start = 0
	while True:
		i = seq.find(motif, start)
		if i < 0:
			break
		hits.append(i)
		start = i + 1
	if len(hits) == 0:
		return None, 0
	if len(hits) > 1:
		return None, len(hits)
	return hits[0], 1

seqs = ReadFasta(Args.fasta)
by_norm = IndexFasta(seqs)
rows = []
seen = set()
n_skip = 0

for Line in open(Args.input):
	Line = Line.rstrip("\n")
	if len(Line) == 0:
		continue
	Fields = Line.split("\t")
	if Fields[0] in ("Label", "label"):
		continue
	if len(Fields) < 8:
		Die("bad motif row: %s" % Line[:80])
	label = Fields[0]
	key = NormLabel(label)
	if key in seen:
		Die("duplicate annot label: %s" % label)
	seen.add(key)
	fasta_id, seq = LookupSeq(label, seqs, by_norm)
	if seq is None:
		Warning("no FASTA sequence for %s" % label)
		continue
	seq = seq.replace("-", "")
	for name, col in MotifCols:
		motif = MotifOrNone(Fields[col])
		if motif is None:
			n_skip += 1
			continue
		idx, n_hits = FindUniqueStart(seq, motif)
		if n_hits == 0:
			Warning("motif %s not found in %s: %s" % (name, label, motif))
			continue
		if n_hits > 1:
			Die("motif %s found %d times in %s: %s" %
				(name, n_hits, label, motif))
		lo = idx + 1
		hi = idx + len(motif)
		rows.append((motif, fasta_id, lo, hi, name))

if len(rows) == 0:
	Die("no motif features")

Out = open(Args.output, "w")
for name, color in FeatureColors:
	Out.write("%s\t%s\n" % (name, color))
for motif, label, lo, hi, name in rows:
	Out.write("%s\t%s\t-1\t%d\t%d\t%s\n" %
		(motif, label, lo, hi, name))
Out.close()

sys.stderr.write("labels=%d features=%d skipped_missing=%d\n%s\n" %
	(len(seen), len(rows), n_skip, Args.output))
