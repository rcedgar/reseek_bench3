#!/usr/bin/python3

import os
import sys

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)
DefaultFasta = os.path.join(PalmDir, "reference_data", "ref.fa")
DefaultAnnot = os.path.join(PalmDir, "reference_annotations", "motifseqs.tsv")

from ref_label_map import BaseOf, ReadSuffixMap

def MotifOrNone(s):
	if s is None:
		return None
	t = s.strip().upper()
	if t == "" or t == ".":
		return None
	return t

def ReadFasta(fn):
	seqs = {}
	label = None
	parts = []
	def Flush():
		nonlocal label, parts
		if label is None:
			return
		seq = "".join(parts)
		if label in seqs:
			if len(seq) > len(seqs[label]):
				seqs[label] = seq
		else:
			seqs[label] = seq
		parts = []
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		if Line.startswith(">"):
			Flush()
			label = Line[1:].split()[0]
		else:
			if label is None:
				continue
			parts.append(Line.strip().replace(" ", "").upper())
	Flush()
	return seqs

def ReadMotifAnnot(fn):
	annot = {}
	for Line in open(fn):
		Line = Line.rstrip("\n")
		if len(Line) == 0:
			continue
		Fields = Line.split("\t")
		if Fields[0] in ("Label", "label"):
			continue
		if len(Fields) < 6:
			continue
		label = Fields[0]
		seqA = MotifOrNone(Fields[3])
		seqB = MotifOrNone(Fields[4])
		seqC = MotifOrNone(Fields[5])
		annot[label] = (seqA, seqB, seqC)
	return annot

def FindMotifStart(seq, motif):
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
	return hits[0]

def MotifStarts(seq, seqA, seqB, seqC):
	starts = {}
	for name, motif in (("A", seqA), ("B", seqB), ("C", seqC)):
		if motif is None:
			return None
		pos = FindMotifStart(seq, motif)
		if pos is None:
			return None
		starts[name] = pos
	return starts

def IsAbcOrder(starts):
	if starts is None:
		return False
	return starts["A"] < starts["B"] < starts["C"]

def IsCabOrder(starts):
	if starts is None:
		return False
	return starts["C"] < starts["A"] < starts["B"]

def LookupSeq(label, seqs):
	seq_labels = list(seqs.keys())
	for seq_label in seq_labels:
		seq_label = seq_label.split()[0].upper()
		if seq_label.startswith(label.upper()):
			return seqs[seq_label]
	return None

def ExcludeLabels(fasta_path=DefaultFasta, annot_path=DefaultAnnot):
	"""Labels lacking complete A/B/C motifs in A-B-C sequence order."""
	seqs = ReadFasta(fasta_path)
	annot = ReadMotifAnnot(annot_path)
	suffix_map = ReadSuffixMap(fasta_path)
	excl = set()
	for label in annot:
		seq = LookupSeq(label, seqs)
		if seq is None:
			excl.add(label)
			continue
		seqA, seqB, seqC = annot[label]
		starts = MotifStarts(seq, seqA, seqB, seqC)
		if not IsAbcOrder(starts):
			excl.add(label)
	for seq_label in seqs:
		if seq_label in annot:
			continue
		if seq_label.endswith("-rdrp") or seq_label.endswith("-decoy"):
			excl.add(seq_label)
			continue
		suffixed = suffix_map.get(seq_label)
		if suffixed is not None and suffixed not in annot:
			excl.add(suffixed)
	return excl

# Backward-compatible alias
CabLabels = ExcludeLabels

def LogExcluded(labels):
	if len(labels) == 0:
		return
	parts = sorted(labels)
	sys.stderr.write("excluded_non_abc=%d labels: %s\n" % (len(parts), " ".join(parts)))

if __name__ == "__main__":
	excl = ExcludeLabels()
	LogExcluded(excl)
	for lab in sorted(excl):
		sys.stdout.write("%s\n" % lab)
