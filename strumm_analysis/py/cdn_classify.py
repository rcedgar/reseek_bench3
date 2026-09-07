#!/usr/bin/env python3

import sys
import fasta

a3m = sys.argv[1]

label2row = {}
def onseq(label, seq):
	row = ""
	for c in seq:
		if c.isupper() or c == "-":
			row += c
	label2row[label] = row

fasta.ReadSeqsOnSeq(a3m, onseq)

for label, row in label2row.items():
	motifI_hGS = row[34:37]
	motifII_DxD = row[53:56]
	motifIII_Dvil = row[114:117]

	critical = motifI_hGS[1] + motifI_hGS[2] + motifII_DxD[0] + motifII_DxD[2] + motifIII_Dvil[0]
	critical_slash = motifI_hGS[1] + motifI_hGS[2] + "/" + motifII_DxD[0] + "-" + motifII_DxD[2] + "/" + motifIII_Dvil [0]

	cat = "(other)"
	if critical_slash == "GS/D-D/D":
		cat = "Cdn-like"
	elif critical_slash == "GS/E-D/D":
		cat = "CCa/cGAS"
	elif critical_slash == "MN/E-E/Q":
		cat = "MAB21"
	elif critical_slash == "GG/D-D/R":
		cat = "Pol-mu"
	elif critical_slash == "GS/D-D/R":
		cat = "Pol-beta"

	s = label
	s += "\t" + critical_slash
	s += "\t" + cat
	print(s)
