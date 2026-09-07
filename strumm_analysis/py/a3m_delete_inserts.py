#!/usr/bin/python3

import sys
import fasta

FN = sys.argv[1]
OneLine = False
if len(sys.argv) > 2:
	if sys.argv[2] == "-oneline":
		OneLine = True
	else:
		assert False

RowLength = None

def OnSeq(Label, Seq):
	global RowLength

	Row = ""
	for c in Seq:
		if c.isupper() or c == "-":
			Row += c
	if RowLength == None:
		RowLength = len(Row)
	else:
		assert len(Row) == RowLength
	if OneLine:
		print(">" + Label)
		print(Row)
	else:
		fasta.WriteSeq(sys.stdout, Row, Label)

fasta.ReadSeqsOnSeq(FN, OnSeq)
