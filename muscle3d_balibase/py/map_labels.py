#!/usr/bin/python3

import glob
import os
import sys

ScriptDir = os.path.dirname(os.path.abspath(__file__))
RepoDir = os.path.dirname(ScriptDir)

MapFn = os.path.join(RepoDir, "map_labels", "original_struct.tsv")
InDir = os.path.join(RepoDir, "original_labels")
OutDir = os.path.join(RepoDir, "subset_labels")

def Die(msg):
	sys.stderr.write(msg + "\n")
	sys.exit(1)

def ReadOriginalSet(fn):
	s = set()
	structs = set()
	for Line in open(fn, newline=""):
		Line = Line.rstrip("\n")
		if Line == "":
			continue
		Fields = Line.split("\t")
		if len(Fields) != 2:
			Die("bad map row: %s" % Line[:80])
		original, struct = Fields
		if original in s:
			Die("duplicate original: %s" % original)
		if struct in structs:
			Die("duplicate struct: %s" % struct)
		s.add(original)
		structs.add(struct)
	return s

def Main():
	Originals = ReadOriginalSet(MapFn)
	if not os.path.isdir(OutDir):
		os.mkdir(OutDir)

	nwritten = 0
	nskipped = 0
	nfound = 0
	nmissing = 0
	Fns = sorted(glob.glob(os.path.join(InDir, "BB*")))
	for Fn in Fns:
		if not os.path.isfile(Fn):
			continue
		Found = []
		for Line in open(Fn, newline=""):
			Label = Line.rstrip("\n")
			if Label == "":
				continue
			if Label in Originals:
				Found.append(Label)
				nfound += 1
			else:
				nmissing += 1
		if len(Found) > 2:
			OutFn = os.path.join(OutDir, os.path.basename(Fn))
			with open(OutFn, "w", newline="") as Out:
				for Label in Found:
					Out.write(Label + "\n")
			nwritten += 1
		else:
			nskipped += 1
	sys.stderr.write("written=%d skipped=%d found=%d missing=%d\n" %
		(nwritten, nskipped, nfound, nmissing))

if __name__ == "__main__":
	Main()
