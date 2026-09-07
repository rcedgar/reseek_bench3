#!/usr/bin/python3

import os
import sys

ScriptDir = os.path.dirname(os.path.abspath(__file__))
RepoDir = os.path.dirname(ScriptDir)

CalFn = os.path.join(RepoDir, "data", "balibase.cal")
MapFn = os.path.join(RepoDir, "map_accs", "struct_original.tsv")
OutFn = os.path.join(RepoDir, "data", "balibase_old_accs.cal")

def Die(msg):
	sys.stderr.write(msg + "\n")
	sys.exit(1)

def ReadMap(fn):
	m = {}
	for Line in open(fn, newline=""):
		Line = Line.rstrip("\n")
		if Line == "":
			continue
		Fields = Line.split("\t")
		if len(Fields) != 2:
			Die("bad map row: %s" % Line[:80])
		new, old = Fields
		if new in m:
			Die("duplicate new acc: %s" % new)
		m[new] = old
	return m

def Main():
	Map = ReadMap(MapFn)
	nfound = 0
	nmissing = 0
	with open(OutFn, "w", newline="") as Out:
		for Line in open(CalFn, newline=""):
			if Line.startswith(">"):
				acc = Line[1:].rstrip("\n")
				if acc in Map:
					Out.write(">" + Map[acc] + "\n")
					nfound += 1
				else:
					Out.write(">missing_" + acc + "\n")
					nmissing += 1
			else:
				Out.write(Line)
	sys.stderr.write("found=%d missing=%d\n" % (nfound, nmissing))

if __name__ == "__main__":
	Main()
