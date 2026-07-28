#!/usr/bin/python3

import re
import sys

regex = ".EPQ.* ...3=[0-9.eE]*"

def logfile_get_acc(fn):
	truth = "MISSING_TRUTH"
	f = open(fn)
	secs = None
	mem_gb = None
	for line in f:
		if line.find("-truth fam") >= 0:
			truth = "fam"
		elif line.find("-truth superfamilyx") >= 0:
			truth = "sfx"
		elif line.find("-truth superfamily") >= 0 or line.find("-truth sf") >= 0:
			truth = "sf"
		elif line.find("-truth fold") >= 0:
			truth = "fold"
		M = re.search(regex, line)
		if not M is None:
			return M.group(), truth
	return None, None

def main():
	for fn in sys.argv[1:]:
		acc = logfile_get_acc(fn)
		print(acc, fn)
	return 0

if __name__ == "__main__":
	sys.exit(main())
