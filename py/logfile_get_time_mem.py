#!/usr/bin/python3

import re
import sys

def parse_hhmmss(s: str):
	parts = s.strip().split(":")
	try:
		if len(parts) == 2:
			mm, ss = int(parts[0]), int(parts[1])
			return mm * 60 + ss
		if len(parts) == 3:
			hh, mm, ss = int(parts[0]), int(parts[1]), int(parts[2])
			return hh * 3600 + mm * 60 + ss
	except ValueError:
		return None
	return None

def mem_to_gb(s: str):
	s = s.strip()
	m = re.fullmatch(r"([0-9.]+)\s*([kKmMgG]?[bB]?)", s)
	if not m:
		return None
	val = float(m.group(1))
	unit = m.group(2).lower()
	if unit in ("", "b"):
		return val / 1e9
	if unit in ("kb", "k"):
		return val * 1e3 / 1e9
	if unit in ("mb", "m"):
		return val * 1e6 / 1e9
	if unit in ("gb", "g"):
		return val
	return None

def logfile_get_time_mem(fn):
	f = open(fn)
	secs = None
	mem_gb = None
	for line in f:
		if line.startswith("Elapsed time "):
			secs = parse_hhmmss(line[len("Elapsed time "):].strip())
		elif line.startswith("Max memory "):
			mem_gb = mem_to_gb(line[len("Max memory "):].strip())
	return secs, mem_gb

def main():
	for fn in sys.argv[1:]:
		secs, mem_gb = logfile_get_time_mem(fn)
	print(secs, mem_gb, fn)
	return 0

if __name__ == "__main__":
	sys.exit(main())
