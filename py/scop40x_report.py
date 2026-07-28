#!/usr/bin/env python3

import sys
from logfile_get_acc import *
from logfile_get_time_mem import *
from argsort import argsort

def read_search_log(fn):
	flds = fn.split('/')[1].split('_')
	if flds[0] == "reseek28":
		assert len(flds) == 3
		ver, os, mode = flds
		method = ver + "." + os
	else:
		assert len(flds) == 4
		ver, os, mode, truth = flds
		method = ver + "." + os + "." + mode + "." + truth
	secs, mem_gb = logfile_get_time_mem(fn)
	return method, secs, mem_gb

def read_ana_log(fn):
	flds = fn.split('/')[1].split('_')
	if flds[0] == "reseek28":
		assert len(flds) == 4
		ver, os, mode, truth1 = flds
		method = ver + "." + os
	else:
		assert len(flds) == 4
		ver, os, mode, truth1 = flds
		method = ver + "." + os + "." + mode + "." + truth1
	acc, truth = logfile_get_acc(fn)
	return method, acc, truth

methods = set()
truths = set()
secsd = {}
mem_gbd = {}
accd = {}

for fn in sys.argv[1:]:
	if fn.startswith("search_log_scop40x/"):
		method, secs, mem_gb = read_search_log(fn)
		if not method is None and not secs is None and not mem_gb is None:
			methods.add(method)
			secsd[method] = secs
			mem_gbd[method] = mem_gb
	elif fn.startswith("ana_log_scop40x/"):
		method, acc, truth = read_ana_log(fn)
		if not method is None and not acc is None and not truth is None:
			truths.add(truth)
			methods.add(method)
			accd[(method, truth)] = acc
	else:
		assert False, "Expected search_log_scop40x/ or ana_log_scop40x/"

def get_acc_float(acc_str):
	regex = "[ST][uo][mp]3=[0-9Ee.]*"
	M = re.search(regex, acc_str)
	if not M is None:
		s = M.group()
		return float(s[5:])

for mode in [ "fast", "sensitive" ]:
	print()
	print("==== time " + mode)
	s = "%7.7s" % "Secs"
	s += "  %7.7s" % "Pct"
	s += "  Method"
	print(s)

	secs_list = []
	methods_list = []
	min_secs = 9e9
	for method in methods:
		if method.find(mode) < 0:
			continue
		secs = secsd.get(method)
		if secs is None:
			continue
		min_secs = min(min_secs, secs)
		methods_list.append(method)
		secs_list.append(secs)
		
	order = argsort(secs_list, reverse=True)
	for k in order:
		secs = secs_list[k]
		method = methods_list[k]
		pct = 100*(secs/min_secs - 1)
		s = "%7u" % secs
		if secs == min_secs:
			s += "  %7.7s" % ""
		else:
			s += "  %+6.2f%%" % pct
		s += "  " + method
		print(s)

for mode in [ "fast", "sensitive" ]:
	print()
	print("==== memory " + mode)
	s = "%7.7s" % "Gb"
	s += "  %7.7s" % "Pct"
	s += "  Method"
	print(s)

	gb_list = []
	methods_list = []
	min_gb = 9e9
	for method in methods:
		if method.find(mode) < 0:
			continue
		gb = mem_gbd.get(method)
		if gb is None:
			continue
		min_gb = min(min_gb, gb)
		methods_list.append(method)
		gb_list.append(gb)
		
	order = argsort(gb_list, reverse=True)
	for k in order:
		gb = gb_list[k]
		method = methods_list[k]
		pct = 100*(gb/min_gb - 1)
		s = "%7u" % gb
		if gb == min_gb:
			s += "  %7.7s" % ""
		else:
			s += "  %+6.2f%%" % pct
		s += "  " + method
		print(s)

for truth in truths:
	print()
	print("====  " + truth)
	method_list = []
	acc_list = []
	xxx3_list = []
	for method in methods:
		mt = (method, truth)
		acc = accd.get(mt)
		if acc is None:
			continue
		acc_float = get_acc_float(acc)
		xxx3_list.append(acc_float)
		method_list.append(method)
		acc_list.append(acc)
	order = argsort(xxx3_list)
	for k in order:
		method = method_list[k]
		acc = acc_list[k]
		print(acc, method)
