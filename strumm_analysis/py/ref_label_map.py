#!/usr/bin/python3

import os
import sys

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)
DefaultRefFa = os.path.join(PalmDir, "reference_data", "ref.fa")

def Die(msg):
	sys.stderr.write(msg + "\n")
	sys.exit(1)

def BaseOf(label):
	if label.endswith("-rdrp"):
		return label[:-len("-rdrp")]
	if label.endswith("-decoy"):
		return label[:-len("-decoy")]
	return label

def ReadSuffixMap(ref_fa_path=DefaultRefFa):
	"""Map base chain ID (7f0s_A) to suffixed label (7f0s_A-rdrp)."""
	base_to_suffixed = {}
	for Line in open(ref_fa_path):
		Line = Line.rstrip("\n")
		if not Line.startswith(">"):
			continue
		suffixed = Line[1:].split()[0]
		base = BaseOf(suffixed)
		if base in base_to_suffixed:
			Die("duplicate base in %s: %s and %s" %
				(ref_fa_path, base_to_suffixed[base], suffixed))
		base_to_suffixed[base] = suffixed
	if len(base_to_suffixed) == 0:
		Die("no sequences in %s" % ref_fa_path)
	return base_to_suffixed

def ApplySuffix(label, suffix_map):
	if label.endswith("-rdrp") or label.endswith("-decoy"):
		return label
	if label not in suffix_map:
		Die("no suffix mapping for %s" % label)
	return suffix_map[label]

PfamGroup680 = "PF00680"
PfamGroup946 = "PF00946"

# Diagnostic families first so similar architectures sit together.
RdrpFamilyOrder = (
	"PF00680",
	"PF00946",
	"PF00972",
	"PF00978",
	"PF00998",
	"PF02123",
	"PF04196",
	"PF07925",
)
DecoyFamilyOrder = (
	"PF00078",
	"PF17984",
	"PF00990",
	"PF08082",
	"PF08083",
)

def ClassOf(label):
	if label.endswith("-rdrp"):
		return "rdrp"
	if label.endswith("-decoy"):
		return "decoy"
	Die("unknown suffix: %s" % label)

def DisplayLabel(label):
	s = label
	if s.endswith("-rdrp"):
		s = s[:-len("-rdrp")]
	elif s.endswith("-decoy"):
		s = s[:-len("-decoy")]
	else:
		Die("unknown suffix: %s" % label)
	return s

def PdbSortKey(label):
	pdb = DisplayLabel(label)
	if "_" in pdb:
		code, chain = pdb.split("_", 1)
	else:
		code, chain = pdb, ""
	return (code.lower(), chain)

def PfamsOf(label, domains):
	return set(d[0] for d in domains.get(label, []))

def PfamArchitecture(label, domains):
	"""N-to-C Pfam accessions, first occurrence of each."""
	accs = []
	seen = set()
	for pfamid, start, end, length in sorted(
			domains.get(label, []), key=lambda d: (d[1], d[2], d[0])):
		if pfamid in seen:
			continue
		seen.add(pfamid)
		accs.append(pfamid)
	return tuple(accs)

def FamilyGroup(label, domains, family_order):
	pfams = PfamsOf(label, domains)
	for i, acc in enumerate(family_order):
		if acc in pfams:
			return i
	return len(family_order)

def RdrpPfamGroup(label, domains):
	return FamilyGroup(label, domains, RdrpFamilyOrder)

def PfamSortKey(label, domains, family_order):
	return (
		FamilyGroup(label, domains, family_order),
		PfamArchitecture(label, domains),
		PdbSortKey(label),
	)

def SortTipOrder(labels, domains):
	rdrp = [L for L in labels if ClassOf(L) == "rdrp"]
	decoy = [L for L in labels if ClassOf(L) == "decoy"]
	rdrp.sort(key=lambda L: PfamSortKey(L, domains, RdrpFamilyOrder))
	decoy.sort(key=lambda L: PfamSortKey(L, domains, DecoyFamilyOrder))
	return rdrp + decoy
