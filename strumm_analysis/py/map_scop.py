#!/usr/bin/python3
"""Crosswalk SCOP 1.75 domains/families to SCOP 2022-06-29 FA/SF IDs.

Hard-coded inputs under C:/data/scop/{v1.75,2022-06-29}.
Writes to palmfinder2/data/ by default:

  old_dom_to_new_dom.tsv       old_dom  FA-DOMID  SF-DOMID
  old_family_to_new_family.tsv old_family  CL.CF.SF.FA
  scop_map_report.txt          parseable match log
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import Counter, defaultdict

ScriptDir = os.path.dirname(os.path.abspath(__file__))
PalmDir = os.path.dirname(ScriptDir)
DefaultOutDir = os.path.join(PalmDir, "data")


def _scop_root():
	candidates = [
		r"C:\data\scop",
		"/mnt/c/data/scop",
		os.path.join(os.path.expanduser("~"), "data", "scop"),
	]
	for root in candidates:
		if os.path.isdir(os.path.join(root, "v1.75")) and os.path.isdir(
				os.path.join(root, "2022-06-29")):
			return root
	return candidates[0]


ScopRoot = _scop_root()
OldDir = os.path.join(ScopRoot, "v1.75")
NewDir = os.path.join(ScopRoot, "2022-06-29")

CoordSlack = 5  # residues for fuzzy boundary match


def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)


def NormName(s):
	s = (s or "").strip().lower()
	s = re.sub(r"[^a-z0-9]+", " ", s)
	return " ".join(s.split())


def TokenSet(s):
	return set(NormName(s).split()) - {"", "like", "domain", "protein", "the", "of", "and"}


def NameSimilarity(a, b):
	"""Return ('exact'|'fuzzy'|None, detail)."""
	na, nb = NormName(a), NormName(b)
	if na == "" or nb == "":
		return None, ""
	if na == nb:
		return "exact", ""
	ta, tb = TokenSet(a), TokenSet(b)
	if len(ta) == 0 or len(tb) == 0:
		return None, ""
	inter = ta & tb
	union = ta | tb
	jacc = float(len(inter)) / float(len(union))
	# Require strong overlap for fuzzy name.
	if jacc >= 0.75 and len(inter) >= 2:
		return "fuzzy", "jaccard=%.2f shared=%s" % (jacc, ",".join(sorted(inter)[:8]))
	if jacc >= 0.9 and len(inter) >= 1:
		return "fuzzy", "jaccard=%.2f shared=%s" % (jacc, ",".join(sorted(inter)[:8]))
	# Substring after normalization for short names.
	if len(na) >= 8 and (na in nb or nb in na):
		return "fuzzy", "substring"
	return None, ""


def ParseRangeField(field):
	"""Parse SCOP range field into list of (chain, lo_or_None, hi_or_None).

	Examples: 'A:', 'A:2-124', 'B:,A:', 'A:37-243,A:353-401'
	"""
	field = (field or "").strip()
	if field == "" or field == "-":
		return []
	out = []
	for part in field.split(","):
		part = part.strip()
		if part == "":
			continue
		if ":" not in part:
			# bare chain?
			out.append((part, None, None))
			continue
		chain, rest = part.split(":", 1)
		chain = chain.strip()
		rest = rest.strip()
		if rest == "":
			out.append((chain, None, None))
			continue
		m = re.match(r"^(-?\d+)\s*-\s*(-?\d+)$", rest)
		if m:
			lo, hi = int(m.group(1)), int(m.group(2))
			if hi < lo:
				lo, hi = hi, lo
			out.append((chain, lo, hi))
		else:
			# insertion codes etc. — treat as whole chain segment
			out.append((chain, None, None))
	return out


def RangesOverlapExact(old_segs, new_segs):
	"""True if any same-chain pair overlaps exactly (numeric) or both whole-chain."""
	for oc, olo, ohi in old_segs:
		for nc, nlo, nhi in new_segs:
			if oc.upper() != nc.upper():
				continue
			if olo is None or nlo is None:
				# whole-chain on either side counts as exact chain co-location
				return True
			if max(olo, nlo) <= min(ohi, nhi):
				return True
	return False


def RangesOverlapFuzzy(old_segs, new_segs, slack=CoordSlack):
	"""Return (ok, slack_used_or_None) for near-overlapping same-chain intervals."""
	best_slack = None
	ok = False
	for oc, olo, ohi in old_segs:
		for nc, nlo, nhi in new_segs:
			if oc.upper() != nc.upper():
				continue
			if olo is None or nlo is None:
				return True, 0
			if max(olo, nlo) <= min(ohi, nhi):
				return True, 0
			if ohi < nlo:
				gap = nlo - ohi
			elif nhi < olo:
				gap = olo - nhi
			else:
				gap = 0
			dlo = abs(olo - nlo)
			dhi = abs(ohi - nhi)
			cand = min(gap, max(dlo, dhi))
			if cand <= slack:
				ok = True
				if best_slack is None or cand < best_slack:
					best_slack = cand
	return ok, best_slack


def FormatSegs(segs):
	parts = []
	for c, lo, hi in segs:
		if lo is None:
			parts.append("%s:" % c)
		else:
			parts.append("%s:%d-%d" % (c, lo, hi))
	return ",".join(parts) if parts else "."


def ParseScopCla(fn):
	"""Return (fa_info, fa_ids, sf_ids).

	fa_info[fa_domid] = dict(pdb, reg, segs, sf_domid, cl, cf, sf, fa, cla)
	"""
	fa_info = {}
	fa_ids = set()
	sf_ids = set()
	for Line in open(fn, encoding="utf-8", errors="replace"):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		Fields = Line.split()
		if len(Fields) < 11:
			continue
		fa_dom = Fields[0]
		fa_pdb = Fields[1].upper()
		fa_reg = Fields[2]
		sf_dom = Fields[5]
		cla = Fields[-1]
		cl = cf = sf = fa = None
		for part in cla.split(","):
			if part.startswith("CL="):
				cl = part[3:]
			elif part.startswith("CF="):
				cf = part[3:]
			elif part.startswith("SF="):
				sf = part[3:]
			elif part.startswith("FA="):
				fa = part[3:]
		fa_ids.add(fa_dom)
		sf_ids.add(sf_dom)
		fa_info[fa_dom] = {
			"pdb": fa_pdb,
			"reg": fa_reg,
			"segs": ParseRangeField(fa_reg),
			"sf_domid": sf_dom,
			"cl": cl or ".",
			"cf": cf or ".",
			"sf": sf or ".",
			"fa": fa or ".",
			"cla": cla,
		}
	return fa_info, fa_ids, sf_ids


def ParseRepresented(fn, fa_ids, sf_ids):
	"""Return (by_pc, by_pdb).

	by_pc[(pdb, chain)] -> {'fa': set, 'sf': set}
	by_pdb[pdb] -> {'fa': set, 'sf': set}
	"""
	by_pc = defaultdict(lambda: {"fa": set(), "sf": set()})
	by_pdb = defaultdict(lambda: {"fa": set(), "sf": set()})
	for Line in open(fn, encoding="utf-8", errors="replace"):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		Fields = Line.split()
		if len(Fields) < 3:
			continue
		dom, pdb, chain = Fields[0], Fields[1].upper(), Fields[2]
		key = (pdb, chain)
		if dom in fa_ids:
			by_pc[key]["fa"].add(dom)
			by_pdb[pdb]["fa"].add(dom)
		if dom in sf_ids:
			by_pc[key]["sf"].add(dom)
			by_pdb[pdb]["sf"].add(dom)
	return by_pc, by_pdb


def ParseOldCla(fn):
	"""List of dicts: sid, pdb, chain_field, segs, family, chains."""
	rows = []
	for Line in open(fn, encoding="utf-8", errors="replace"):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		Fields = Line.rstrip("\n").split("\t")
		if len(Fields) < 4:
			continue
		sid = Fields[0].strip()
		pdb = Fields[1].strip().upper()
		rng = Fields[2].strip()
		family = Fields[3].strip()
		segs = ParseRangeField(rng)
		chains = sorted(set(c for c, _lo, _hi in segs))
		if len(chains) == 0 and ":" in rng:
			chains = [rng.split(":", 1)[0]]
		rows.append({
			"sid": sid,
			"pdb": pdb,
			"range": rng,
			"segs": segs,
			"family": family,
			"chains": chains,
		})
	return rows


def ParseOldDesFamilies(fn):
	"""old_family sccs -> description (fa level only)."""
	names = {}
	for Line in open(fn, encoding="utf-8", errors="replace"):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		Fields = Line.rstrip("\n").split("\t")
		if len(Fields) < 5:
			continue
		if Fields[1] != "fa":
			continue
		sccs = Fields[2].strip()
		desc = Fields[4].strip()
		if sccs not in names:
			names[sccs] = desc
	return names


def ParseNewDes(fn):
	"""node_id -> name; also fa_nodes set (ids starting with 4 and length 7)."""
	names = {}
	fa_nodes = {}
	for Line in open(fn, encoding="utf-8", errors="replace"):
		if Line.startswith("#") or len(Line.strip()) == 0:
			continue
		parts = Line.split(None, 1)
		if len(parts) < 2:
			continue
		nid, name = parts[0], parts[1].strip()
		names[nid] = name
		# Family nodes are 4xxxxxx in SCOP2
		if nid.startswith("4") and len(nid) == 7 and nid.isdigit():
			fa_nodes[nid] = name
	return names, fa_nodes


def NewFamilyString(info):
	return "%s.%s.%s.%s" % (info["cl"], info["cf"], info["sf"], info["fa"])


def MapDomains(old_rows, fa_info, represented, represented_pdb, report):
	"""Return list of (old_sid, fa_dom, sf_dom) and fill report lines."""
	out_rows = []
	n_exact = n_fuzzy = n_none = 0
	for row in old_rows:
		sid = row["sid"]
		pdb = row["pdb"]
		old_segs = row["segs"]
		hits = []  # (tier, fa_dom, sf_dom, match, detail)
		# Collect candidate FA-DOMIDs from all chains of this old domain
		cands = set()
		for ch in row["chains"]:
			rec = represented.get((pdb, ch))
			if rec:
				cands |= rec["fa"]
		# PDB-level fallback if chain key missing
		if len(cands) == 0:
			rec = represented_pdb.get(pdb)
			if rec:
				cands |= rec["fa"]

		for fa_dom in sorted(cands):
			info = fa_info.get(fa_dom)
			if info is None:
				continue
			sf_dom = info["sf_domid"] if info["sf_domid"] else "."
			new_segs = info["segs"]
			if RangesOverlapExact(old_segs, new_segs):
				hits.append((0, fa_dom, sf_dom, "exact",
					"old_range=%s new_range=%s" % (
						FormatSegs(old_segs) or row["range"],
						FormatSegs(new_segs) or info["reg"])))
				continue
			ok, slack = RangesOverlapFuzzy(old_segs, new_segs)
			if ok:
				hits.append((1, fa_dom, sf_dom, "fuzzy_coord",
					"old_range=%s new_range=%s slack=%s" % (
						FormatSegs(old_segs) or row["range"],
						FormatSegs(new_segs) or info["reg"],
						slack if slack is not None else ".")))
				continue
			# chain-only fallback when old has whole-chain or no numeric range
			whole = all(lo is None for _c, lo, _hi in old_segs) if old_segs else True
			if whole and len(cands) > 0:
				hits.append((2, fa_dom, sf_dom, "fuzzy_coord",
					"old_range=%s new_range=%s detail=chain_only" % (
						FormatSegs(old_segs) or row["range"],
						FormatSegs(new_segs) or info["reg"])))

		if len(hits) == 0:
			n_none += 1
			report.append(
				"DOM old=%s match=none detail=no_pdb_hit pdb=%s range=%s" %
				(sid, pdb, row["range"]))
			continue

		best_tier = min(h[0] for h in hits)
		kept = [h for h in hits if h[0] == best_tier]
		seen_fa = set()
		for tier, fa_dom, sf_dom, match, detail in kept:
			if fa_dom in seen_fa:
				continue
			seen_fa.add(fa_dom)
			out_rows.append((sid, fa_dom, sf_dom if sf_dom else "."))
			ch = row["chains"][0] if row["chains"] else "."
			report.append(
				"DOM old=%s new_fa=%s new_sf=%s match=%s pdb=%s chain=%s detail=%s" %
				(sid, fa_dom, sf_dom if sf_dom else ".", match, pdb, ch, detail))
			if match == "exact":
				n_exact += 1
			else:
				n_fuzzy += 1
	return out_rows, n_exact, n_fuzzy, n_none


def MapFamilies(old_fam_names, dom_map_rows, old_rows, fa_info, new_fa_names,
		report):
	"""Return list of (old_family, new_family_str)."""
	# sid -> family
	sid_fam = {r["sid"]: r["family"] for r in old_rows}
	# old_family -> Counter of new FA node ids from domain maps
	from_doms = defaultdict(Counter)
	for sid, fa_dom, _sf in dom_map_rows:
		fam = sid_fam.get(sid)
		if fam is None:
			continue
		info = fa_info.get(fa_dom)
		if info is None or info["fa"] == ".":
			continue
		from_doms[fam][info["fa"]] += 1

	# Build reverse: FA node -> example fa_info for CL.CF.SF.FA
	fa_to_info = {}
	for fa_dom, info in fa_info.items():
		fa_node = info["fa"]
		if fa_node != "." and fa_node not in fa_to_info:
			fa_to_info[fa_node] = info

	out = []
	n_exact = n_fuzzy = n_none = 0
	all_old_fams = sorted(set(old_fam_names.keys()) | set(from_doms.keys()))

	for old_fam in all_old_fams:
		old_name = old_fam_names.get(old_fam, "")
		chosen = []  # (new_fam_str, match, detail)

		# Domain-derived FAs (prefer majority)
		if old_fam in from_doms and len(from_doms[old_fam]) > 0:
			counts = from_doms[old_fam]
			best_n = max(counts.values())
			for fa_node, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
				if n < best_n and len(counts) > 1:
					# keep only majority unless unique
					continue
				info = fa_to_info.get(fa_node)
				if info is None:
					continue
				new_s = NewFamilyString(info)
				# Check name agreement
				new_name = new_fa_names.get(fa_node, "")
				kind, det = NameSimilarity(old_name, new_name) if old_name and new_name else (None, "")
				if kind == "exact":
					chosen.append((new_s, "exact_name",
						"name=%s support=%d" % (old_name.replace(" ", "_"), n)))
				elif kind == "fuzzy":
					chosen.append((new_s, "fuzzy_name",
						"old_name=%s new_name=%s %s support=%d" % (
							old_name.replace(" ", "_"),
							new_name.replace(" ", "_"), det, n)))
				else:
					chosen.append((new_s, "from_domains",
						"fa=%s support=%d" % (fa_node, n)))

		# If nothing from domains, try name-only search
		if len(chosen) == 0 and old_name:
			name_hits = []
			for fa_node, new_name in new_fa_names.items():
				kind, det = NameSimilarity(old_name, new_name)
				if kind is None:
					continue
				info = fa_to_info.get(fa_node)
				if info is None:
					continue
				name_hits.append((0 if kind == "exact" else 1, NewFamilyString(info),
					kind + "_name",
					"old_name=%s new_name=%s %s" % (
						old_name.replace(" ", "_"),
						new_name.replace(" ", "_"), det)))
			if name_hits:
				name_hits.sort()
				best = name_hits[0][0]
				for tier, new_s, match, detail in name_hits:
					if tier == best:
						chosen.append((new_s, match, detail))

		if len(chosen) == 0:
			n_none += 1
			report.append(
				"FAM old=%s match=none detail=no_map name=%s" %
				(old_fam, old_name.replace(" ", "_") if old_name else "."))
			continue

		seen = set()
		for new_s, match, detail in chosen:
			if new_s in seen:
				continue
			seen.add(new_s)
			out.append((old_fam, new_s))
			report.append(
				"FAM old=%s new=%s match=%s detail=%s" %
				(old_fam, new_s, match, detail))
			if match.startswith("exact"):
				n_exact += 1
			elif match.startswith("fuzzy"):
				n_fuzzy += 1
			else:
				n_exact += 1  # from_domains counted as resolved
	return out, n_exact, n_fuzzy, n_none


def Main():
	AP = argparse.ArgumentParser(description=__doc__)
	AP.add_argument("--outdir", default=DefaultOutDir,
		help="Output directory (default: palmfinder2/data)")
	Args = AP.parse_args()
	outdir = Args.outdir
	if not os.path.isdir(outdir):
		os.makedirs(outdir)

	old_cla = os.path.join(OldDir, "dir.cla.scop.txt")
	old_des = os.path.join(OldDir, "dir.des.scop.txt")
	new_cla = os.path.join(NewDir, "scop-cla-latest.txt")
	new_des = os.path.join(NewDir, "scop-des-latest.txt")
	new_rep = os.path.join(NewDir, "scop-represented-structures-latest.txt")
	for fn in (old_cla, old_des, new_cla, new_des, new_rep):
		if not os.path.isfile(fn):
			Die("missing %s" % fn)

	sys.stderr.write("parsing new cla...\n")
	fa_info, fa_ids, sf_ids = ParseScopCla(new_cla)
	sys.stderr.write("  FA-DOMIDs=%d SF-DOMIDs=%d\n" % (len(fa_ids), len(sf_ids)))
	sys.stderr.write("parsing represented structures...\n")
	represented, represented_pdb = ParseRepresented(new_rep, fa_ids, sf_ids)
	sys.stderr.write("  pdb-chain keys=%d pdbs=%d\n" %
		(len(represented), len(represented_pdb)))
	sys.stderr.write("parsing old cla/des...\n")
	old_rows = ParseOldCla(old_cla)
	old_fam_names = ParseOldDesFamilies(old_des)
	_new_names, new_fa_names = ParseNewDes(new_des)
	sys.stderr.write("  old domains=%d old families=%d new FA names=%d\n" %
		(len(old_rows), len(old_fam_names), len(new_fa_names)))

	report = []
	sys.stderr.write("mapping domains...\n")
	dom_rows, d_ex, d_fz, d_no = MapDomains(
		old_rows, fa_info, represented, represented_pdb, report)
	sys.stderr.write("mapping families...\n")
	fam_rows, f_ex, f_fz, f_no = MapFamilies(
		old_fam_names, dom_rows, old_rows, fa_info, new_fa_names, report)

	dom_path = os.path.join(outdir, "old_dom_to_new_dom.tsv")
	fam_path = os.path.join(outdir, "old_family_to_new_family.tsv")
	rep_path = os.path.join(outdir, "scop_map_report.txt")

	with open(dom_path, "w", encoding="utf-8") as fh:
		fh.write("old_dom\tFA-DOMID\tSF-DOMID\n")
		for sid, fa_dom, sf_dom in sorted(dom_rows):
			fh.write("%s\t%s\t%s\n" % (sid, fa_dom, sf_dom))

	with open(fam_path, "w", encoding="utf-8") as fh:
		fh.write("old_family\tnew_family\n")
		for old_f, new_f in sorted(fam_rows):
			fh.write("%s\t%s\n" % (old_f, new_f))

	with open(rep_path, "w", encoding="utf-8") as fh:
		fh.write("# SCOP 1.75 -> 2022-06-29 map report\n")
		fh.write("# DOM old=... new_fa=... new_sf=... match=exact|fuzzy_coord|none ...\n")
		fh.write("# FAM old=... new=... match=exact_name|fuzzy_name|from_domains|none ...\n")
		for line in report:
			fh.write(line + "\n")

	sys.stderr.write("wrote %s (%d rows)\n" % (dom_path, len(dom_rows)))
	sys.stderr.write("wrote %s (%d rows)\n" % (fam_path, len(fam_rows)))
	sys.stderr.write("wrote %s (%d lines)\n" % (rep_path, len(report)))
	sys.stderr.write(
		"DOM exact=%d fuzzy=%d none=%d | FAM exact/from=%d fuzzy=%d none=%d\n" %
		(d_ex, d_fz, d_no, f_ex, f_fz, f_no))


if __name__ == "__main__":
	Main()
