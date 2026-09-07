#!/usr/bin/python3
"""Shared mmCIF polymer sequence construction.

Sequences come from _pdbx_poly_seq_scheme when present (full polymer,
including unobserved residues), else from CA atoms in _atom_site.
Chains are keyed by auth strand id (pdb_strand_id), matching *_cif_chain.files.

Only amino-acid and nucleic-acid polymer residues are kept; other chains
(ligands, water, …) are skipped. Unknown-residue warnings are emitted only
when that chain was explicitly requested via chain_filter.
"""

import sys
from pathlib import Path

MISSING = {".", "?", ""}
ALTLOC_OK = {"", " ", ".", "?", "A"}

AA3_TO_1 = {
	"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B",
	"CYS": "C", "GLN": "Q", "GLU": "E", "GLX": "Z", "GLY": "G",
	"HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
	"MSE": "M", "CAS": "C", "PHE": "F", "PRO": "P", "PYL": "O",
	"SEC": "U", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y",
	"UNK": "X", "VAL": "V",
}

# Modified AA parents still count as amino acids for chain typing.
STANDARD_COMP = {
	"ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX",
	"GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "PYL",
	"SEC", "SER", "THR", "TRP", "TYR", "UNK", "VAL",
}

NT_TO_1 = {
	"A": "A", "C": "C", "G": "G", "U": "U", "I": "I", "T": "T",
	"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DI": "I",
	"ADE": "A", "CYT": "C", "GUA": "G", "THY": "T", "URA": "U",
}


def Die(msg):
	sys.stderr.write("*ERROR* " + msg + "\n")
	sys.exit(1)


def Present(value):
	if value is None:
		return False
	return value.strip() not in MISSING


def Comp1(resname):
	"""Map a component id to (letter, kind, is_gap) or (None, None, False).

	kind is 'aa' or 'nt'. Non-polymer components return (None, None, False).
	"""
	name = (resname or "").strip().upper()
	if name in AA3_TO_1:
		return AA3_TO_1[name], "aa", name not in STANDARD_COMP
	if name in NT_TO_1:
		return NT_TO_1[name], "nt", False
	return None, None, False


def Aa1(resname):
	"""Back-compat: 1-letter code and is_gap; unknown polymer → ('X', True)."""
	letter, kind, is_gap = Comp1(resname)
	if letter is None:
		return "X", True
	return letter, is_gap


def MakeQuery(stem, chain, chainchar="_"):
	if chain is None or chain == "":
		return stem
	return stem + chainchar + chain


def ReadPathList(fn):
	"""Return list of (path, chain_or_None, lineno)."""
	entries = []
	for lineno, line in enumerate(open(fn), 1):
		line = line.rstrip("\n")
		if len(line) == 0 or line.startswith("#"):
			continue
		fields = line.split("\t")
		path = fields[0].strip()
		if path == "":
			continue
		chain = None
		if len(fields) > 1:
			chain = fields[1].strip()
			if chain == "":
				chain = None
		entries.append((path, chain, lineno))
	return entries


class TokenStream:
	def __init__(self, fh):
		self._gen = _cif_tokens(fh)
		self._peek = None

	def peek(self):
		if self._peek is None:
			self._peek = next(self._gen, None)
		return self._peek

	def next(self):
		if self._peek is not None:
			tok = self._peek
			self._peek = None
			return tok
		return next(self._gen, None)


def _tokenize_line(line):
	i = 0
	n = len(line)
	while i < n:
		c = line[i]
		if c.isspace():
			i += 1
			continue
		if c == "#":
			return
		if c in ("'", '"'):
			quote = c
			i += 1
			start = i
			while i < n and line[i] != quote:
				i += 1
			yield line[start:i]
			if i < n:
				i += 1
			continue
		start = i
		while i < n and not line[i].isspace() and line[i] != "#":
			i += 1
		yield line[start:i]
		if i < n and line[i] == "#":
			return


def _cif_tokens(fh):
	while True:
		line = fh.readline()
		if line == "":
			return
		line = line.rstrip("\r\n")
		if not line:
			continue
		if line[0] == ";":
			chunks = [line[1:]]
			while True:
				more = fh.readline()
				if more == "":
					yield "\n".join(chunks)
					return
				more = more.rstrip("\r\n")
				if more.startswith(";"):
					yield "\n".join(chunks)
					rest = more[1:]
					if rest.strip():
						yield from _tokenize_line(rest)
					break
				chunks.append(more)
			continue
		stripped = line.lstrip()
		if stripped.startswith("#"):
			continue
		yield from _tokenize_line(line)


def SplitTag(tag):
	if tag.startswith("_"):
		tag = tag[1:]
	if "." not in tag:
		return tag, ""
	cat, item = tag.split(".", 1)
	return cat, item


def IsCifTag(tok):
	return tok is not None and tok.startswith("_") and "." in tok


def IsCifStop(tok):
	if tok is None:
		return True
	if tok in ("loop_", "stop_"):
		return True
	if tok.startswith("data_") or tok.startswith("save_") or tok.startswith("global_"):
		return True
	return IsCifTag(tok)


def ReadLoopRows(ts, fields):
	items = [SplitTag(f)[1] for f in fields]
	n = len(fields)
	rows = []
	while not IsCifStop(ts.peek()):
		row_vals = []
		for _ in range(n):
			if IsCifStop(ts.peek()) and len(row_vals) == 0:
				return rows
			val = ts.next()
			if val is None:
				Die("unexpected EOF in mmCIF loop")
			row_vals.append(val)
		if len(row_vals) != n:
			Die("incomplete mmCIF loop row")
		rows.append(dict(zip(items, row_vals)))
	return rows


def SkipLoop(ts):
	while not IsCifStop(ts.peek()):
		ts.next()


def IterCifLoops(cif_path):
	"""Yield (category, rows) for each loop_ in the file."""
	with open(cif_path, encoding="utf-8", errors="replace") as fh:
		ts = TokenStream(fh)
		while True:
			tok = ts.next()
			if tok is None:
				break
			if tok != "loop_":
				if IsCifTag(tok):
					ts.next()
				continue
			fields = []
			while IsCifTag(ts.peek()):
				fields.append(ts.next())
			if not fields:
				continue
			cat = SplitTag(fields[0])[0]
			yield cat, ReadLoopRows(ts, fields)


def _strand_matches_filter(strand, asym, chain_filter, label_to_auth):
	if chain_filter is None:
		return True
	if strand == chain_filter or asym == chain_filter:
		return True
	auth = label_to_auth.get(chain_filter)
	if auth is not None and strand == auth:
		return True
	return False


def ResolveChain(chain_id, chains, label_to_auth):
	"""Map a label or auth chain id to an auth key present in chains."""
	if chain_id is None:
		return None
	if chain_id in chains:
		return chain_id
	auth = label_to_auth.get(chain_id)
	if auth is not None and auth in chains:
		return auth
	return None


def _finalize_chains(chains, label_to_auth, chain_kinds, unknown_by_strand,
	cif_path, chain_filter):
	"""Drop non-polymer chains; warn on unknown comps only for explicit chains."""
	keep = {}
	for strand, residues in chains.items():
		if len(residues) == 0:
			continue
		kinds = chain_kinds.get(strand, set())
		if not kinds:
			continue
		keep[strand] = residues
	for strand in keep:
		keep[strand].sort(key=lambda t: t[0])

	if chain_filter is not None:
		auth = ResolveChain(chain_filter, keep, label_to_auth)
		warn_key = auth if auth is not None else chain_filter
		unknown = unknown_by_strand.get(warn_key)
		if unknown is None and auth != chain_filter:
			unknown = unknown_by_strand.get(chain_filter)
		if unknown:
			sys.stderr.write(
				"warning: unknown residue %s in %s chain %s\n"
				% (unknown, cif_path, chain_filter))
		if auth is not None and auth in keep:
			keep = {auth: keep[auth]}
		else:
			keep = {}
	return keep, label_to_auth


def _parse_poly_seq_scheme(rows, cif_path, chain_filter):
	"""Return (chains, label_to_auth).

	chains: auth_strand -> [(label_seq_id, aa, is_gap), ...]
	label_to_auth: label asym_id -> auth strand id
	"""
	chains = {}
	seen = {}
	label_to_auth = {}
	chain_kinds = {}
	unknown_by_strand = {}
	for row in rows:
		asym = (row.get("asym_id") or "").strip()
		strand = (row.get("pdb_strand_id") or "").strip()
		if not Present(strand):
			strand = asym
		if Present(asym) and Present(strand) and asym not in label_to_auth:
			label_to_auth[asym] = strand

	for row in rows:
		asym = (row.get("asym_id") or "").strip()
		strand = (row.get("pdb_strand_id") or "").strip()
		if not Present(strand):
			strand = asym
		if not Present(strand):
			continue
		if not _strand_matches_filter(strand, asym, chain_filter, label_to_auth):
			continue
		seq_id_s = (row.get("seq_id") or "").strip()
		if not Present(seq_id_s):
			continue
		try:
			seq_id = int(seq_id_s)
		except ValueError:
			continue
		resname = (row.get("pdb_mon_id") or "").strip()
		if not Present(resname):
			resname = (row.get("auth_mon_id") or "").strip()
		if not Present(resname):
			resname = (row.get("mon_id") or "").strip()
		if not Present(resname):
			resname = "UNK"
		letter, kind, is_gap = Comp1(resname)
		if letter is None:
			if strand not in unknown_by_strand:
				unknown_by_strand[strand] = resname
			letter, is_gap = "X", True
			kind = None
		key = (strand, seq_id)
		if key in seen:
			continue
		seen[key] = True
		chains.setdefault(strand, []).append((seq_id, letter, is_gap))
		if kind is not None:
			chain_kinds.setdefault(strand, set()).add(kind)
	return _finalize_chains(
		chains, label_to_auth, chain_kinds, unknown_by_strand,
		cif_path, chain_filter)


def _parse_atom_site_ca(rows, cif_path, chain_filter):
	chains = {}
	seen = {}
	first_model = None
	chain_kinds = {}
	unknown_by_strand = {}
	label_to_auth = {}
	for row in rows:
		atom = (row.get("label_atom_id") or row.get("auth_atom_id") or "").strip()
		if atom != "CA":
			continue
		alt = (row.get("label_alt_id") or ".").strip()
		if alt not in ALTLOC_OK:
			continue
		model = row.get("pdbx_PDB_model_num")
		if Present(model):
			m = str(model).strip()
			if first_model is None:
				first_model = m
			elif m != first_model:
				continue
		label_chain = (row.get("label_asym_id") or "").strip()
		chain = (row.get("auth_asym_id") or "").strip()
		if not Present(chain):
			chain = label_chain
		if not Present(chain):
			continue
		if Present(label_chain) and label_chain not in label_to_auth:
			label_to_auth[label_chain] = chain
		if not _strand_matches_filter(chain, label_chain, chain_filter, label_to_auth):
			continue
		seq_id_s = (row.get("label_seq_id") or "").strip()
		if not Present(seq_id_s):
			continue
		try:
			seq_id = int(seq_id_s)
		except ValueError:
			continue
		resname = row.get("label_comp_id") or row.get("auth_comp_id") or "UNK"
		letter, kind, is_gap = Comp1(resname)
		if letter is None:
			if chain not in unknown_by_strand:
				unknown_by_strand[chain] = (resname or "").strip() or "UNK"
			letter, is_gap = "X", True
			kind = None
		key = (chain, seq_id)
		if key in seen:
			continue
		seen[key] = True
		chains.setdefault(chain, []).append((seq_id, letter, is_gap))
		if kind is not None:
			chain_kinds.setdefault(chain, set()).add(kind)
	for c in list(chains.keys()):
		label_to_auth.setdefault(c, c)
	return _finalize_chains(
		chains, label_to_auth, chain_kinds, unknown_by_strand,
		cif_path, chain_filter)


def ParseChainSequences(cif_path, chain_filter=None):
	"""Return (chains, label_to_auth).

	chains: auth_chain -> [(label_seq_id, aa, is_gap), ...]
	label_to_auth: label asym_id -> auth_chain

	If chain_filter is set, only that chain (label or auth id) is returned,
	and unknown-residue warnings apply only to that chain.
	"""
	cif_path = str(cif_path)
	scheme_rows = []
	atom_rows = []
	with open(cif_path, encoding="utf-8", errors="replace") as fh:
		ts = TokenStream(fh)
		while True:
			tok = ts.next()
			if tok is None:
				break
			if tok != "loop_":
				if IsCifTag(tok):
					ts.next()
				continue
			fields = []
			while IsCifTag(ts.peek()):
				fields.append(ts.next())
			if not fields:
				continue
			cat = SplitTag(fields[0])[0]
			if cat == "pdbx_poly_seq_scheme":
				scheme_rows.extend(ReadLoopRows(ts, fields))
			elif cat == "atom_site":
				if scheme_rows:
					SkipLoop(ts)
				else:
					atom_rows.extend(ReadLoopRows(ts, fields))
			else:
				SkipLoop(ts)
	if scheme_rows:
		chains, label_to_auth = _parse_poly_seq_scheme(
			scheme_rows, cif_path, chain_filter)
		if chains or chain_filter is not None:
			return chains, label_to_auth
	return _parse_atom_site_ca(atom_rows, cif_path, chain_filter)


def ChainSequence(residues):
	"""Contiguous 1-letter sequence from residue list."""
	return "".join(aa for _sid, aa, _gap in residues)


def SequenceSlice(residues, start, end):
	"""1-letter sequence for inclusive label_seq_id range; missing -> X."""
	by_id = {sid: aa for sid, aa, _gap in residues}
	return "".join(by_id.get(sid, "X") for sid in range(start, end + 1))


def FormatFasta(header, sequence, width=80):
	lines = [">%s" % header]
	if width and width > 0:
		for i in range(0, len(sequence), width):
			lines.append(sequence[i : i + width])
	else:
		lines.append(sequence)
	return "\n".join(lines)
