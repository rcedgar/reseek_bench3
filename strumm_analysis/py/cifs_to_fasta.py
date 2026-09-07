#!/usr/bin/python3

import argparse
import sys
from pathlib import Path

from cif_seq import (
	ChainSequence,
	Die,
	FormatFasta,
	MakeQuery,
	ParseChainSequences,
	ReadPathList,
)

Usage = \
(
"Convert mmCIF files to FASTA using the shared cif_seq polymer sequence "
"construction (same as map_motifs_to_cifs / cif_make_sifts_annotations). "
"Input is a text file with one CIF pathname per line, optionally followed "
"by TAB and chain id. Writes FASTA to stdout."
)

AP = argparse.ArgumentParser(description=Usage)
AP.add_argument("input",
	help="Text file: one CIF path per line; optional TAB + chain id")
AP.add_argument("--chainchar", default="_",
	help="Separator between stem and chain in FASTA id (default '_')")
AP.add_argument("--width", type=int, default=80,
	help="FASTA line width (default 80; 0 = no wrap)")
Args = AP.parse_args()


def ProcessFile(path, chain_filter, chainchar, width):
	cif_path = Path(path)
	if not cif_path.is_file():
		Die("not a file: %s" % path)
	chains, label_to_auth = ParseChainSequences(
		cif_path, chain_filter=chain_filter)
	if len(chains) == 0:
		if chain_filter is not None:
			sys.stderr.write(
				"warning: chain %s not in %s\n" % (chain_filter, path))
		else:
			sys.stderr.write("warning: no polymer chains in %s\n" % path)
		return
	stem = cif_path.stem
	chain_ids = sorted(chains.keys())
	for chain in chain_ids:
		residues = chains[chain]
		seq = ChainSequence(residues)
		if len(seq) == 0:
			sys.stderr.write(
				"warning: empty sequence for %s chain %s\n" % (path, chain))
			continue
		header = MakeQuery(stem, chain, chainchar)
		sys.stderr.write(header + '\n')
		sys.stdout.write(FormatFasta(header, seq, width))
		sys.stdout.write("\n")


def main():
	entries = ReadPathList(Args.input)
	if len(entries) == 0:
		Die("no input paths in %s" % Args.input)
	width = Args.width if Args.width > 0 else 0
	for path, chain_filter, _lineno in entries:
		sys.stderr.write(path + '\n')
		ProcessFile(path, chain_filter, Args.chainchar, width)


if __name__ == "__main__":
	main()
