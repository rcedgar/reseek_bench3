Synthetic regression fixtures for homval
========================================

Files
-----
  test.lookup              5-domain lookup (2 families in SF a.1.1; fold a.1 has 2 SF)
  test_hits.tsv            Hand-crafted hits with --score (higher better)
  derived_info.json        python ../py/build_derived_info.py --lookup test.lookup \
                               --output derived_info.json
  make_expected_edf.py     Brute-force reference (no homval imports); writes expected.*.edf
  expected.superfamily.edf Golden .edf for --truth superfamily
  expected.superfamilyx.edf Golden .edf for --truth superfamilyx

Workflow
--------
  cd homval/test
  python ../py/build_derived_info.py --lookup test.lookup --output derived_info.json
  python make_expected_edf.py
  python ../py/test_homval.py

make_expected_edf.py and hits_to_edf.py must agree on expected.*.edf bodies,
denominators, and Sum3.  If they agree, both implement the same spec.

What test_hits.tsv exercises
----------------------------
  da da              self-hit -> ignored
  da du              unknown target -> ignored
  da db, da dc, db da same superfamily -> TP
  da dd, dc dd, db dd same fold, different superfamily -> FP (superfamily)
                       -> ignored (superfamilyx)
  da de, dd de       different fold -> FP (both truths)
  score 80 tie       da dc (TP) and db dd (FP) share one score bucket

superfamily
===========
reseek -fast_bench_hits test_hits.tsv -lookup test.lookup --qfield 1 -tfield 2 -scorefield 3 -truth sf
SEPQ0.1=0.333 SEPQ1=0.667 SEPQ10=0.667 Sum3=2.333 sf test_hits

# expected.superfamily.edf
# SEPQ0.1=0.333 SEPQ1=0.500 SEPQ10=0.500 Sum3=1.917

superfamilyx
============
reseek -fast_bench_hits test_hits.tsv -lookup test.lookup --qfield 1 -tfield 2 -scorefield 3 -truth sfx
SEPQ0.1=0.000 SEPQ1=1.000 SEPQ10=1.000 Sum3=2.500 sfx test_hits

# expected.superfamilyx.edf
# SEPQ0.1=0.500 SEPQ1=0.500 SEPQ10=0.500 Sum3=2.250