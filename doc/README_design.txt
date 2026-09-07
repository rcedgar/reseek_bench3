Homology search validation
==========================

Self-contained scripts for comparative validation of protein homology search
algorithms on SCOP or CATH all-vs-all benchmarks.  Designed for large hit
lists (~10^8 pairs) by condensing raw hits into small summary files that
support fast plot regeneration.

Reference databases
-------------------
  CATH40    Unedited v.4.4.0 https://download.cathdb.info/cath/releases/all-releases/v4_4_0/non-redundant-data-sets/cath-dataset-nonredundant-S20-v4_4_0.pdb.tgz
  SCOP40    Unedited SCOP40 v1.75 http://scop.mrc-lmb.cam.ac.uk/legacy/parse/dir.*.scop.txt_1.75
  SCOP40x   Release 0fc3d7 of curated SCOP40 v1.75 https://github.com/rcedgar/scop40c

Overview
--------
Processing is split into three stages:
  Stage 0   build_derived_info.py   lookup -> derived_info JSON (once)

  Stage 1   hits_to_edf.py          hits -> pair-list summary (.edf)
            hits_to_topcat.py       hits -> leave-category-out summary (.tcat)

  Stage 2   plot_curve.py           summaries -> one CVE / ROC / PR panel (SVG)

Stage 1 reads the large hits file (slow).
Stage 2 reads only the small summary files (fast), so that axis scales,
algorithm subsets, and plot ranges can be adjusted for publication without
re-processing hits.

Summary statistics (Sum3, Top3, etc.) are reported separately; plot_curve.py
generates figures only.


Input files
-----------

Lookup TSV (no header)
  Column 1: domain_id          e.g. d1asha_
  Column 2: family_id          e.g. a.1.1.2
                               (class.fold.superfamily.family)

Hits TSV (no header)
  One row per reported pair.  Column layout varies by algorithm.
  Three fields are used (1-based indices, default 1,2,3):
    query    domain label (may include extra text after '/', which is stripped)
    target   domain label
    score    alignment score or E-value

  Self-hits (query == target after stripping) are always ignored.

  Domains not present in the lookup file are silently ignored.

  Score direction (required):
    --evalue   lower score = stronger hit (E-values, P-values)
    --score    higher score = stronger hit (bit scores, similarity)
    Exactly one of --evalue or --score must be set; no auto-detection.


Directory contents
------------------
  README.txt                    This file (operation, file formats, workflow)
  truth_standards.md            Authoritative definitions of truth standards
                                and all metrics (coverage, error, precision,
                                recall, Sum3, Top3, ROC).  This README does not
                                repeat those definitions.

  info/                         Lookup TSV files (domain_id, family_id)
  derived_info/                 Stage 0 JSON output (denominators, query sets)
  test/                         Synthetic regression fixtures

  bash/                         Bash scripts

  py/                           Python scripts
    lookup.py                   Shared helper to load domain/fold/superfamily/family maps from lookup TSV
    common.py                   Shared helpers (score comparison, field parsing)
    build_derived_info.py       Stage 0: precompute denominators and query sets
    hits_to_edf.py              Stage 1a: pair-list condensation
    hits_to_topcat.py           Stage 1b: leave-category-out top-hit condensation
    plot_curve.py               Stage 2: one plot panel from summary files
    plot_epq.py                 E-value calibration plot (FPEPQ vs E-value)


Stage 0: build_derived_info.py
------------------------------

  python build_derived_info.py --lookup ../info/scop40x.lookup \
      --output ../derived_info/scop40x.json

Or run bash/stage0_generate_derived_info.bash for all reference databases.

Reads the lookup file and writes a JSON summary containing:

  ndom                    number of domains in lookup
  query_set_topfold       domain ids eligible for topfold analysis
  query_set_topsf         domain ids eligible for topsf analysis
  N_possible_tp.<truth>   per-truth-standard quantities consumed by Stage 1
  N_possible_fp.<truth>
  N_possible_tp.topfold
  N_possible_tp.topsf

The truth standards, query sets, and the meaning of these denominators are
defined in truth_standards.md.


Stage 1a: hits_to_edf.py  (pair-list analysis)
-----------------------------------------------
Examples:
  python hits_to_edf.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json \
      --truth superfamily --score --output reseek.superfamily.edf

  python hits_to_edf.py --hits foldseek.scop40.tsv --fields 1,2,3 \
      --derived-info ../derived_info/scop40x.json \
      --evalue --truth superfamilyx --output foldseek.superfamilyx.edf

Reads denominators (ndom, N_possible_tp, N_possible_fp) from derived_info JSON.
Loads domain annotations from the lookup TSV for per-hit TP/FP classification.

Single streaming pass over the hits file.  Each non-ignored hit is classified
as TP or FP according to the truth standard.  Counts are accumulated per
distinct score value (not per pair), producing a compact histogram.

Output: .edf file (TSV with comment header)

  # algo=...
  # reference=...
  # truth=...
  # fields=1,2,3
  # score_direction=higher_better | lower_better
  # ndom=10457
  # N_possible_tp=...
  # N_possible_fp=...
  # SEPQ0.1=... SEPQ1=... SEPQ10=... Sum3=... SFFP=...
  score    n_tp    n_fp    cum_tp    cum_fp

Rows are sorted best-first (descending score, or ascending E-value).
cum_tp and cum_fp are running totals over rows at or above the current
score threshold (i.e. hits at least as good as this score).

Memory: O(unique scores), typically far smaller than the number of hits.


Stage 1b: hits_to_topcat.py  (leave-category-out analysis)
----------------------------------------------------------

  python hits_to_topcat.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json \
      --score --truth topfold --output reseek.topfold.tcat

  python hits_to_topcat.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json \
      --score --truth topsf --output reseek.topsf.tcat

Reads query sets and N_possible_tp from derived_info JSON.

Single streaming pass.  For each domain in the query set, the two per-query
best scores defined in truth_standards.md are tracked (S_tp and S_xf).  See
truth_standards.md for the held-out experiments, the query sets, and how the
tracked scores are turned into TP / FN / TN / FP.

Output: .tcat file (TSV with comment header)

  # truth=topfold | topsf
  # N_possible_tp=...
  # score_direction=...
  # TEPQ0.001=... TEPQ0.01=... TEPQ0.1=... Top3=...
  domain    S_tp    S_xf

  '.' means no qualifying hit was seen for that score.


Stage 2: plot_curve.py
----------------------

One invocation produces one SVG panel.  Multiple --input files overlay
algorithms on the same plot.  All inputs must share the same truth standard
and the same curve type.  Do not mix .edf and .tcat inputs in one run.

  python plot_curve.py --type roc --input reseek.superfamily.edf \
      --input foldseek.superfamily.edf --output roc_superfamily.svg

  python plot_curve.py --type cve --input reseek.superfamily.edf \
      --xscale log --xlim 1e-4,0.1 --output cve_superfamily.svg

  python plot_curve.py --type pr --input reseek.superfamily.edf \
      --output pr_superfamily.svg

  python plot_curve.py --type cve --input reseek.topfold.tcat \
      --output cve_topfold.svg

  python plot_curve.py --type roc --input reseek.topsf.tcat \
      --output roc_topsf.svg

  python plot_curve.py --type pr --input reseek.topsf.tcat \
      --output pr_topsf.svg

Required flags:
  --type cve|roc|pr
  --input PATH          (repeatable)
  --output PATH         (.svg)

Optional flags:
  --xscale linear|log   (default linear)
  --yscale linear|log   (default linear)
  --xlim LO,HI
  --ylim LO,HI
  --ylabel TEXT         override y-axis label
  --nolegend            omit the panel legend
  --legend-only         write a horizontal legend SVG (no curves)
  --decimal-xticks      decimal tick labels (0.001 not 1e-3)
  --decimal-yticks      decimal tick labels (0.001 not 1e-3)

Supported curve types by input kind:
  .edf   CVE, ROC, PR
  .tcat  CVE, ROC, PR

Reads comment headers and data rows from summary files.  Each .edf row yields
one curve point from cum_tp and cum_fp; .tcat curves sweep score thresholds
using S_tp and S_xf (same logic as hits_to_topcat.py).  Legend labels come
from the algo= header field (reference= is metadata only).

Axis definitions (coverage, error, precision, recall, FPR) are in
truth_standards.md.


Stage 2b: plot_epq.py  (E-value calibration)
---------------------------------------------

Plots FPEPQ (mean false positives per query, cum_fp / ndom) against E-value
on log-log axes.  A dotted y=x line shows ideal E-value calibration
(E = FPEPQ).

  python plot_epq.py --input reseek.superfamily.edf --output epq_sf.svg

  python plot_epq.py --input reseek.superfamily.edf --output epq_sf.svg \
      --xlim 1e-10,10 --ylim 1e-10,10

Required: --input (a single .edf file with E-value scores) and --output (SVG).
Optional:
  --xlim LO,HI
  --ylim LO,HI
  --dbsize N      Scores are P-values; E-value = N * P

If --dbsize is set, scores are treated as P-values and converted to E-values
by multiplying by N.  Otherwise, scores must be E-values
(score_direction=lower_better).


Typical workflow
----------------

  # Once per lookup (or bash/stage0_generate_derived_info.bash)
  python build_derived_info.py --lookup ../info/scop40x.lookup \
      --output ../derived_info/scop40x.json

  # Once per algorithm (pair-list, each truth standard of interest)
  python hits_to_edf.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json --score \
      --truth superfamily --output reseek.superfamily.edf
  python hits_to_edf.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json --score \
      --truth superfamilyx --output reseek.superfamilyx.edf

  # Once per algorithm (leave-category-out)
  python hits_to_topcat.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json --score \
      --truth topfold --output reseek.topfold.tcat
  python hits_to_topcat.py --hits reseek.scop40.tsv --lookup ../info/scop40x.lookup \
      --derived-info ../derived_info/scop40x.json --score \
      --truth topsf --output reseek.topsf.tcat

  # Fast plot iteration (repeat as needed)
  python plot_curve.py --type cve --input reseek.superfamily.edf \
      --input foldseek.superfamily.edf --input dali.superfamily.edf \
      --xscale log --xlim 1e-4,0.2 --output fig_cve_sf.svg


Performance notes
-----------------

  ~121M hit pairs: stage 1 takes minutes to tens of minutes depending on
  disk I/O.  RAM stays bounded by the number of distinct score values
  (typically 10^3-10^6), not the hit count.

  If distinct scores approach the hit count (rare continuous scores), use
  the --spill option documented in hits_to_edf.py to write a temporary
  binary file and sort on disk.

  Stage 2 completes in seconds.


Domain label parsing
--------------------

  Hit file label:     d1asha_/a.1.1.2
  Parsed domain_id:   d1asha_         (text before first '/')

  Lookup domain_id:   d1asha_

  Matching is on domain_id only.


Regression testing
------------------

  cd homval/test
  python make_expected_edf.py          # brute-force golden .edf files
  python make_expected_tcat.py         # brute-force golden .tcat files
  python ../py/test_homval.py          # compare with hits_to_edf.py / hits_to_topcat.py
