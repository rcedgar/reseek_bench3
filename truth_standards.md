# Truth standards and benchmark metrics

This document explains, from the ground up, how `homval` decides whether a
hit reported by a structure- or sequence-search algorithm is "right" or
"wrong", and how those decisions are turned into the numbers and plots we use
to compare algorithms. It is written for a reader who knows what a protein
domain and a homology search are, but who has not seen this particular style
of benchmarking before. No prior exposure to ROC curves, precision/recall, or
"errors per query" is assumed.

If you only remember one thing: **a "truth standard" is just a rule that looks
at two domains and says "these should match", "these should not match", or
"don't count this pair". Everything else is bookkeeping on top of that rule.**


## 1. The raw material

We run an *all-vs-all* search: every domain in a reference set (e.g. SCOP40)
is searched against every other domain. The algorithm reports **hits**, one
per line, in a tab-separated file:

```
query    target    score
d1asha_  d1mbaa_   312.5
d1asha_  d2gdma_   88.1
...
```

- **query** and **target** are domain identifiers. They sometimes carry an
  annotation suffix after a slash (`d1asha_/a.1.1.2`); we strip everything
  from the slash onward and keep `d1asha_`.
- **score** measures how strong the hit is. For most methods a **higher score
  is better**. For methods that report an **E-value** (or P-value), a **lower
  value is better**. `homval` is told which convention applies (`--score` or
  `--evalue`); it never guesses.
- A domain matched against itself (a **self-hit**) is always ignored.
- These files are large: a full all-vs-all on SCOP40 is on the order of
  10^8 (hundreds of millions) of hit lines.

Separately, a **lookup** file tells us the classification of each domain:

```
d1asha_    a.1.1.2
```

The second column is a SCOP/CATH label of the form
`class . fold . superfamily . family`. So `a.1.1.2` means class `a`, fold
`a.1`, superfamily `a.1.1`, family `a.1.1.2`. The hierarchy is nested: domains
in the same family are also in the same superfamily, fold, and class.

The biological intuition behind the hierarchy:

- **Family**: clearly related, usually obvious from sequence alone.
- **Superfamily**: related (common ancestor inferred from structure/function)
  but the sequence similarity may be too weak to see directly. This is the
  interesting regime for structure-based search.
- **Fold**: same overall shape. May reflect common ancestry, or may be
  convergent ("analogous") — shape alone does not prove homology.
- **Class**: very broad (e.g. "all alpha"); not used as a truth standard here.


## 2. The single most important idea: a truth standard

When the algorithm reports a hit between query `q` and target `t`, is that hit
correct? That depends on what we are trying to detect. A **truth standard** is
a rule that classifies the pair `(q, t)` into one of three buckets:

- **True positive (TP)** — the pair *should* be found (they are genuinely
  related at the level we care about).
- **False positive (FP)** — the pair should *not* be found (unrelated at the
  level we care about); reporting it is an error.
- **Ignored** — we deliberately refuse to score this pair, because it sits in
  a gray zone where calling it right or wrong would be unfair.

Self-hits are always ignored, and any pair involving a domain missing from the
lookup is ignored.

The reason we need *several* truth standards is that "related" is not a single
yes/no question — it depends on whether you mean same family, same
superfamily, or same fold. A method that is excellent at finding same-family
relationships may be mediocre at the harder same-superfamily problem.


## 3. Two different questions, two analysis styles

`homval` answers two genuinely different questions, and it is worth being
clear about which is which, because they use the words "true positive",
"sensitivity" and "error" in related but distinct ways.

1. **Pairwise retrieval (the `.edf` / "pair-list" analysis).**
   *Treating every reported pair as a yes/no decision, how well does the score
   separate related pairs from unrelated pairs?* This is the classic database
   search benchmark.

2. **Category classification (the `.tcat` / "leave-category-out" analysis).**
   *If I use a domain's single best hit to predict its category, how often is
   the prediction right, and how often does the method confidently assign a
   category when it should have stayed silent?* This is a model of how the
   tool is actually used for annotation.

Sections 4-6 cover the pairwise analysis; sections 7-9 cover classification.


## 4. Pairwise truth standards

For the pairwise analysis, each non-ignored ordered pair `(q, t)` is a TP or
an FP. The four standards differ only in where they draw the TP/FP line:

| Standard        | TP (should match)   | FP (should not match)   | Ignored (beyond self/unknown) |
|-----------------|---------------------|-------------------------|-------------------------------|
| `fold`          | same fold           | different fold          | —                             |
| `superfamily`   | same superfamily    | different superfamily   | —                             |
| `family`        | same family         | different family        | —                             |
| `superfamilyx`  | same superfamily    | **different fold**      | **same fold, different superfamily** |

The first three are the obvious ones. `superfamilyx` deserves explanation.

**Why `superfamilyx` exists.** Suppose two domains share a fold but belong to
different superfamilies. Are they "unrelated"? Honestly, we don't know: they
might be distant homologs whose relationship SCOP simply hasn't committed to,
or they might be analogous (same shape, no shared ancestor). Counting such a
hit as a false positive would punish a method for finding something that may
well be real; counting it as a true positive would reward it for something
that may be noise. `superfamilyx` sidesteps the dilemma: a same-superfamily
hit is a TP, a *different-fold* hit is a clear FP, and the ambiguous
same-fold/different-superfamily middle ground is **ignored**. This gives a
cleaner read on superfamily-level sensitivity without the gray zone polluting
the error count.

A note on "ordered pairs". Because an all-vs-all file generally contains both
`q→t` and `t→q` (often with slightly different scores, since most aligners are
not perfectly symmetric), `homval` treats every reported line as its own pair.
The denominators in section 6 are counted the same way (ordered pairs,
excluding self).


## 5. The confusion-matrix vocabulary (gentle version)

All the metrics below are built from four counts. Imagine we pick a **score
threshold** and *accept* every hit at least as strong as the threshold:

- **TP** — accepted pairs that are truly related.
- **FP** — accepted pairs that are not related (errors we let through).
- **FN** — related pairs we *failed* to accept (missed; either not reported,
  or scored too weakly).
- **TN** — unrelated pairs we correctly did not accept.

As we lower the threshold (accept more), TP and FP both go up while FN and TN
go down. Every metric in this document is some ratio of these four numbers,
evaluated as the threshold is swept from strict to lenient.

Two denominators we precompute once from the lookup (they do not depend on any
algorithm):

- **N_possible_tp** — the total number of related (ordered) pairs that *could*
  be found under the truth standard. This is the denominator for "what
  fraction did we recover?"
- **ndom** — the number of domains (queries). Used for the error axis below.


## 6. Pairwise metrics, with a worked example

Let us define each metric, note its synonyms, and compute it on a tiny
example. Suppose, under the `superfamily` standard, there are
**N_possible_tp = 6** truly-related pairs in total and **ndom = 5** domains. We
sweep the hits best-score-first and, at some threshold, we have accepted
`cum_tp = 3` true and `cum_fp = 2` false pairs.

**Sensitivity = Recall = Coverage = True Positive Rate (TPR).**
The fraction of all findable true pairs that we have recovered:

```
coverage = cum_tp / N_possible_tp = 3 / 6 = 0.50
```

These four words mean the *same* number in this analysis. We mostly say
**coverage**.

**Precision (Positive Predictive Value).**
Of the pairs we accepted, the fraction that were correct:

```
precision = cum_tp / (cum_tp + cum_fp) = 3 / (3 + 2) = 0.60
```

**Error = Errors Per Query (EPQ).**
The average number of false hits a user would see per query domain at this
threshold:

```
error = cum_fp / ndom = 2 / 5 = 0.40
```

This is **not** a rate between 0 and 1; it can exceed 1 (more than one false
hit per query, on average). It answers the practical question "if I lower my
cutoff this far, how much junk do I wade through per query?"

**Specificity** = TN / (TN + FP), and the **False Positive Rate (FPR)** =
FP / (FP + TN) = 1 − specificity. We *can* compute these, but see the warning
below about why we do not feature them.

### Why we plot coverage-vs-error, and are wary of ROC

A natural instinct is to plot the ROC curve: FPR on the x-axis, TPR on the y.
The problem is **class imbalance**. In an all-vs-all search the number of
unrelated pairs (the FPR denominator, `N_possible_fp`) is enormous compared to
the number of related pairs. So even thousands of false positives produce a
*tiny* FPR, the ROC curve hugs the left-hand axis, and very different methods
look almost identical and almost perfect. ROC flatters everyone.

`homval` therefore makes **coverage-vs-error (CVE)** the primary plot:

- x-axis = **error** = cum_fp / ndom (errors per query), and
- y-axis = **coverage** = cum_tp / N_possible_tp.

Because the denominator of the x-axis is the number of *queries* (thousands),
not the number of *unrelated pairs* (tens of millions), a handful of false
positives moves the curve visibly. The low-error region (say, fewer than 0.1
or 1 false hits per query) is exactly where users operate, so that is where we
zoom in.

(ROC and precision-recall are still produced for completeness — we hide
nothing — but CVE is what we read.)

### E-value calibration: FPEPQ vs E-value

When scores are E-values, we can ask *how well-calibrated* they are.
FPEPQ = cum_fp / ndom is the observed mean number of false positives per
query at a given score cutoff.  An ideal E-value satisfies E = FPEPQ: the
E-value *is* the expected number of false positives.

`plot_epq.py` plots FPEPQ (y) against E-value (x) on log-log axes with a
dotted y=x reference line.  A curve close to the diagonal means the E-values
are well-calibrated; a curve above the diagonal means the method is
conservative (fewer FPs than predicted); below means it is anti-conservative.

If scores are P-values instead of E-values, pass `--dbsize N` to convert via
E = N * P.

### Summarising a CVE curve with one number: Sum3

Comparing whole curves by eye is hard, so we summarise each curve with a
single score that rewards high coverage *at low error*. We read the coverage
off the curve at three error levels and combine them with decreasing weight:

```
SEPQ0.1  = coverage at error <= 0.1   (very strict: <1 false hit per 10 queries)
SEPQ1    = coverage at error <= 1     (1 false hit per query)
SEPQ10   = coverage at error <= 10

Sum3 = 2 * SEPQ0.1  +  1.5 * SEPQ1  +  1 * SEPQ10
```

Higher Sum3 is better. The strict, low-error point carries the most weight
because that is the regime that matters in practice. (If a threshold is never
reached, the final coverage value is used in its place.)

### SFFP: sensitivity above first false positive

A complementary single-number summary that does not depend on threshold
weights:

```
SFFP = (number of TPs with score strictly better than the first FP)
       / N_possible_tp
```

The EDF rows are sorted best-first. SFFP records the fraction of all
possible TPs that appear above (strictly better score than) the first false
positive. If no FPs exist at all, all observed TPs qualify. Higher SFFP is
better; it measures how deeply the method can go before the first mistake.

SFFP applies only to the pairwise analysis (.edf), not to the
leave-category-out analysis (.tcat).


## 7. The second question: category classification

The pairwise analysis treats every pair independently. But in real use you
often have *one query* and ask: "what is it? what is its superfamily / fold?"
A simple model: **take the single best hit above some threshold, and copy its
category onto the query.** The leave-category-out analysis measures how good
that model is.

The trick is that we must not let the query find *itself* or near-trivial
relatives, or the test is too easy and leaks the answer. So we **hold out**
(remove from the searchable targets) part of the query's own annotation, and
we do this in two controlled experiments.

Throughout, the **query set** is restricted to domains for which the test is
well-posed (defined per variant below).

There are two variants, which are exact analogues of each other one level
apart:

- **topsf** — can we recover the **superfamily** of a query?
- **topfold** — can we recover the **fold** of a query?

We walk through `topsf` in detail; `topfold` is identical with
family→superfamily and superfamily→fold.


## 8. `topsf` step by step

**Query set.** Only domains whose superfamily contains **at least two
families** are queries. Why: in the positive control below we remove the
query's own family; if that were the *only* family in the superfamily, there
would be nothing left to find and the query would be unfairly doomed. Requiring
≥2 families guarantees a correct answer is *reachable*. (Note: "reachable" does
not mean "found" — the algorithm still has to score it well.)

We track exactly **two** best-scores per query `q`:

- **S_tp** = best score to a target in the **same superfamily but a different
  family** (the reachable correct answer).
- **S_xf** = best score to a target in a **different superfamily** (a wrong
  answer).

A key simplification: S_xf is the *same* quantity in both experiments below,
because removing the query's family and removing the query's whole superfamily
both leave every different-superfamily target untouched. So two numbers per
query are all we need.

### Experiment A — positive control (family held out)

We remove from the target pool every domain in the query's **own family**. A
correct same-superfamily answer is still present (it lives in a *different*
family of the same superfamily). At a given threshold T, look at the query's
**top hit** (the better of S_tp and S_xf):

1. **No hit reaches T** (both S_tp and S_xf are below T or absent) → **FN**.
   The method stayed silent when it could have answered.
2. **Top hit is same-superfamily** (S_tp ≥ T and S_tp ≥ S_xf) → **TP**.
   Correct category recovered. (By construction it is a different family, since
   the family was held out.)
3. **Top hit is different-superfamily** (S_xf ≥ T and S_xf > S_tp) → **FP**.
   A confident *wrong* answer: a wrong-superfamily hit outranked the correct
   one.

These three cases are exhaustive, so the positive control partitions the whole
query set: **TP + FP + FN = |Q|**.

### Experiment B — negative control (superfamily held out)

We remove from the target pool every domain in the query's **entire
superfamily**. Now a correct answer is *impossible*; every remaining target is
a different superfamily. At threshold T:

1. **There is a top hit** (S_xf ≥ T) → **FP**. The method confidently assigned
   a superfamily when the honest answer was "I don't know" — necessarily wrong.
2. **No hit reaches T** (S_xf < T or absent) → **TN**. Correctly silent.

So the negative control partitions the query set: **FP + TN = |Q|**.

### So there are two ways to produce a false positive

Yes — and this is deliberate:

- in the **positive** control, by ranking a wrong-superfamily hit above the
  reachable correct one (case A3), and
- in the **negative** control, by answering at all (case B1).

These are correlated: whenever A3 happens (S_xf ≥ T and S_xf > S_tp), B1 also
happens for that query (S_xf ≥ T). The converse is not true — a query can be a
clean TP in A while still being an FP in B (when S_tp ≥ S_xf ≥ T).


## 9. `topsf` metrics

**Coverage (sensitivity / TPR).** Fraction of answerable queries we classified
correctly:

```
coverage = TP / |Q|
```

Because the positive control partitions `|Q|` into TP/FP/FN, a case-A3
misclassification *lowers coverage* (it is not a TP) — exactly as it should.

**Error.** The error axis uses the **negative-control** false positives only:

```
error = FP_negative / |Q|         (this equals the negative-control FPR)
```

Interpretation: the rate at which the method makes a confident, necessarily
wrong call when the correct answer has been removed — a clean "false discovery
on a decoy" rate, bounded in [0, 1]. Under this definition the positive-control
case-A3 false positives still **reduce coverage** (they are not TPs) but are not
plotted on the error axis. Counting them on the error axis as well would
partially double-count, because A3 implies B1 for the same query (whenever a
wrong-superfamily hit outranks the correct one in the positive control, that
same wrong hit is also above threshold in the negative control). Keeping the
error axis as the negative-control FPR therefore keeps both axes clean,
bounded, and free of double counting.

**Specificity** = TN / (TN + FP_negative) = 1 − error. High specificity means
the method knows when to stay silent.

### Summarising with Top3

Exactly parallel to Sum3, but the x-axis is the negative-control FPR (the error
defined above), and we read coverage at three low error levels:

```
TEPQ0.001 = coverage at error >= 0.001
TEPQ0.01  = coverage at error >= 0.01
TEPQ0.1   = coverage at error >= 0.1

Top3 = 2 * TEPQ0.001  +  1.5 * TEPQ0.01  +  1 * TEPQ0.1
```

Higher Top3 is better; unreached thresholds back-fill with the final coverage.

**Why the Top3 thresholds (0.001, 0.01, 0.1) differ from the Sum3 thresholds
(0.1, 1, 10).** This is deliberate, and it follows directly from the very
different class imbalance in the two analyses:

- In the **pairwise** analysis the error axis is *errors per query*
  (cum_fp / ndom). The pool of unrelated pairs is enormous, so it is entirely
  normal — and still useful — to operate where a query accumulates on the order
  of 0.1, 1, even 10 false hits; those are the practically relevant cutoffs, so
  Sum3 samples there.
- In the **classification** analysis the error axis is a *bounded false
  positive rate* over queries (FP_negative / |Q|), where each query contributes
  at most one decision. A method that fires a wrong superfamily on 10% of decoy
  queries is already poor, so the practically interesting regime is much
  stricter — fractions of a percent up to ~10% — and Top3 samples 0.001, 0.01,
  0.1 accordingly.

In both cases the three thresholds span the range a practitioner would
typically consider reasonable for that quantity. Because the two summaries live
on different x-axes (an unbounded errors-per-query for Sum3, a bounded FPR for
Top3), their numeric values must **not** be compared across the two analyses.


## 10. `topfold` (the analogue)

Identical machinery, shifted one level up:

- **Query set**: domains whose **fold** contains at least two superfamilies.
- **S_tp** = best hit to a **same-fold, different-superfamily** target.
- **S_xf** = best hit to a **different-fold** target.
- **Experiment A (positive, superfamily held out)**: top hit same-fold → TP;
  top hit different-fold → FP; no hit → FN.
- **Experiment B (negative, fold held out)**: any top hit → FP; none → TN.
- Metrics and Top3 exactly as in section 9.


## 11. Glossary of synonyms

Different communities use different words for the same quantities. In this
project:

| Term we use | Same as | Definition (pairwise) |
|-------------|---------|------------------------|
| coverage    | sensitivity, recall, true positive rate (TPR) | cum_tp / N_possible_tp |
| precision   | positive predictive value (PPV) | cum_tp / (cum_tp + cum_fp) |
| error       | errors per query (EPQ) | cum_fp / ndom |
| specificity | true negative rate | TN / (TN + FP) |
| FPR         | fall-out, 1 − specificity | cum_fp / N_possible_fp |

In the classification (topcat) analysis, **coverage** = TP / |Q| and **error**
= negative-control FP / |Q| (= the negative-control FPR). The word "coverage"
is used in both analyses for the y-axis; the word "error" is used for the
x-axis in both, but it is *errors per query* in the pairwise analysis and a
*bounded false-positive rate* in the classification analysis. Keep the two
straight.


## 12. One-page summary

| | Pairwise (`.edf`) | Classification (`.tcat`) |
|---|---|---|
| Unit | one reported pair | one query domain |
| Question | are related pairs scored above unrelated ones? | does the top hit predict the right category? |
| Truth standards | fold, superfamily, family, superfamilyx | topsf, topfold |
| y-axis (coverage) | cum_tp / N_possible_tp | TP / \|Q\| |
| x-axis (error) | cum_fp / ndom (errors per query) | negative-control FP / \|Q\| (FPR) |
| Summary stat | Sum3 (error 0.1, 1, 10); SFFP | Top3 (FPR 0.001, 0.01, 0.1) |
| Also produced | precision-recall, ROC | — |
