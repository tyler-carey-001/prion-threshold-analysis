# T0-2 Results — PRNP penetrance against gnomAD v4.1.1

> **These are population-level estimates. They are not clinical guidance for any
> individual, and no variant is reclassified here.** Reclassification is a
> clinical-genetics judgement requiring segregation and functional data, and
> belongs to expert panels.

> **This is a control-side refresh.** Case allele counts are frozen at the 2016
> surveillance cohorts (10,460 sequenced cases) and are *not* updated. It is not
> a full re-run of the 2016 analysis.

> **ExAC ⊂ gnomAD v2 ⊂ (largely) gnomAD v4.** These are nested samples, not
> independent replications, and nothing below is framed as two measurements
> disagreeing.

Retrieved 2026-08-11. gnomAD **v4.1.1**, exomes, GRCh38. Raw responses cached in
`data/gnomad_prnp_v4_raw.json`, `data/gnomad_prnp_v2_raw.json`,
`data/gnomad_prnp_v4_coverage.json`. Analysis: `t02_penetrance.py`.

v4.1.1 was confirmed from its release notes to share v4.1's sample set — it
changes constraint metrics, LOFTEE flags and annotations only, so allele counts
are the v4.1 counts.

---

## 1. Headline

**The zeros resolved for three of four variants, and E200K now has a finite
penetrance estimate for the first time.**

| variant | 2016 (ExAC, n=60,706) | gnomAD v4.1.1 crude (n≈730,939) | v4 AC | status |
|---|---|---|---|---|
| **E200K** | 100% [15.9, 100] | **61.4% [33.1, 100]** | 13 | RESOLVED |
| **P102L** | 100% [5.9, 100] | 100% [**37.1**, 100] | 2 | RESOLVED |
| **D178N** | 100% [5.5, 100] | 100% [**45.1**, 100] | 1 | RESOLVED |
| **A117V** | 100% [0.7, 100] | 100% [8.6, 100] | 0 | zero persists |

E200K's crude point estimate of **61.4%** sits at the lower edge of the published
survival-analysis range (~60–90%) that the 2016 paper noted its interval
contained. P102L and D178N remain clamped at a 100% point estimate — their
control frequencies are still too low to pull the ratio below 1 — but their
**lower bounds rose roughly six- and eight-fold**, which is where the information
gain actually lands.

## 2. The 2016 zeros were never surprising — they were underpowered

The most important framing result, and it argues against reading this as a
correction to the 2016 paper.

If v4's allele frequencies are taken as the truth, how many alleles *should* ExAC
have seen in 2 × 60,706 = 121,412 chromosomes?

| variant | v4 AF | expected in ExAC | P(observe 0) | verdict |
|---|---|---|---|---|
| P102L | 1.37e-06 | 0.17 | 0.85 | unremarkable |
| A117V | 0 | 0.00 | 1.00 | unremarkable |
| D178N | 6.84e-07 | 0.08 | 0.92 | unremarkable |
| E200K | 8.89e-06 | 1.08 | 0.34 | unremarkable |

**Every 2016 zero is entirely consistent with the v4 frequencies at ExAC's sample
size.** Even E200K, the most common of the four, had an expectation of only ~1
allele in ExAC; seeing zero had probability 0.34. The 2016 analysis was not
wrong. It was underpowered for exactly the reason its authors stated, and the
refresh supplies the power rather than correcting an error.

## 3. Re-calling vs. new samples (A7)

The pre-registered check was whether any resolution is an artifact of
re-processing largely the same samples rather than new sequencing.

| variant | ExAC | v2.1.1 (AN=251,466) | v4.1.1 | v2 expected under v4 AF | P(obs ≤ observed) |
|---|---|---|---|---|---|
| P102L | 0 | 0 | 2 | 0.34 | 0.71 |
| A117V | 0 | 0 | 0 | 0.00 | 1.00 |
| D178N | 0 | 0 | 1 | 0.17 | 0.84 |
| E200K | 0 | **1** | 13 | 2.24 | 0.35 |

**All three call sets are mutually consistent with a single underlying frequency
per variant.** No count is anomalous given the sample size at which it was
observed, so **no re-calling explanation is required or supported.** The
progression 0 → 1 → 13 for E200K is sample-size scaling.

An earlier draft of this analysis labelled E200K's resolution "re-calling on
overlapping samples" purely because it appeared in v2. That was wrong: v2 is
~2× ExAC, so a v2 detection is just as likely to be a new sample. The label was
corrected before publication of these results and the consistency test above is
what actually discriminates.

## 4. Detection quality (A6)

Filtered records were queried explicitly, not PASS-only. The response contains
`AS_VQSR` (25 exome, 3 genome) and `AC0` (11 exome, 3 genome) records at other
PRNP sites, confirming the query surfaces filtered variants.

| variant | v4 filters | flags | exome coverage (mean) | ≥20× |
|---|---|---|---|---|
| P102L | `[]` PASS | none | 50.2× | 100% |
| A117V | *no record* | — | 50.9× | 100% |
| D178N | `[]` PASS | none | 34.9× | 100% |
| E200K | `[]` PASS | none | 34.8× | 100% |

All three detections are **PASS with no flags** — none is a filtered-only
detection.

**A117V's zero is a true absence, not a coverage gap.** There is no record at
c.350C>T in either call set. The surrounding sequence is densely represented —
c.340, c.342, c.344, c.345, c.351, c.357, c.359, c.362, c.365 all carry records,
*including* `AC0`- and `AS_VQSR`-filtered ones at c.345 and c.362 — which
demonstrates that zero-count and low-quality sites at this locus do get reported
when they exist. Coverage in the bracketing bins is ~50× mean with 100% of
samples above 20×, higher than at D178N or E200K.

## 5. Ancestry composition (A4)

gnomAD v4 exomes are **76.1% NFE** by allele number; ExAC's inferred composition
was **54.1% NFE**. The pooled control frequency is a composition-weighted
average, and the weights changed substantially between epochs.

Direct standardization to ExAC's composition, using weights frozen before the
pull (nfe 54.15%, sas 13.86%, amr 9.90%, afr 8.40%, eas 7.14%, fin 6.55%);
6.3% of v4's allele number sits in `asj`/`mid`/`remaining`, which have no 2016
counterpart and are excluded:

| variant | crude AF | standardized AF | crude penetrance | standardized |
|---|---|---|---|---|
| P102L | 1.368e-06 | 9.740e-07 | 100% | 100% |
| D178N | 6.841e-07 | 4.870e-07 | 100% | 100% |
| **E200K** | 8.893e-06 | **1.052e-05** | **61.4%** | **51.9%** |

Per-ancestry detail (v4 exomes):

```
P102L   nfe=2/1,111,886
D178N   nfe=1/1,112,012
E200K   nfe=10/1,112,012   amr=2/44,724   fin=1/53,404
```

E200K's frequency is ~5× higher in `amr` and ~2× higher in `fin` than in `nfe`.
Because v4 is far more NFE-weighted than ExAC was, the crude v4 frequency
*understates* what an ExAC-composed cohort would show — so **standardizing lowers
E200K's penetrance from 61.4% to 51.9%.** That ~10-point gap **is** the
composition effect, and it is reported rather than absorbed into either number.

> **Caveat, stated because it materially limits the standardized figure.** The
> standardized estimate converts a weighted allele frequency into an effective
> count at the full allele number and reuses the binomial interval machinery.
> That propagates sampling error in the total but **not** the extra variance from
> reweighting very small per-group counts (`amr` AC=2, `fin` AC=1). The
> standardized interval is therefore **anti-conservative** and should be read as
> indicating the direction and rough size of the composition effect, not as a
> calibrated interval. **The crude interval is the calibrated one.**

## 6. Pre-registered predictions — outcomes

| id | prediction | outcome |
|---|---|---|
| P1 | three single-variant rows tighten ≥25% | **VOID** — no ancestry-matched cohort exists (see `T0-2-gate-decision.md`) |
| P2 | realised tightening within ±15 pp of projection | **VOID** — same reason |
| P3 | no row exceeds its case-only ceiling | **HELD** — all four intervals remain wider than their ceilings (0.162–0.679 log-width) |
| P4 | Mendelian zero resolves for ≥1 variant | **HELD** — resolved for 3 of 4 |

**A5's carve-out fired as designed.** A117V's zero persists, giving a lower bound
of only **8.6%** — pre-declared *uninformative* because 33 case alleles cap the
bound regardless of control cohort size. It is reported here as uninformative and
is **not** presented alongside P102L's 37.1% or D178N's 45.1% as though
comparable. Had this not been fixed in advance, A117V's "zero persists at 730,000
exomes" would have been an easy and misleading headline.

## 7. Limitation inherited from the 2016 analysis (A8)

The 2016 M232R and V180I estimates used ExAC individuals assigned to a **JPT**
label, n = 663. gnomAD's own ancestry inference assigns only **76** individuals
to `jpn` in v2.1.1 — and ExAC is nested inside v2, so those ~587 individuals are
still present, classified otherwise.

**The published 2016 estimates for M232R and V180I therefore rest on a population
assignment that gnomAD itself no longer endorses.** This is a property of the
published analysis, independent of any refresh, and it was established without
retrieving any PRNP data.

Stated as a difference between two ancestry-inference methods, **not** as an error
in the 2016 work. Establishing which assignment is better would require
individual-level genetic data, which this project does not download.

This is also why M232R, V180I and V210I are absent from the tables above: no
match-preserving control cohort exists in any current gnomAD release, so their
interval changes are uninterpretable rather than merely uncertain.

## 8. Figure

`fig_t02_penetrance_shift.png` — 2016 ExAC (grey) vs gnomAD v4.1.1 (teal;
A117V in red as the persisting zero), crude estimates with 95% corner-method
intervals, shaded band showing the published E200K survival range.

**What would falsify the claim this figure depicts:** if the v4 allele counts for
P102L, D178N or E200K were substantially inflated by contamination, sample
duplication, or misannotation — or if a material fraction of gnomAD v4
contributors were ascertained for neurodegenerative disease, violating the
control assumption — the control frequencies would be overstated and every
penetrance estimate here would be biased downward. Concretely: if E200K's true
control allele count in an unselected cohort were ≤4 rather than 13, its
penetrance would return to the 100%-clamped regime and the headline result would
disappear.

## 9. What this does and does not support

**Supported:** a first finite penetrance estimate for E200K against a large
control cohort (61.4% crude, 51.9% ancestry-standardized, 95% CI [33.1, 100]
crude); substantially raised lower bounds for P102L and D178N; and the finding
that the 2016 zeros were consistent with sampling, not evidence of a discrepancy.

**Not supported:** any reclassification; any statement about an individual's
risk; any claim that penetrance "changed" — the case numerator is frozen at 2016
and the control cohorts are nested, so what changed is precision, not the
underlying biology; and any comparison for M232R, V180I or V210I.
