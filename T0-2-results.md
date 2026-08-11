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

E200K's crude point estimate is **61.4%**, at the lower edge of the published
survival-analysis range (~60–90%) that the 2016 paper noted its interval
contained. **This overlap should not be read as independent methods converging.**
Both corrections identifiable here — control-cohort composition (§5) and
age-related depletion of manifested carriers (§5b) — push the estimate *below*
61.4%, so proximity to the survival range is more plausibly coincidence than
consilience. P102L and D178N remain clamped at a 100% point estimate — their
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

> **The specific figure 51.9% is not defensible and is withdrawn as a point
> estimate.** The standardization rests on three non-NFE alleles. `amr` alone
> contributes **42%** of the standardized allele frequency from **AC = 2**.
> Perturbing those counts by one allele:
>
> | perturbation | standardized penetrance |
> |---|---|
> | `amr` 2 → 1 | 65.7% |
> | observed | 51.9% |
> | `amr` 2 → 3 | 42.9% |
> | `amr` 1, `fin` 0 | 77.1% |
> | `amr` 3, `fin` 2 | 39.1% |
>
> **Report: the composition effect moves E200K downward from its crude 61.4%,
> plausibly into the 40–65% region, with a span of roughly 39–77% attributable to
> Poisson noise on three alleles alone.** The direction is well determined —
> v4's NFE over-representation understates an ExAC-composed control frequency —
> but the magnitude is not estimable from three non-NFE alleles, and calling it
> "anti-conservative" understated that. **The crude estimate and interval are the
> calibrated ones; the standardized figure is a direction with a wide range.**

## 5b. Onset-age check (pre-registered in `T0-2-prediction-prereg.md` §5)

The pre-registration committed to checking whether gnomAD v4 exposes age at
per-variant resolution and reporting the finding either way. **It does** —
`age_distribution` is available on the single-variant query (though not on the
gene-level variant list). Answered here rather than declined.

| variant | carriers with an age value | ages (5-year bins) |
|---|---|---|
| E200K | 5 of 13 | 40–45 ×1, 55–60 ×1, 60–65 ×1, 65–70 ×2 |
| P102L | 2 of 2 | 40–45 ×1, 50–55 ×1 |
| D178N | 1 of 1 | 55–60 ×1 |

**Three of the five aged E200K carriers are at or beyond the ~59–60 year mean
onset age for E200K, and were not affected at ascertainment.**

### Direction of the bias — derived, and it is not the one I was handed

The estimator needs `P(A)` = the allele frequency in the *general population,
including people who will later die of prion disease* (the paper states this
assumption explicitly in Methods). So counting pre-onset carriers in the
denominator is **correct**, not inflationary.

The real deviation runs the other way. gnomAD v4 is UK Biobank-heavy and UKB
recruits at **40–69**, which straddles E200K's onset. Prion disease kills within
about a year of onset, so carriers who manifested before recruitment cannot be
enrolled. The cohort is therefore **depleted** of exactly the highest-penetrance
carriers, making `af_control` too *low* and penetrance biased **upward**.

I was told this bias runs downward. I derive it as upward and am recording the
disagreement rather than adopting the framing. I hold this with genuine
uncertainty — it turns on whether the target denominator is the birth-cohort or
surviving-adult allele frequency, and the paper's own wording is what tips it.

### Consequence for the headline

Both identified corrections push the same way: composition (§5) and age
depletion both suggest **61.4% is an overestimate**, so the true value likely
sits *below* the published 60–90% survival range rather than at its lower edge.

That matters for how §1 should be read. "E200K lands at the bottom of the
published survival range" invites a consilience reading — independent methods
agreeing. **That reading is not supported.** Landing near 60% is what a control
set with residual pre-onset and depletion effects would produce, and the two
corrections we can identify both push below it. The agreement is more plausibly
coincidence than convergence.

Noting against my own inclination: I found it easy to read the survival-range
overlap as validation, and had to be pushed off it. That is the same pattern
recorded in the session's other reversals.

### Verdict

**Checked, and partially informative — not quantifiable.** Five ages on 13
carriers, in 5-year bins, cannot support a numerical correction. The
qualitative finding stands: unaffected E200K carriers past mean onset age exist
in the control set, which is direct evidence of incomplete penetrance, and is
consistent with a finite estimate without pinning its value. **No age correction
is applied to any number in this document.**

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
penetrance estimate here would be biased downward.

Concretely, for E200K. Clamping begins when
`af_control ≤ af_case × baseline = (571/20,920) × 2e-4 = 5.459e-6`, which at
AN = 1,461,878 is **AC = 7.98**:

| E200K control AC | penetrance |
|---|---|
| 13 (observed) | 61.4% |
| 10 | 79.8% |
| 9 | 88.7% |
| **8** | **99.8%** |
| **7** | **100% (clamped)** |

**The finite estimate survives losing 5 of the 13 alleles and disappears on
losing 6.** An earlier draft of this section claimed the threshold was AC ≤ 4,
which overstated the margin by roughly a factor of two and did so in the
direction that flatters the result. The margin is real but materially tighter
than that claim implied.

## 9. What this does and does not support

**Supported:** a first finite penetrance estimate for E200K against a large
control cohort (61.4% crude, 51.9% ancestry-standardized, 95% CI [33.1, 100]
crude); substantially raised lower bounds for P102L and D178N; and the finding
that the 2016 zeros were consistent with sampling, not evidence of a discrepancy.

**Not supported:** any reclassification; any statement about an individual's
risk; any claim that penetrance "changed" — the case numerator is frozen at 2016
and the control cohorts are nested, so what changed is precision, not the
underlying biology; and any comparison for M232R, V180I or V210I.
