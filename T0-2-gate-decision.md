# T0-2 — Gate decision: ancestry availability and scope

**Written after the availability check, before any PRNP allele counts were
retrieved from gnomAD.** Evaluates the gate defined in
`T0-2-prereg-amendment-1.md` §A3.

**Outcome: REDUCED SCOPE for M232R / V180I / V210I. Full scope for the unpooled
Mendelian variants, which are unaffected by the constraint.**

---

## 1. R cross-check: the primary estimand's baseline is now verified

Amendment 1 made interval width the gated quantity, and §6 of
`T0-2-reproduction.md` recorded that the three ExAC interval widths — the
denominators of every tightening percentage — were never printed in the paper and
so had never been checked against anything.

That gap is now closed by running the authors' own R.

**Environment:** conda-forge `r-base` 4.5.3 + `r-sqldf`; `binom` from CRAN.
`r-binom` is not available on conda-forge for osx-arm64.

**One documented patch.** R 4.5 rejects the 2016 script's
`fileEncoding='utf8'` (line 38) — modern R requires `'UTF-8'`. A copy,
`generate_figures_patched.r`, changes that one string and nothing else; the
`diff` is a single line and is reproduced in the session log. Without it the
script halts at line 38, before `forest_data` is built.

| forest row | R log-width | Python log-width | abs diff |
|---|---|---|---|
| M232R_ExAC | 2.1694870802 | 2.1694870802 | 4.4e-11 |
| V180I_ExAC | 2.8373072064 | 2.8373072064 | 3.3e-11 |
| V210I_ExAC | 2.8740439664 | 2.8740439664 | 3.2e-11 |
| Mendelians_ExAC | 1.2228465941 | 1.2228465941 | 2.9e-11 |

Agreement to ~1e-11, i.e. floating-point noise. The R run also confirms the
input derivation (`table_s01`/`s07`/`s08` → `forest_data`), not just the interval
function, since the script builds those rows itself.

**Every quantity the tightening analysis rests on is now verified** — point
estimates and two CIs against the paper, and all four interval widths against the
authors' code.

## 2. Ancestry resolution: checked, not assumed

Queried the gnomAD GraphQL API using a **non-PRNP gene (HBB)**, so the PRNP data
the gate governs remained unseen.

**gnomAD v4 exposes 9 groups, all continental:**
`afr, amr, asj, eas, fin, mid, nfe, remaining, sas`

**No subcontinental subdivisions.** Schema introspection confirms
`VariantPopulation` carries only `{id, ac, an}`; there is no subpopulation field.
The only ancestry-refinement field in the schema,
`local_ancestry_populations`, applies to admixed-American local ancestry in v4
genomes and does not provide Japanese or Southern European groupings.

**gnomAD v2.1.1 does have them, but only in the download VCF, not the API.**
Reading the v2.1.1 exomes VCF header directly confirms the INFO fields:

```
AN_eas_jpn  "Total number of alleles in samples of Japanese ancestry"
AN_eas_kor, AN_eas_oea
AN_nfe_seu  "Total number of alleles in samples of Southern European ancestry"
AN_nfe_nwe, AN_nfe_est, AN_nfe_swe, AN_nfe_bgr, AN_nfe_onf
```

So tier 1 does not exist; tier 2 exists. The remaining question is whether tier 2
cohorts are actually larger than 2016's.

## 3. They are not. Gate evaluation.

Subpopulation sample counts from gnomAD's own v2.1 release post:
`jpn` **76**, `kor` 1,909, `oea` 7,212 (EAS total 9,197);
`seu` **5,752**, `nwe` 21,111, `swe` 13,067, `onf` 15,499, `bgr` 1,335, `est` 121.

| row | 2016 matched cohort | tier-2 matched cohort | multiplier | projected tightening | gate |
|---|---|---|---|---|---|
| M232R | ExAC JPT n=663 | v2 `jpn` n=76 | **0.11×** | **−119%** (widens) | FAIL |
| V180I | ExAC JPT n=663 | v2 `jpn` n=76 | **0.11×** | **−115%** (widens) | FAIL |
| V210I | ExAC TSI n=4,795 | v2 `seu` n=5,752 | **1.20×** | **+7.0%** | FAIL |

**0 of 3 rows clear the ≥25% bar from a match-preserving cohort.**

The `jpn` result is counter-intuitive and worth stating precisely: the data did
not shrink. ExAC is a subset of gnomAD v2 exomes, so those 663 individuals are
still present — gnomAD's stricter ancestry inference simply assigns only 76 to
`jpn` and distributes the rest into `kor`/`oea`. This is a **labelling-granularity
mismatch, not a data-availability problem**, and it is not fixable by choosing a
different release.

Combining `jpn`+`kor` (n=1,985, 3.0×) was considered and **rejected**: Korean and
Japanese allele frequencies for M232R are not interchangeable, so that cohort
answers a different question. Substituting it would be precisely the estimand
substitution A1 exists to prevent.

### Decision for these three rows

Per A1 and the amended gate: **tightening is reported as uninterpretable.** No
gnomAD-based tightening figure is produced for M232R, V180I or V210I.

A v4 continental-group estimate (EAS-wide, NFE-wide) may still be reported, but
**only** as an explicitly different estimand — "lifetime risk in an East Asian
heterozygote" is not "lifetime risk in a Japanese heterozygote" — and is never
differenced against the 2016 value as though it were the same quantity.

## 4. The Mendelian row is unaffected, and is now the primary deliverable

The constraint in A1 applies to rows whose ancestry match degrades. **The
Mendelian row never had an ancestry match to degrade.** Its 2016 control cohort
is the full ExAC dataset:

```r
forest_data$n_control[6] = sum(exac_pop_summary$n_exac)   # 60,706 — all of ExAC
```

Unlike the three intermediate variants, this row was always a global,
unselected-population comparison. gnomAD v4's full cohort is **the same estimand
construct at ~12× the size**, so no substitution occurs and the tightening logic
applies without the A1 caveat.

Combined with the A2 unpooling rule — declared before any of this was known — the
reduced-scope deliverable is:

> **First per-variant penetrance estimates for P102L, A117V, D178N and E200K
> against a control cohort of ~730,000 exomes**, replacing a single pooled
> zero-count upper-bound statement from 2016.

Two things follow:

- A per-variant **zero** at n≈730k remains a substantive result (a tight lower
  bound on high penetrance), per A2. It is not a null.
- Ancestry composition still matters for **interpreting** the point estimate,
  since gnomAD v4's mix differs from ExAC's and E200K in particular has founder
  populations. The per-ancestry breakdown is reported alongside every estimate.
  This is a caveat on interpretation, not an estimand substitution, because the
  2016 quantity was itself a global mixture.

## 5. Revised prediction set

P1 and P2 are **void** — they were conditional on match preservation (A1), and no
row preserves its match. They are recorded as void, not as passed or failed, and
are excluded from the hit rate.

P3 (no row exceeds its ceiling) is retained for the Mendelian variants.

P4 is retained and assessed per variant per A2: for each of P102L, A117V, D178N,
E200K, does the control count resolve off zero? Both outcomes are informative.

**Net honest summary:** the pre-registration did its job by voiding the primary
estimand rather than letting a diluted cohort produce a narrow interval around
the wrong quantity. What survives is the part with the most scientific value, and
it survives because A2 unpooled it before the data was seen.
