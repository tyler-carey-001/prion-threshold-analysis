# T0-2 — Pre-registration

**Committed before any gnomAD data is retrieved.** Verify with `git log` that
this commit precedes the results commit. Companion to `T0-2-power-check.md`,
which contains the arithmetic this document freezes.

---

## 1. Frozen variant list

The analysis covers exactly these four rows and no others. The list is taken
as-is from the 2016 forest plot (`generate_figures.r:135-161`), restricted to
rows whose control cohort is ExAC. Selection was made in 2016 on prior
pathogenicity reports, which is a defensible basis independent of anything
visible in gnomAD.

| # | row | control cohort (2016) | ac_case | n_case | ac_ctrl | n_ctrl |
|---|---|---|---|---|---|---|
| 1 | M232R | ExAC JPT | 67 | 1,533 | 5 | 663 |
| 2 | V180I | ExAC JPT | 222 | 1,533 | 2 | 663 |
| 3 | V210I | ExAC TSI | 171 | 1,054 | 2 | 4,795 |
| 4 | Mendelians pooled (P102L, A117V, D178N, E200K) | ExAC all | 1,034 | 10,460 | 0 | 60,706 |

**No variant may be added to this list because it looks interesting in gnomAD
v4.** Regression to the mean has a direction here: variants drew attention in
2016 partly because their ExAC counts were surprisingly high, and a larger sample
pulls those back down, which *raises* estimated penetrance and would read as "the
2016 paper overstated its case." That artifact is manufactured precisely by
v4-driven selection. The list above is the guard against it.

The 23andMe forest rows are excluded: a gnomAD refresh does not touch them.

## 2. Primary deliverable: interval precision, not point movement

**Primary estimand — change in 95% interval width** (log scale), per row.

This is the quantity the 2016 paper's headline was about, it is nearly guaranteed
to move for control-dominated rows, and it is clean: it does not require the
allele frequency to be stable between cohorts, only that more controls were
sequenced.

**Secondary estimand — change in point estimate.** Reported, but labelled
explicitly as confounded, for reasons that are structural and known in advance:

```
shift = penetrance_now / penetrance_2016 = af_control_ExAC / af_control_gnomAD
```

- ExAC is **nested** inside gnomAD v4 and contributes roughly 8% of its exomes,
  so this ratio is driven by the newly-added samples versus ExAC, not by two
  independent measurements of the same quantity.
- The 2016 denominators are **2 to 5 alleles**. Poisson noise on a count of 2 is
  enormous; a large part of any observed shift is that noise resolving.
- Neither component is "penetrance changed."

Any movement claim is therefore reported **only** in its ancestry-matched form
(§4), and the unmatched version is not reported as a result at all.

## 3. What this pre-registration can and cannot fail on

The power model in `T0-2-power-check.md` is **arithmetic and cannot fail.** Given
allele counts, the corner-method interval width is determined. Predicting that
control-dominated rows tighten is not a risky prediction about the world.

**The risky content — the part that can actually fail — is allele-frequency
stability between cohorts.** The projections assume gnomAD's matched-ancestry
allele frequency for each variant equals ExAC's, so that `ac_control` scales with
`n_control`. If gnomAD's frequencies differ materially, realised tightening will
depart from projection, and for rows where frequency rose the interval could even
widen. A "hit" on the projection table must be read as *the frequencies were
stable*, not as *the power model was validated*.

### Pre-registered predictions

| # | prediction | fails if |
|---|---|---|
| P1 | M232R, V180I, V210I each realise ≥25% log-width tightening | any of the three realises <25% |
| P2 | Realised tightening for each of the three lands within ±15 percentage points of the projection at that row's measured matched multiplier | any row departs by >15 pp |
| P3 | No row exceeds its ceiling (78.2% / 91.1% / 90.0%) | any does — would indicate a coding error, since the ceiling is a hard limit given frozen case counts |
| P4 | The pooled Mendelian row resolves off zero, i.e. gnomAD shows ac_control > 0 for at least one of P102L, A117V, D178N, E200K in a matched cohort | it does not, in which case the row collapses to a degenerate [100%, 100%] and is reported as uninformative rather than as a finding |

Misses are reported with the same prominence as hits, including the hit rate.

## 4. Ancestry matching is a design constraint, not a correction

The 2016 estimator is already ancestry-matched per variant (JPT for the Japanese
variants, TSI for the Italian one). The refresh **preserves that matching**:
each row's gnomAD control cohort is the closest available genetic-ancestry group
to its 2016 cohort, and the matched multiplier for that row is measured, not
assumed to be 12×.

Where gnomAD's ancestry labels do not map cleanly onto the 2016 1000-Genomes-style
subpopulation labels (JPT and TSI are population-level; gnomAD v4 reports
continental groups such as EAS and NFE with finer subdivisions of varying
availability), the mapping used is stated explicitly per row, along with what was
available. Any row where a defensible match cannot be made is reported as such
and excluded from the movement claim, not silently matched to a continental group.

## 5. Freed-parameter treatment, applied with the structure in mind

Three levers exist: assumed prevalence (`baseline_risk = 2e-4`), case
ascertainment, and an onset-age correction on controls.

**For absolute penetrance:** full sweep over prevalence and ascertainment,
reported as a band on every point estimate. No point estimate without a CI, and
no CI without its parameter band.

**For the shift (secondary estimand):** prevalence and ascertainment **cancel
exactly** — prevalence is the same multiplier in both epochs, ascertainment acts
on `af_case`, which is identical in both by construction (frozen numerator).
Sweeping them against the shift would inflate the band on the one quantity that
is well determined. They are therefore not swept there, and this document is the
justification.

**The onset-age correction is not applied, and the shift is reported unbanded on
that lever.** Reasons, pre-registered:

- It is **not in the original analysis.** The 2016 code has no age term; it
  assumes controls are below onset age. Introducing one would be my addition, not
  a reproduction.
- A lever I design myself, applied to the headline quantity, with a functional
  form I choose after seeing the data, is exactly where the pinned-parameter
  failure would live. There is no published per-variant onset distribution for
  these four rows at the resolution the correction would need.
- The 2016 control counts are **2 to 5 alleles**. An age-stratified correction on
  a denominator of 2 is noise amplification, not a correction.

**Phase C must still check** whether gnomAD v4.1.1 exposes age at per-variant
resolution. If it does not — or if carrier counts are in the 5–50 range where age
data is uninformative — that finding is reported as the reason the correction
remains inapplicable, which is a stronger statement than declining it a priori.
If per-variant age data turns out to be both available and well-populated, the
correction is added **as a clearly-labelled sensitivity that is specified in
writing and committed before it is computed**, never as an adjustment to a
headline number.

## 6. Kill/resize gate

> **Full scope** if ≥2 of the 4 frozen rows show projected log-width tightening
> ≥25% at the ancestry-matched multiplier measured in Phase C.
> **Reduced scope** otherwise: a note reporting only the Mendelian
> zero-resolution, tightening claim dropped.
