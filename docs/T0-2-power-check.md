# T0-2 Phase 0 — Power check

**Scope statements that apply to every T0-2 output:**

1. **Control-side refresh with a frozen 2016 numerator.** Case allele counts come
   from the 2016 surveillance cohorts and are not being updated. This is not a
   full update of the 2016 analysis.
2. **ExAC is nested inside gnomAD v4.** Not two independent samples. Nothing here
   may be framed as two measurements disagreeing.
3. **Population-level estimates, not clinical guidance for any individual.**

**This document is not the reproduction gate.** It is an independent Python
recomputation from the 2016 supplementary tables, used to size the task before
installing R or pulling any new data. The reproduction gate (Phase A) runs the
original R and checks it against numbers printed in the paper.

Source: `github.com/ericminikel/prnp_penetrance` @ `da681e0` (2016-10-17).
Script: `t02_power_check.py`.

---

## 1. The interval method — checked, not inherited

The plan previously asserted in one line that the 2016 interval was "Wilson on
doubled allele counts," and built a falsifiable prediction on top of that. That
assertion was unverified. It has now been read from source
(`src/generate_figures.r:120-129`):

```r
penetrance_confint = function (ac_case, n_case, ac_control, n_control, ...) {
  case_confint    = binom.confint(x=ac_case,    n=2*n_case,    method='wilson')
  control_confint = binom.confint(x=ac_control, n=2*n_control, method='wilson')
  lower_bound   = penetrance(case_confint$lower, control_confint$upper, ...)
  best_estimate = penetrance(case_confint$mean,  control_confint$mean,  ...)
  upper_bound   = penetrance(case_confint$upper, control_confint$lower, ...)
}
```

**It is not case-only.** Both cohorts propagate, so control-side uncertainty is
genuinely in the published interval and a larger control cohort can genuinely
tighten it. The concern that gnomAD could not move the published interval no
matter what does not apply.

**But it is not the delta method either.** This is the Kirov et al. 2014
opposite-corner construction: the interval ends are formed by pairing each
cohort's CI bound against the other's opposite bound. On the log scale that makes
total width **additive** in the two component widths:

```
corner:  W ≈ W_case + W_control
delta:   W ≈ sqrt(W_case² + W_control²)
```

The corner method is always the wider of the two. A prediction table computed
under the delta formula would therefore have been predicting tightening under a
different variance model than the intervals are actually built with. **All
numbers below are computed under the repo's own corner method**, for both epochs,
so the two epochs are compared like with like.

## 2. Two structural facts that reshape the task

Reading the forest-plot construction (`generate_figures.r:134-161`) turned up two
things not visible from the paper abstract or the repo README.

### 2a. There is no E200K row. The four Mendelians are pooled, at zero controls.

```r
mendelians = c("P102L","A117V","D178N","E200K")
forest_data$ac_case[6:7]  = sum(ac_case["TOTAL",mendelians])   # 1034
forest_data$ac_control[6] = sum(ac_exac$ac[ac_exac$variant %in% mendelians]) # 0
```

P102L, A117V, D178N and E200K are **absent from `table_s03` entirely** (verified:
zero matching rows), so the sum is a true zero, not a missing value. They appear
in the forest plot only as a single pooled row.

This retires the E200K question we were both arguing. I claimed E200K would be
among the largest tighteners because it is control-dominated; that reasoning is
correct as arithmetic but E200K has no separately estimated interval in the 2016
analysis for gnomAD to tighten. The question was not well posed against the actual
paper, and neither answer was right.

### 2b. The 2016 controls are tiny, and already ancestry-matched.

| forest row | control cohort | ac_control | n_control |
|---|---|---|---|
| M232R | ExAC **JPT** | 5 | **663** |
| V180I | ExAC **JPT** | 2 | **663** |
| V210I | ExAC **TSI** | 2 | **4,795** |
| Mendelians (pooled) | ExAC all | 0 | 60,706 |

The single-variant rows do not use the 60,706-person cohort. They use one
1000-Genomes-style subpopulation each — 663 Japanese individuals for the two
Japanese-enriched variants, 4,795 Tuscans for the Italian-enriched one.

Two consequences:

- **The ancestry confounder is largely pre-handled by the original design.** The
  concern that a moved point estimate could be pure cohort composition is real,
  but the 2016 estimator already matches ancestry per variant. The refresh must
  preserve that matching — which makes it a design constraint, not a post-hoc
  correction.
- **"12× more controls" is the wrong framing per row.** The relevant multiplier
  is not gnomAD-total over ExAC-total. It is gnomAD's *East Asian* count over
  663, and gnomAD's *Southern European* count over 4,795. gnomAD v4 is UK
  Biobank-dominated and therefore European-skewed, so these two multipliers
  differ from each other and from 12×. Both are measured in Phase C, not assumed.

## 3. Recomputed 2016 intervals and the tightening ceiling

| forest row | penetrance (95% CI) | 1/x_case | 1/x_ctrl | dominated by | **ceiling** |
|---|---|---|---|---|---|
| M232R | 0.12% [0.04%, 0.34%] | 0.0149 | 0.2000 | control | **78.2%** |
| V180I | 0.96% [0.23%, 3.97%] | 0.0045 | 0.5000 | control | **91.1%** |
| V210I | 7.78% [1.85%, 32.71%] | 0.0058 | 0.5000 | control | **90.0%** |
| Mendelians (pooled) | 100% [29.4%, 100%] | 0.0010 | ∞ | control | n/a (censored) |

"Ceiling" is the maximum achievable reduction in log-interval-width as control
uncertainty goes to zero — the case-only interval. It is a hard limit set by the
frozen 2016 numerator, and it does **not** depend on guessing how big gnomAD's
matched subpopulations are.

**The premise holds, strongly.** All three single-variant rows are
control-dominated by a factor of 13–110×, with 78–91% of their interval width
removable in principle. This is not a marginal call.

Minimum control multiplier needed to realise a given tightening:

| forest row | 25% | 50% | 75% |
|---|---|---|---|
| M232R | 2.3× | 8.2× | 636× |
| V180I | 2.0× | 5.5× | 36.7× |
| V210I | 2.1× | 5.7× | 41.4× |

A 2× matched-ancestry cohort already buys 25%. Even ExAC JPT's 663 individuals
should be beaten by more than that in gnomAD v4.

## 4. The pooled Mendelian row behaves qualitatively differently

Its upper end is clamp-censored at 100% (`pmin(1, ...)`, `generate_figures.r:117`)
because `af_control = 0` makes the ratio infinite. Log-width is undefined, so the
row cannot "tighten." What more zero-allele controls do instead is push the
**lower** bound up:

| n_control (ac still 0) | interval |
|---|---|
| 60,706 (ExAC, 2016) | [29.4%, 100%] |
| 150,000 | [72.7%, 100%] |
| 303,530 | [100%, 100%] |

If gnomAD v4 still showed zero Mendelian alleles, this row would collapse to a
degenerate [100%, 100%] — the 2016 headline inverted. It almost certainly will
not: at 730,947 exomes some E200K and D178N carriers are expected. **So this row's
result is driven entirely by whether and how the zero resolves, which is a
different scientific question from interval tightening** and is reported
separately, never pooled into a "how much did intervals tighten" summary.

## 5. External cross-check against the paper

`generate_figures.r:288-289` annotates a sentence from the paper: the
Mendelians_ExAC interval contains published E200K survival-analysis penetrance
estimates of ~60–90%.

- Computed here, independently in Python: **[29.4%, 100.0%]** → contains 60–90%: **yes**.

This is a genuine external check — a paper-quoted claim recovered from the
supplementary tables without running the authors' R. It is one check, not the
reproduction gate, and it is weak (a wide interval trivially contains a
sub-range). Phase A does the real thing against printed point estimates.

## 6. Kill/resize gate — decided

> **Full scope** if ≥2 of the 4 frozen rows show projected log-width tightening
> **≥25%** at the *ancestry-matched* control multiplier actually observed in
> gnomAD v4.1.1 (measured in Phase C, step 14).
> **Reduced scope** otherwise: a short note reporting only the Mendelian
> zero-resolution, with the tightening claim dropped.

On present evidence the full-scope branch is expected to trigger: 3 of 4 rows
clear 25% at a matched multiplier of ~2.1×. The gate is retained anyway because
the matched multipliers are not yet measured, and the East Asian one is the
plausible failure case.
