# T0-4 — Joint fit of the PrP dose–response; tighten `x_crit`

> **STATUS: PREMISE FALSIFIED (see `T0-4-identifiability.md`).** The central
> argument below — that the Tga20 overexpression arm would tighten `x_crit` by
> constraining `r(x)` from the opposite side — was tested before extraction and
> does not hold. Removing Tga20 entirely changes the interval by ~3 percentage
> points. The binding constraint is **between-study scatter (σ ≈ 35% at fixed PrP
> level)**, not the number of anchors, so adding points along the PrP axis cannot
> narrow the range.
>
> **The task is not cancelled, it is redirected:** from *"fit jointly to tighten"*
> to *"reduce σ, or establish that it cannot be reduced."* The revised priority
> order is endpoint harmonisation → time-to-plateau → dispersion, with the strain
> axis demoted to a test of prediction Q1 rather than a source of constraint.
> Baseline and predictions are frozen in `T0-4-prereg.md`.
>
> The original text is kept below unedited, because the argument it makes is
> plausible and the record of it failing is the useful part.

---

**Replaces T0-3 (fibril structures) in the queue.** T0-3 is interesting and not
decision-relevant. This is.

## Why this and not more genetics

T0-2 produced a real result, but penetrance moves counselling and trial
*eligibility*, not treatment. The therapeutic content of this whole project sits
in one number: the **65–90% required-knockdown range**, which survived T0-1
untouched because it lives in the replication layer, not the toxicity layer that
T0-1 rebuilt.

That range is embarrassingly wide, and it is wide for an avoidable reason:
`β` was fitted to essentially **one anchor** — the survival ratio at ~50%
lowering. The literature contains a dose–response spanning roughly two orders of
magnitude in PrP level that nobody appears to have fitted jointly.

Ionis reopened PrProfile for a third, higher dosing regimen without knowing where
the bar is. So did everyone else. This is the number every PrP-lowering programme
currently guesses at.

## The anchors

Approximate, **every one to be verified against primary sources before fitting**
— no anchor enters the fit on the strength of this list. Per `CLAUDE.md` rule 4,
any value not confirmed from the paper is marked `[UNVERIFIED]` and excluded.

| genotype / intervention | PrP level | incubation effect |
|---|---|---|
| `Prnp-/-` | 0× | never develops disease |
| `Prnp+/-` | ~50% | ~2× |
| wild-type | 1× | 1× (reference) |
| divalent siRNA, chronic presymptomatic | 49% residual | 2.7× |
| PrP-lowering ASO, chronic early | ~50% | up to ~3× |
| AAV base editing (R37X), humanised | ~50% | +52% lifespan |
| **Tga20 overexpressor** | **8–10×** | **incubation cut to ~62 d** |

**The overexpression arm is the point.** Every therapeutic anchor sits between
0.5× and 1×, so they constrain one side of `r(x)` and are nearly collinear. Tga20
constrains the curve from the opposite side at 8–10×, which is where the leverage
on curvature — and therefore on `x_crit` — actually lives. Fitting the two arms
jointly should tighten `x_crit` substantially versus fitting a single point.

This is the depth-titration experiment §6 of `RESEARCH_PLAN.md` says to run,
approximated by meta-analysis across published genotypes and doses. Laptop-
tractable, uses everything already built.

## Method

1. **Verify every anchor** against its primary source: PrP level (as a fraction
   of wild-type, with the assay stated), incubation time, n, strain, route,
   mouse background. Anchors that cannot be verified are dropped, not estimated.
2. **Pre-register before fitting** (this is a hard gate, see below).
3. Express each anchor as `(x_i, T_i)` — relative PrP level, incubation time —
   with uncertainty on both. Fit `β` (and the growth-rate map `r(x)`) jointly by
   weighted least squares or a simple hierarchical model, propagating anchor
   uncertainty into `x_crit`.
4. Compare the joint-fit `x_crit` interval against the current single-anchor
   64.3–91.9%. **Report the interval, never a point.**
5. Leave-one-out: refit dropping each anchor in turn, especially Tga20. If
   `x_crit` depends materially on a single study, say so prominently.

## The methodological risk, and the pre-registration that guards it

**Cross-study pooling can manufacture a curve.** Strain (RML vs ME7 vs 263K),
inoculation route (i.c. vs i.p.), dose, and mouse background all shift incubation
time independently of PrP level. A fit that ignores this will produce a
confident, wrong `x_crit`, and it will look clean.

Pre-register, before any fitting:

- **Stratify by strain** where anchor counts allow; report strain-stratified
  `x_crit` alongside the pooled one.
- **Report between-study heterogeneity explicitly** (τ² or equivalent), not just
  the fitted parameter.
- **Pre-commit to the null as a first-class deliverable, not a fallback.** The
  Tga20 overexpression data is old, from multiple labs, with strain and
  background differences. **The prior probability that a joint fit yields a clean
  `x_crit` is not high**, and going in expecting otherwise is how the previous
  three tasks generated retractions.

  > If between-study heterogeneity dominates the PrP-level signal, the deliverable
  > is: **"the published dose–response cannot constrain the threshold, and here is
  > the titration experiment that would"** — with the experiment specified,
  > powered, and costed. That is a finding the field can act on, and it has the
  > same shape as T0-1's nesting result: a well-established impossibility that
  > redirects effort rather than a number that invites false confidence.

  State the heterogeneity threshold **numerically before fitting**. A null
  outcome under that threshold is written up at full length, not as a footnote.
- State which anchors are fitted and which are held out.

### Do the literature check first

T0-2's process failure was finding Minikel's April 2024 gnomAD v4 post *after*
completing the analysis it duplicated. **Before any fitting, search for existing
joint fits of the PrP dose–response and for any published `x_crit` estimate.**
If one exists, this task becomes a replication or an extension, and that has to
be known at the start.

This is the same discipline that made T0-2's A5 carve-out work: the honest
failure mode has to be declared while it is still cheap.

## Second deliverable — prevention-trial power at 61% penetrance

The bridge from T0-2 to what the field is actually blocked on.

Given E200K penetrance of ~61% (crude; and lower under the corrections in
`T0-2-results.md` §5 and §5b), lifetime risk in carriers, and a plausible
biomarker-triggered enrolment window, compute: how many carriers, followed for
how long, to detect a given reduction in conversion-to-disease?

Report across the penetrance range T0-2 actually supports, not at a point — the
whole argument of T0-2 is that the point estimate is soft. A prevention trial
sized on 100% penetrance is sized wrong, and that is the practical consequence of
the genetics work.

## Deliverables

- `T0-4-prereg.md` — anchors, stratification plan, heterogeneity threshold,
  null-outcome commitment. Committed before fitting.
- `T0-4-anchors.md` — every anchor with primary source, verified values, and
  anything dropped with the reason.
- `t04_dose_response.py`
- `T0-4-results.md` — joint-fit `x_crit` interval vs the current 64.3–91.9%,
  strain-stratified, leave-one-out, heterogeneity reported.
- `t04_prevention_power.py` + power curves across the penetrance range.

## Acceptance

- No anchor used without a verified primary source.
- Pre-registration commit precedes the fitting commit in `git log`.
- `x_crit` reported as an interval, with the heterogeneity statistic beside it.
- Leave-one-out reported, Tga20 specifically.
- Null outcome reported if the pre-registered heterogeneity threshold is crossed.
