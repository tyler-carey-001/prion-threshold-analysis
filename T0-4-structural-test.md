# The threshold may not exist: a structural test of `x_crit`

> Model output, not evidence. Anchor values remain unverified against primary
> sources (`T0-4-prereg.md` rule 1). Not clinical guidance.

**Summary.** The literature check that gated T0-4 turned up a paper already cited
in this repo's own README whose measured kinetics have the **opposite structure**
to the model. The load-bearing consequence is not a model-selection verdict but an
extrapolation one:

> **`x_crit` sits at 33% residual PrP. The lowest data point in existence is 49%.
> The competing structures agree within ~2× everywhere data exists and diverge
> 13-fold below 40%, where nobody has measured.** The threshold is extrapolated
> beyond every observation — and so is its absence.

**The headline claim — "65–90% knockdown required" — is therefore a property of an
assumed model structure rather than a measured quantity.** That argument is
airtight and does not depend on any model comparison. §3 reports a comparison
anyway, but it is *supporting* evidence and weaker than it first appeared — see
the correction there.

`RESEARCH_PLAN.md` §6 named this in advance as the falsification condition:

> "If there is no inflection — if survival scales smoothly all the way down — the
> nucleated-polymerization framing is wrong and something else limits propagation."

---

## 1. The model is built the wrong way round

Sandberg et al. 2011 (*Nature* 470:540) measured prion titre over time in mice at
three PrP expression levels. Verified quotes (from Sandberg 2014, PMC4104459,
restating the result):

> "In phase 1, prions propagate exponentially, **not rate limited by cellular
> prion protein (PrP^C) concentration**, rapidly rising by ~10⁶-fold to reach a
> maximal prion titre, which is also independent of PrP^C concentration **over the
> range we studied**."

> "This is followed by a plateau phase (phase 2), **which determines time to
> clinical onset of disease, the duration of which is inversely proportional to
> PrP^C concentration**."

> "Remarkably, total PrP levels ... were essentially unchanged throughout phase 1
> in wild-type, Prnp+/o and Tg20 mice, a period during which titres rose in each
> case by ~10⁶-fold."

This repo's `relative_survival_analytic` assumes the reverse assignment:

| | where the PrP dependence lives | PrP-independent part |
|---|---|---|
| **this model** | replication phase, via `r(x)` | toxic phase (`t_tox_frac`) |
| **Sandberg (measured)** | plateau phase, ∝ 1/[PrP] | exponential phase |

Both end up with roughly a quarter of the course PrP-independent — the model's
default `t_tox_frac = 0.25`, and the fitted `p1 = 0.251` below. **They agree on
the fraction and disagree on which phase it is** — and that disagreement is
exactly what determines whether a threshold exists.

## 2. Sandberg's mechanism, written as a model

Phase 1 duration independent of `x`; phase 2 duration ∝ 1/`x`. Normalising so
untreated survival is 1:

```
survival_ratio(x) = p1 + (1 - p1) / x        [one free parameter]
```

**This has no threshold.** As `x → 0`, survival rises smoothly as 1/`x` and
diverges only at `x = 0` — which is consistent with `Prnp-/-` mice being
completely resistant. It reproduces the knockout without any critical
concentration.

## 3. Model selection

| model | k | SSR | AIC | AICc |
|---|---|---|---|---|
| **Sandberg, 1/x, no threshold** | **1** | 0.5484 | **−15.83** | **−15.03** |
| NPM (this repo) | 5 | 0.3590 | −10.79 | +49.21 |

ΔAIC = **5.0**, ΔAICc = **64.2**, both nominally favouring the threshold-free
model.

### Correction: this comparison is much weaker than ΔAICc = 64 suggests

An earlier version of this document led with that number. Decomposing it
(`t04_bootstrap.py`):

| term | free (k=1) | NPM (k=5) | difference |
|---|---|---|---|
| `n·ln(SSR/n)` — **the fit** | −17.83 | −20.79 | **−2.97** |
| `2k` — AIC penalty | 2.00 | 10.00 | +8.00 |
| `2k(k+1)/(n−k−1)` — AICc correction | 0.80 | 60.00 | **+59.20** |
| **total** | | | **+64.23** |

**92% of ΔAICc is an arithmetic penalty, not evidence.** At n=7, k=5 the
denominator `n−k−1` equals 1, so the correction term explodes independently of
how well anything fits. **The fit term itself favours the NPM**, by −2.97.

A 400-sample bootstrap over the seven anchors gives a stable sign (99.1% favour
the threshold-free model, 2.5–97.5 percentile +42 to +67) — but **that stability
is an artifact of the same deterministic penalty**, not a property of the data.

Two further corrections to the earlier version:

- The NPM sits at SSR 0.3590, essentially **at** the noise floor (0.3583, a
  constant fitted to the five ~50% anchors). The threshold-free form sits at
  0.5484, **above** it. On the cluster carrying most of the data, the NPM matches
  a constant and the threshold-free model does worse than one.
- The claim that the threshold-free form fits Tga20 *better* was **wrong**. It
  used the repo's old calibrated `β` rather than a joint refit. Refitted on the
  same footing the NPM fits Tga20 **exactly** (0.42 predicted, 0.42 observed)
  while the threshold-free form misses by −0.197 in log units.

**What the comparison actually supports:** seven anchors — five of them at a single
x-value — cannot support five parameters. That is a statement about this dataset's
resolving power. **It is not a statement that the threshold-free structure
describes the data better.** By fit alone, the NPM is ahead.

## 4. The threshold lies outside the data

| residual PrP | NPM | Sandberg | ratio | data exists? |
|---|---|---|---|---|
| 8.00 | 0.19 | 0.34 | 0.5 | **yes** |
| 1.00 | 1.00 | 1.00 | 1.0 | **yes** |
| 0.50 | 3.27 | 1.75 | 1.9 | **yes** |
| 0.40 | 7.74 | 2.12 | 3.6 | no |
| 0.35 | 30.8 | 2.39 | 12.9 | no |
| **0.334** | **divergent** | 2.49 | ∞ | no |
| 0.25 | divergent | 3.25 | ∞ | no |

**The calibrated `x_crit` is 0.334 residual PrP. The lowest anchor is 0.49.**

The two models are within ~2× of each other at every point where data exists, and
diverge 13-fold by 35% residual. Both are fitted on, and validated in, precisely
the region where they cannot be told apart.

Sandberg's own hedge is the same one — *"over the range we studied"* — and that
range was Prnp+/o (50%), wild-type, and Tg20 (8×). **Nobody has measured below
~50% residual PrP.** The threshold is an extrapolation, and so is its absence.

### The hedge cuts both ways, and this is the right framing

It would be easy to read Sandberg as *refuting* the threshold. It does not.

The NPM's own prediction is that phase-1 PrP-independence **must break down as
`x → x_crit`** — that is what a threshold means. Sandberg measured
PrP-independence over 50–800% and found it held. **A threshold lying below 50% is
entirely compatible with that observation.** The measurement was taken in the
region where the NPM also predicts near-independence.

So these are not two competing theories with a winner. They are **the same
unmeasured region seen from two sides**: one model extrapolating a threshold down
into it, the other extrapolating smooth 1/x behaviour down into it, and no
measurement in between. That is a stronger and more honest statement than either
"the threshold exists" or "the threshold is refuted."

## 5. What this does and does not establish

**Established:**

- The model's phase assignment contradicts the directly measured kinetics in a
  paper this repo already cites.
- A one-parameter model derived from those measured kinetics accounts for the
  same data adequately. Formal model selection nominally prefers it, but 92% of
  that margin is a small-sample arithmetic penalty and the fit term favours the
  NPM (§3) — so this is weak support, not a result.
- `x_crit` at 33% residual is extrapolated beyond every existing observation.
- Therefore **the 65–90% figure cannot be presented as a data-supported
  requirement.** It is what one assumed structure implies.

**Not established:**

- That there is no threshold. Sandberg's result covers 0.5–8× PrP; a threshold
  below 50% residual is untested, not excluded. `Prnp-/-` resistance proves
  *something* changes at low PrP.
- That the threshold-free structure is better supported. By fit it is worse; the
  model comparison turns on a parameter-count penalty (§3).
- That the NPM is wrong as biophysics. It is a well-founded model; the claim here
  is that **this dataset cannot support it over a simpler alternative**, which is
  a statement about the evidence, not the mechanism.
- Any revised knockdown target. The honest position is that the required depth is
  **unknown**, not that it is lower.

## 6. What would discriminate — and it sharpens the experiment

| residual PrP | NPM predicts | Sandberg predicts |
|---|---|---|
| 40% | 7.7× | 2.1× |
| 30% | divergent | 2.7× |
| 25% | divergent | 3.2× |
| 15% | divergent | 5.2× |

**A depth titration must go below 40% residual PrP to say anything at all.** Above
that the models are within a factor of ~2 and the 35% between-study scatter
(`T0-4-identifiability.md`) swamps the difference.

This refines the design in `RESEARCH_PLAN.md` §6, which already proposes 60 / 40 /
25 / 15% residual. That is the right ladder. What this analysis adds:

- **The 60% arm is nearly uninformative** for the structural question — both
  models predict ~1.5–1.7× there.
- **The 25% and 15% arms carry essentially all the discriminating power**, where
  the predictions differ by a factor of ∞ versus 3–5×.
- The readout should include **titre-vs-time**, not just survival, because the
  two models make different predictions about *which phase* lengthens — and that
  is a more direct test than survival, which confounds both phases.

## 7. Consequence for the project

The 65–90% range was this project's one therapeutically decision-relevant output.
Two analyses have now landed on it in sequence:

1. `T0-4-identifiability.md` — the range is the genuine likelihood width, and no
   amount of additional published data narrows it, because between-study scatter
   at fixed PrP is 35%.
2. This document — the range's *existence* depends on a model structure whose
   defining feature sits entirely outside the measured region, and the field's
   own kinetic data are equally compatible with a structure that has no threshold
   at all.

Together these say the same thing the toxicity work said in T0-1: **the question
cannot be settled from published data, and the specific experiment that would
settle it is identifiable.** That is a negative result, and it is the useful kind
— it redirects effort rather than adding a number to guess with.

The README's headline claim has been amended accordingly.
