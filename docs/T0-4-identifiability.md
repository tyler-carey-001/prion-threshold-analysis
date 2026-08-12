# Why the required-knockdown range is 65–90%, and what would narrow it

> Model output, not evidence. The kinetic model is a hypothesis generator. Anchor
> values are literature figures not yet verified against primary sources — see
> `T0-4-prereg.md` rule 1. Nothing here is clinical guidance.

**Summary.** The single therapeutically decision-relevant number this project
produces is the PrP knockdown required to push prion replication below
self-sustaining: **~65–90%**. That range is wide, and the obvious move was to
narrow it by fitting the published dose–response jointly rather than to one
anchor — using in particular the Tga20 overexpression arm, which sits at 8–10×
PrP and so constrains the growth-rate curve from the opposite side.

Running the identifiability check first shows two things:

1. **The published range is real.** A profile likelihood gives **62–92% required
   knockdown**, reproducing the reported 64.3–91.9%. The range was not an
   artifact of which anchor was chosen or of the nuisance sweep.
2. **The joint fit will not narrow it.** Deleting the Tga20 arm entirely moves
   the interval by ~3 percentage points. The binding constraint is not how many
   studies there are; it is that studies at the **same** PrP level disagree by
   **35%**.

Adding data along the PrP-level axis cannot fix a scatter problem. What would
narrow the range is reducing the scatter — and that requires **σ ≤ 0.10** for a
±5-point answer, a 3× reduction from today.

---

## 1. Setup

Under nucleated polymerisation the prion growth rate at substrate level `x` is
the dominant eigenvalue of the linearised 2×2 system, and the threshold where it
crosses zero is analytic:

```
x_crit = [ a(a + b(2n-1)) + b² n(n-1) ] / (b · β)
```

Required knockdown is `1 − x_crit/x0`. Survival is modelled two-phase after
Sandberg 2011: a replication phase scaling as `1/r(x)`, plus a toxic phase taken
as a fixed fraction `t_tox_frac` of untreated survival.

The published anchors are survival ratios at known residual PrP:

| residual PrP | survival ratio | source |
|---|---|---|
| 0.00 | ∞ (never) | `Prnp-/-`, Büeler 1993 |
| 0.49 | 2.70 | divalent siRNA, chronic pre-sx |
| 0.50 | 3.00 | ASO, chronic early (upper) |
| 0.50 | 2.00 | `Prnp+/-` heterozygote |
| 0.50 | 1.52 | base editing, prophylactic |
| 0.49 | 1.64 | di-siRNA, single dose at onset |
| 1.00 | 1.00 | wild-type reference |
| 8.00 | 0.42 | Tga20 overexpressor (approx) |

**Method.** Profile the sum of squared residuals in log survival ratio, holding
`x_crit` on a grid and minimising out every genuinely unknown nuisance
(`a`, `b`, `n`, `t_tox_frac`) at each point with L-BFGS-B from 12 starts. `β` is
solved analytically to place `x_crit` exactly. Interval from a χ²(1) threshold
scaled by σ².

## 2. The noise floor is the whole story

Five anchors sit at essentially the same PrP level (49–50%). Their survival
ratios are **1.52, 1.64, 2.00, 2.70, 3.00** — a two-fold spread at a fixed value
of the independent variable.

```
SSR of a single constant fitted to those five   :  0.3583
σ = 0.299 log units  ≈  35% multiplicative scatter
```

That is the floor. No model that treats these five as exchangeable observations
of one quantity can do better than 0.3583 on them.

**The profile minimum over all seven anchors is 0.3590.**

The model does buy 2.34 SSR against a constant fitted to all seven — but that
gain comes almost entirely from reaching wild-type (x=1) and Tga20 (x=8), which
a constant cannot. **On the 50% cluster it does no better than a constant**, and
that cluster is the only region where the data could discriminate one `x_crit`
from a neighbouring one.

## 3. The profile, and the interval

| | |
|---|---|
| minimum | SSR 0.3590 at **25% residual PrP** = **75% knockdown** |
| 95% interval | **7–38% residual** = **62–92% required knockdown** |
| previously reported | 64.3–91.9% |

The profile is **not** flat — it rises to 1.4 by 40% residual and the fit fails
entirely by 50%. But the width is large, and it is set by σ.

**This validates the headline number on better grounds than it previously had.**
The 65–90% figure had been justified as "the spread across published anchors and
across the `t_tox_frac` sweep," which is a sensitivity range, not an inference.
It turns out to coincide with the likelihood width. That is worth reporting, and
it was not guaranteed.

## 4. The overexpression arm carries no leverage

| anchor set | required knockdown | min SSR |
|---|---|---|
| all anchors | **62 – 92%** | 0.3590 |
| **Tga20 removed** | **60 – 95%** | 0.3587 |
| wild-type removed | 62 – 92% | 0.3590 |

Removing the entire overexpression arm — the one data point on the far side of
the curve, the point that motivated this whole task — **widens the interval by
about 3 percentage points**.

The reason is structural. Tga20 is a single observation carrying σ = 0.30 of
noise. One point at 35% uncertainty cannot pin a curve whose other seven points
are scattered by the same amount, however favourably it is placed.

**This generalises**, and it is the useful part: *any* addition along the
PrP-level axis inherits the same limit. More studies at 40–60% residual, or a
strain axis at fixed host PrP, add observations carrying the same 35% scatter —
and the strain axis adds free parameters (`β`, `b`, `a` per strain) at least as
fast as it adds constraint. Predictions Q1 and Q2 in `T0-4-prereg.md` record this
as testable rather than assumed.

## 5. What would actually work

Holding the anchor values fixed and shrinking only the scatter:

| σ | scatter | required knockdown | width |
|---|---|---|---|
| **0.299 (current)** | 35% | **62 – 92%** | 30 pts |
| 0.250 | 28% | 63 – 89% | 26 pts |
| 0.200 | 22% | 66 – 86% | 20 pts |
| 0.150 | 16% | 68 – 82% | 14 pts |
| **0.100** | **11%** | **70 – 80%** | **10 pts** |
| 0.050 | 5% | 73 – 77% | 4 pts |

**To reach ±5 points on required knockdown, between-study scatter must fall from
35% to ~11%.**

Three routes attack σ rather than adding points, and they are worth attempting in
this order:

1. **Endpoint harmonisation.** Studies differ on whether "incubation period"
   means first clinical sign, confirmed diagnosis, or terminal cull. In a ~150-day
   course those differ by weeks — **10–20%, against a total scatter of 35%**. This
   is the largest single reducible component and requires only careful reading.
2. **Time-to-plateau instead of survival.** Sandberg's two-phase result means
   prion titre rises to a plateau before the toxic phase begins. Fitting `r(x)`
   to plateau timing **bypasses the toxicity layer entirely** — removing
   `t_tox_frac` as a nuisance and removing the toxic phase's contribution to
   scatter. Highest ceiling, conditional on titre-vs-time curves actually being
   published.
3. **Incubation-time dispersion as a second observable.** As `x → x_crit`,
   `r → 0` and incubation-time variance diverges. Dispersion growing sharply
   between 50% and 40% residual localises the threshold in a way means cannot.
   This is genuinely orthogonal information, and it is currently being discarded
   by fitting means to data that report distributions.

## 6. What this means for the field

Ionis reopened PrProfile in March 2026 for a third, higher dosing regimen. The
model says the target is somewhere in **62–92% lowering**, and this analysis says
**the published mouse literature cannot say where in that range**, because
studies at nominally identical PrP levels disagree by a factor of two on survival.

That is not a call for more of the same experiments. It is an argument for the
one experiment that would settle it — a **within-study depth titration**, where
the scatter that dominates here is controlled by construction:

> Infect wild-type mice with RML. Titrate to stable residual PrP of ~60%, ~40%,
> ~25%, ~15%, chronic dosing from before inoculation, **one laboratory, one
> strain, one route, one endpoint definition, individual incubation times
> reported**. Measure survival and, at matched timepoints, brain prion seeding
> activity by RT-QuIC.

The design already appears in `RESEARCH_PLAN.md` §6. What this analysis adds is
**why it is necessary rather than merely desirable**: the cross-study route to the
same number is blocked, and it is blocked by a quantity (σ = 0.35) that no amount
of additional published data will reduce.

## 7. Status and honesty notes

- This is a **negative result about method**, established before doing the
  extraction work it evaluates. The check cost minutes; the extraction it
  redirects would have cost days.
- **T0-4's stated premise was falsified**, not quietly dropped. The task
  specification has been revised accordingly.
- The baseline is frozen in `T0-4-prereg.md` so any later claim of narrowing is
  measured against a committed number.
- σ is estimated from five anchors (4 df) and is itself uncertain. The noise model
  attributes all within-cluster spread to between-study noise; if some is real
  intervention difference, the anchors are not exchangeable and pooling was wrong
  to begin with. Both readings support §4's conclusion.
- Anchor values remain **unverified against primary sources**. The conclusions
  here are about the *structure* of the inference, which is why they can be stated
  ahead of that verification — but no fitted number should be reported until it is
  done.
