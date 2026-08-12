# T0-4 — Pre-registration: baseline before any literature extraction

**Committed before any anchor extraction, endpoint harmonisation, or joint
fitting.** Verify with `git log` that this commit precedes any commit containing
newly extracted dose–response data.

Its purpose is to fix the comparator. Any later claim that "the range narrowed"
must be measured against the numbers below, not against a baseline reconstructed
afterwards.

Script: `t04_identifiability.py`. Model output, not evidence.

---

## 1. Baseline, frozen

Profile likelihood on `x_crit` over the 7 finite anchors in
`prion_model.ANCHORS`, with `β` solved analytically to hit each `x_crit` and the
unknown nuisances (`a`, `b`, `n`, `t_tox_frac`) minimised out by L-BFGS-B from 12
starts.

| quantity | value |
|---|---|
| profile minimum | SSR **0.3590** at **25% residual PrP** (75% knockdown) |
| **95% profile interval** | **7–38% residual PrP** |
| **as required knockdown** | **62% – 92%** |
| currently reported in `README` / `REPRODUCTION.md` | 64.3% – 91.9% |
| between-study σ | **0.299** log units (σ² = 0.0896), ≈ **35% multiplicative** |
| noise floor (constant fit to the five ~50% anchors) | SSR 0.3583 |
| null (one constant ratio, all 7 anchors) | SSR 2.7017 |

**Pre-registered comparator: 62–92% required knockdown, σ = 0.299.**

### Noise model, stated as the choice it is

σ² is estimated from disagreement among the five anchors at the *same* nominal
PrP level (~50%): survival ratios **1.52, 1.64, 2.00, 2.70, 3.00**. Residuals are
treated as i.i.d. log-normal. This attributes all within-cluster spread to
between-study noise rather than to real differences in intervention. If some of
that spread is real (different modalities genuinely differing at equal PrP
lowering), σ is overestimated and the true interval is narrower — but then the
anchors are not exchangeable and pooling them was wrong to begin with. Either way
the conclusion in §2 holds.

## 2. Falsified: the premise T0-4 was built on

T0-4 argued that the **Tga20 overexpression arm** (8–10× PrP) would constrain
`r(x)` from the opposite side of the curve and tighten `x_crit`.

| anchor set | required knockdown | min SSR |
|---|---|---|
| all anchors | **62% – 92%** | 0.3590 |
| **Tga20 removed entirely** | **60% – 95%** | 0.3587 |
| wild-type reference removed | 62% – 92% | 0.3590 |

**Deleting the overexpression arm widens the interval by ~3 percentage points.**
The premise does not survive. Recorded as falsified rather than quietly dropped.

The general form: **the binding constraint is σ, not the number of anchors.**
Adding points along the PrP-level axis cannot narrow `x_crit` while studies at the
same PrP level disagree by 35%.

### Fit quality, stated fairly

The model buys 2.34 SSR over a constant across all 7 anchors — but essentially
all of that is from fitting wild-type (x=1) and Tga20 (x=8), which a constant
cannot reach. **On the ~50% cluster it does no better than a constant** (0.3590 vs
0.3583), and that cluster is the only part of the data that could discriminate
`x_crit` from its neighbours.

## 3. Pre-registered predictions

| id | prediction | fails if |
|---|---|---|
| **Q1** | The **strain axis** (RML vs ME7 vs 22L at fixed host PrP) will not narrow the interval by more than 5 percentage points, because each strain introduces its own `β`, `b`, `a` — adding parameters at least as fast as constraint. | it narrows by >5 pp |
| **Q2** | Adding any further **survival-ratio** anchors at 40–60% residual PrP will not narrow the interval by more than 5 pp, unless they also reduce σ. | it narrows by >5 pp without σ falling |
| **Q3** | **Endpoint harmonisation** is the single largest reducible component of σ. Endpoint definitions (first clinical sign / confirmed diagnosis / terminal cull) differ by weeks in a ~150-day course, i.e. 10–20%, against a total scatter of 35%. | harmonising leaves σ ≥ 0.25 |
| **Q4** | Reaching **±5 pp on required knockdown** requires **σ ≤ 0.10** (scatter ≤ 11%), a 3× reduction from current. | achieved at σ > 0.15 |

Q4 is arithmetic from the frozen profile and cannot fail; it is recorded as the
**target**, not as a risky prediction. Q1–Q3 are risky and are assessed with
misses reported as prominently as hits.

### σ-to-interval map, frozen

| σ | scatter | required knockdown |
|---|---|---|
| 0.299 (current) | 35% | 62–92% |
| 0.250 | 28% | 63–89% |
| 0.200 | 22% | 66–86% |
| 0.150 | 16% | 68–82% |
| **0.100** | **11%** | **70–80%** |
| 0.050 | 5% | 73–77% |

## 4. Rules fixed in advance for any extraction

1. **Every anchor verified against a primary source** before entering a reported
   fit — PrP level as a fraction of wild-type with the assay named, incubation
   time, n, strain, route, mouse background, **and stated endpoint**. Unverified
   anchors are **dropped, not estimated**. The current `ANCHORS` table is
   unverified and the Tga20 entry is marked "approx"; the §2 conclusion is stated
   to be robust to its value only because a single point at x=8 cannot overcome
   σ=0.30, and that reasoning is what must hold, not the number.
2. **Endpoint exclusion rule, fixed now:** a study whose endpoint definition
   cannot be determined from its text is **excluded**. This is committed before
   seeing which direction exclusions push the estimate.
3. **Stratify by strain** where counts allow; report strain-stratified alongside
   pooled.
4. **Report σ explicitly** with every fit. A narrowed interval accompanied by an
   unchanged σ is an error, not a result.
5. **Null outcome is a first-class deliverable.** If σ cannot be brought below
   ~0.15, the write-up is *"the published dose–response cannot constrain `x_crit`
   beyond 62–92%; the limiting factor is 35% between-study scatter at fixed PrP
   level, not the number of studies; here is the titration experiment that would
   settle it"* — specified, powered, costed. Written at full length, not as a
   footnote.
6. **Literature check first.** T0-2's process failure was discovering an existing
   publication of the same numbers *after* completing the analysis. Before any
   fitting, search for existing joint fits of the PrP dose–response and for any
   published `x_crit` estimate.

## 5. Known limitations of this baseline

- `n` is treated as continuous (2–12) although it is a nucleus size. The model's
  formulae are smooth in `n`; restricting to integers would quantise the profile.
- σ² is estimated from 5 anchors (4 df) and is itself noisy.
- The profile interval uses a χ²(1) threshold scaled by σ², which assumes the
  log-residual noise model in §1.
- Anchors are literature values not yet verified (rule 1).
