# The experiment this repository argues for

> Model output, not evidence. This is a design argument derived from a
> lumped-parameter model and published anchors, not a protocol. Anyone running it
> should treat the numbers as starting points. Prion work is BSL-3 and there are
> confirmed occupational fatalities — see the biosafety note at the end.

**One sentence:** infect mice with RML, hold PrP at ~60 / 40 / 30 / 25% of normal
in separate arms, and measure **prion titre over time** as well as survival — in
one laboratory, one strain, one route, one endpoint definition, with individual
animal times reported.

That experiment answers a question the published literature structurally cannot,
and it is answerable with about ten animals per arm.

---

## 1. The question

**Does lowering PrP below ~40% of normal stop prion replication, or merely slow
it further?**

Two structures are compatible with everything published:

| | mechanism | prediction below 40% PrP |
|---|---|---|
| **Threshold** (nucleated polymerisation) | replication rate falls with substrate; below a critical concentration `x_crit` polymer creation cannot outrun clearance | **disease stops entirely** |
| **No threshold** (Sandberg two-phase) | exponential phase is not PrP-limited; the plateau phase sets onset, ∝1/[PrP] | **disease still occurs, ~3× slower** |

These are not fringe alternatives. The first is this repository's model. The
second follows directly from Sandberg et al. 2011 (*Nature* 470:540), whose
measured kinetics have the opposite structure to it.

The distinction decides the target product profile for every PrP-lowering
programme. Under a threshold, the goal is "clear a specific, knowable bar."
Without one, the goal is "lower as much as tolerable, forever."

## 2. Why the published data cannot answer it

Three analyses in this repository, each a negative result:

1. **`T0-4-structural-test.md` — the threshold sits outside the data.**
   `x_crit` calibrates to **33% residual PrP**. The lowest data point in the
   entire published dosage series is **49%**. The two structures agree within
   ~2× everywhere data exists and diverge **13-fold** below 40%. Both are
   validated exactly where they cannot be told apart.
2. **`T0-4-identifiability.md` — more published data will not help.** Studies at
   the *same* nominal PrP level report survival ratios of 1.52, 1.64, 2.00, 2.70
   and 3.00 — **35% between-study scatter**. Removing the Tga20 overexpression
   arm entirely moves the interval by ~3 points. The binding constraint is the
   scatter, not the number of studies.
3. **`T0-1-findings.md` — the toxicity question is separately blocked.**
   Flux-driven and neuronal-load-driven toxicity are *nested*: flux is the
   fast-clearance limit of neuronal load. No reversal-timing or knockdown-depth
   experiment separates them; only a direct clearance measurement does.

## 3. Design

### Arms

| arm | target residual PrP | purpose |
|---|---|---|
| untreated control | 100% | anchor the course; ~150 d for RML |
| A | ~60% | reproduce the known regime; sanity check |
| B | ~40% | first point where the structures separate (7.7× vs 2.1×) |
| C | ~30% | threshold model predicts **no disease** |
| D | ~25% | threshold model predicts **no disease**; still within lifespan |

**Chronic dosing from before inoculation**, so PrP is at target throughout.
Divalent siRNA now reaches ~17% residual with a single dose, so these depths are
achievable in a way they were not when the two-phase work was done.

### Why the ladder stops at 25% — a constraint that is easy to miss

Predicted incubation, if the no-threshold structure is correct:

| residual PrP | predicted incubation | vs ~730 d mouse lifespan |
|---|---|---|
| 60% | 225 d | fine |
| 40% | 319 d | fine |
| 30% | 412 d | fine |
| **25%** | **487 d** | **feasible** |
| 20% | 599 d | marginal |
| **15%** | **787 d** | **exceeds lifespan** |

**An arm whose predicted incubation approaches natural lifespan cannot
distinguish "very slow disease" from "no disease"** — the animals die of other
causes either way, and a negative result means nothing.

`RESEARCH_PLAN.md` §6 proposes a 15% arm. On these numbers that arm is
**uninterpretable**, and the 60% arm is nearly uninformative (both structures
predict 1.5–2.1×). The discriminating window is **40–25% residual PrP**. That is
the single most useful refinement this analysis offers to the design.

### Readouts

**Primary: prion titre versus time**, by RT-QuIC seeding activity or bioassay, at
matched timepoints — not survival alone.

This matters more than it sounds. The two structures disagree about *which phase*
lengthens:

- Threshold model: the **exponential** phase slows and eventually cannot start.
- Sandberg model: the exponential phase is unchanged; the **plateau** lengthens.

Survival confounds both phases into one number and therefore discards exactly the
information that discriminates. A titre timecourse separates them directly. This
is also the observable that would let a future meta-analysis bypass the toxicity
layer entirely.

**Secondary: survival**, with **individual animal times**, not just mean ± SD.
Near a threshold, incubation-time variance diverges as `r → 0`; dispersion
growing sharply between arms B and C is itself evidence of a threshold, and it is
information routinely discarded by reporting means.

**Tertiary, and it answers a different question:** in a subset, deplete neuronal
PrP mid-infection (the Mallucci 2003 design) and measure **neuronal** PrP-Sc
clearance against behavioural recovery. `T0-1-findings.md` shows this is the only
measurement that separates flux-driven from neuronal-load-driven toxicity, and
that it is decisive if neuronal PrP-Sc clears with a half-life longer than roughly
**8–37 days** (a model-derived band that should be rebuilt on the real assay).

### What must be held constant

The 35% scatter that blocks the meta-analytic route comes from cross-study
variation. Controlling it is most of the design's value:

- **one laboratory**, one mouse background, one prion strain (RML), one
  inoculation route and dose;
- **one endpoint definition**, stated explicitly — first clinical sign vs
  confirmed diagnosis vs terminal cull differ by weeks in a 150-day course, which
  is comparable to the effect being measured;
- **PrP level measured, not assumed** — assay stated, measured repeatedly through
  the course, since the whole design depends on residual PrP actually being
  stable at target.

### Size

For the **structural** question, this is cheap. If the threshold model is right,
arms C and D give 0 events by day 600; if it is wrong, essentially all animals
convert by ~410–490 d. Distinguishing 0/n from n/n needs only n ≈ 5 for
p < 0.01; **n ≈ 10–15 per arm** gives comfortable margin and allows for
intercurrent deaths over a 500-day study.

For **locating `x_crit` precisely** — as opposed to establishing whether it
exists — more animals and more arms are needed, and the required number depends
on incubation-time dispersion, which is exactly the quantity current papers
discard. Establishing existence first is the right order.

## 4. What each outcome means

| observation | conclusion |
|---|---|
| Arms C and D convert at ~410–490 d, exponential phase unchanged | **No threshold.** Nucleated-polymerisation framing is wrong; PrP lowering buys time proportionally and the target is "as deep as tolerable, indefinitely." |
| Arms C and D show no disease by day 600, with titre falling or flat | **Threshold exists**, and lies between arms B and C. Target product profile becomes a knowable bar. |
| Arm B (40%) converts near 319 d but C does not | Threshold between 40% and 30%; the range this repository could not narrow gets located directly. |
| Exponential-phase rate falls with PrP in any arm | Contradicts Sandberg's PrP-independence **below** the range they studied — which their own hedge, "over the range we studied," leaves open. |
| Dispersion grows sharply from B to C | Independent evidence of threshold behaviour, from data usually thrown away. |

**Every outcome is publishable and every outcome changes what the field should do
next.** That is the test for a well-posed experiment, and it is why this one is
worth the 18 months it would take.

## 5. Honest limitations

- The predicted incubation times come from a **lumped-parameter model fitted to
  unverified literature anchors**. The *ordering* and the *lifespan constraint*
  are robust; the specific day counts are not.
- The `x_crit` = 33% figure is itself model-derived and is the thing being
  tested. If the threshold exists but sits at 20%, arms C and D would convert and
  the result would be misread as "no threshold." A follow-on arm at 15% is
  therefore worth running **in parallel with a longer-lived mouse strain**, or
  the study should be read as bounding rather than locating.
- Achieving *stable* residual PrP at 25% for 500 days is the hard experimental
  problem, and nothing here addresses it.
- Sandberg's two-phase result is the basis for the no-threshold alternative, and
  it was established over 50–800% PrP. Extending its mechanism below 50% is an
  extrapolation on that side too — which is the point.

## 6. Biosafety — not a formality

Prions resist autoclaving, formalin and most disinfectants. **Émilie Jaumain**, an
INRAE technician, pricked her thumb through two gloves in 2010 handling
prion-infected humanised mouse brain and died of vCJD in 2019, aged 33 — the first
pathologically confirmed occupational prion transmission. A second French case led
five institutions to impose a moratorium in 2021, and a Spanish prion researcher
has also died of CJD.

This experiment belongs in a properly equipped institution with trained staff and
institutional approval. Nothing in this repository should be read as encouraging
anyone to improvise it.
