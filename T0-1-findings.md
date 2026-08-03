# T0-1 findings — is the flux toxicity assumption load-bearing?

**Headline.** The two-compartment model **falsifies standing-total-load
toxicity** robustly. It **cannot cleanly separate flux-driven from
neuronal-load-driven toxicity** on recovery *latency* — the discriminating
quantity there is non-identifiable from the available training data. It *can*
separate them by whether Mallucci's toxicity/burden *dissociation* recurs at
partial knockdown depths — a parameter-robust readout — and that points to a
concrete two-readout successor experiment (behaviour + RT-QuIC seeding) tied to
`x_crit`.

This reframes the original T0-1 question. The task asked whether a load-driven
variant can reproduce Mallucci's reversal after fair refitting. The answer
depends critically on *which* load model: total-load cannot; neuronal-load can,
and is not distinguishable from flux by the observation the task proposed to use.

**Which observations were fitted vs held out.**
- *Fitted:* the untreated course (emergent onset ~76 d, terminal fixed at 150 d)
  and the `ANCHORS` 50%-lowering survival dose-response (`prion_model.py:272`).
- *Held out:* the Mallucci reversal and its recovery timecourse (Mallucci 2007,
  fetched only after the power-check prediction intervals were committed —
  `STEP5_FINDINGS.md`, `CITATIONS.md`).

---

## 1. load_total is falsified — robustly, without any fit

Under neuron-specific knockdown (`k_n≈0.9`, `k_g=0`) the untreated glial
compartment keeps converting, so total burden `z_n+z_g` rises. A toxicity term
driven by that total therefore cannot fall: across a full ρ scan with κ refit,
`load_total` never satisfies the reversal criterion (D-drop >50% while burden
rises). This is exactly Mallucci 2003 — "the propagation of nonneuronal PrPSc is
not pathogenic" (PMID 14593181, quoted in `CITATIONS.md`) — and it holds
independent of the free parameters. **The intuitive "toxicity tracks plaque
load" model is dead.**

## 2. flux vs neuronal-load: latency cannot separate them

The plan intended recovery *latency* to be the discriminator: flux collapses
within the drug ramp (days), neuronal load decays on the clearance timescale
`~1/a`. The Step 5 power check (committed before the held-out value was seen)
showed this fails, for a structural reason:

- The untreated course is quasi-steady, `D ≈ (κ/ρ)·flux`, so onset pins only the
  **combination κ/ρ**. Profiling ρ over 50× leaves onset flat (76–81 d) while
  the flux recovery latency swings **6→96 d**. Flux is therefore *nearly
  unfalsifiable on latency*.
- The held-out value (Mallucci 2007, PMID 17270731): after depletion at ~9 wpi,
  synaptic function and novel-object recognition recover at 9 wpi and burrowing
  by 10 wpi — a recovery latency of **~7–14 days** (untreated terminal ~13 wpi).
- flux reproduces ~7–14 d trivially, at ρ≈0.05–0.10 (the default ρ was 0.06).
- **neuronal-load can also reproduce it** — but only with fast polymer clearance
  (`a≳0.12`, i.e. PrP-Sc clearing in ≲8 d) plus high repair (ρ≈2). Crucially,
  the survival anchors all cluster at ~50% knockdown (depths −7, 0, 0.5, 0.51),
  so **training does not pin `a`**. My initial [30,124] d floor for neuronal-load
  was partly an artifact of capping `a` at 0.08; with `a` genuinely free (per the
  reviewer's rule that any parameter driving a mode's held-out prediction must be
  free), neuronal-load reaches the observed range.

So the observed fast recovery **disfavours** neuronal-load — it requires
PrP-Sc clearing within days, which sits at the edge of biological plausibility
(PrP-Sc is protease-resistant and accumulates over the disease) — but it does
**not** cleanly falsify it, because training cannot exclude that corner. Latency
is the wrong axis. Reporting "flux wins" here would be circular.

## 3. The discriminator is knockdown *depth* — but resolved by clause, and it
needs a second readout

`t01_depth_discriminator.py` / `fig_t01_depth_discriminator.png`. An earlier
draft reported a "~55% residual recovery threshold" for flux. That was
**wrong to call recovery**: it measured only whether D nets *any* drop 60 d
after knockdown — the weakest of the three pre-registered clauses, and not even
the >50% version. Resolved against the actual clauses (A: D falls >50% within
40 d; B: total load rising; C: survival >450 d), knockdown at day 90:

| depth (residual PrP) | flux | neuronal-load |
|---|---|---|
| clause A (>50% D-drop in 40 d) | met at ≤ **25%** | **never** met (a=0.02: clearance too slow for a 40 d window) |
| clause C (survival > 450 d) | met at ≤ **40%** | met at ≤ **32%** |
| full reversal (A∧B∧C) | ~**25%** residual (~75% KD) | not achieved on a 40 d window |

Two consequences the earlier draft obscured:

- **Survival rescue does not discriminate the modes.** Clause C turns on near
  x_crit (33% residual) for *both* — rescue requires driving replication
  subcritical regardless of the toxicity mechanism. A survival titration
  **locates x_crit**; it does not reveal the mechanism. (At 55% residual, flux
  satisfies *none* of the three clauses cleanly — D-drop only 16%, survival
  270 d — so "flux recovers at 55%" was conflating transient dysfunction relief
  with rescue. Corrected.)
- **The mechanism shows in the *dissociation*.** Does behaviour improve at a
  depth where replication is still supercritical (residual > x_crit, seeding
  still rising)? Max behavioural improvement while supercritical: **flux 38%,
  neuronal-load 2%.** Flux permits Mallucci's dissociation to recur at partial
  knockdown — the animal improves while total prion load is still climbing;
  neuronal-load cannot, because its D only falls once `z_n` clears, which needs
  subcritical replication.

**The successor experiment is underdetermined with one readout.** An observed
behavioural-recovery threshold at 30% residual is equally consistent with flux
at x_crit=15% and neuronal-load at x_crit=30%: two unknowns (recovery depth and
x_crit, which the original analysis pins only loosely to 10–35% residual), one
measurement. The fix — half-present in `RESEARCH_PLAN.md` §6 — is to **add an
RT-QuIC seeding readout**: the depth at which seeding activity starts declining
gives x_crit from *replication alone*, independent of any toxicity assumption;
the depth at which behaviour recovers gives the toxicity threshold.

- recovery at a **shallower** depth than seeding-decline → **flux**;
- recovery **coinciding** with seeding-decline → **neuronal-load**.

Equivalently, and as a yes/no in a single cohort: **at a partial knockdown that
leaves replication supercritical (seeding still rising), does behaviour
recover?** Yes → flux; no → neuronal-load. This replaces the non-actionable
"measure ρ" (ρ is a lumped phenomenological rate, not a wet-lab observable) and
ties the discriminator to x_crit — the number the whole project turns on.

---

## Deviations from the pre-registration (declared, per V7 substance-not-just-order)

- **The Latin-hypercube escape-hatch sweep was not run.** Its purpose was to
  search for load-model parameters reproducing Mallucci. That purpose is already
  served analytically: `load_total` provably cannot reverse (§1), and
  `load_neuronal` provably *can* (§2), so there is no escape region left to
  discover by sampling. Running it would spend hours confirming what the closed
  analysis shows. This is a deliberate deviation, recorded here rather than
  omitted silently.
- **The latency criterion was demoted, not applied as a pass/fail gate.** The
  power check (committed before the fetch) showed it lacks the power to gate on;
  applying it as if it discriminated would misrepresent the evidence. The
  committed prediction intervals stand; the conclusion drawn from them is the
  honest asymmetric one above.

## What the power check bought

Catching the κ/ρ non-identifiability *before* Step 8 meant the ~2–3 h multi-start
fit + sweep was never spent on an underpowered test, and the naive
"flux-matches-Mallucci-so-flux-wins" conclusion — which the original
single-compartment figures gesture toward — was avoided. The negative result on
latency is itself the substantive finding, alongside the robust death of
load_total and the depth-dependence discriminator that replaces it.

## Consequence for the repo's claims

Nothing here weakens the **replication-threshold** result (the ~65–90% knockdown
headline), which lives in the validated replication layer (`test_growth_rate.py`
confirms the determinant condition; `REPRODUCTION.md` confirms the calibration).
It sharpens the **toxicity** story: standing-load toxicity is falsified, flux is
supported but not proven against neuronal-load, and the experiment that would
settle it is a depth titration — reinforcing, not competing with, the project's
central `x_crit` question.
