# Step 5 — Power check findings (committed before the held-out value is fetched)

**One-line result:** the toxicity comparison can *robustly falsify* standing
total-load toxicity, but it can only *weakly and asymmetrically* separate
flux-driven from neuronal-load-driven toxicity — because the discriminating
quantity (recovery latency) is largely non-identifiable from the training data.
This is committed **before** the Mallucci 2007 recovery timecourse is looked up.

Reproduce: `python step5_power_check.py`.

## What was profiled

Training pins: `beta` from the ANCHORS survival dose-response; `kappa` by
root-finding to put untreated terminal at 150 d (now correctly done in the
two-compartment model with each mode's own toxicity form — see the bug note
below); and the emergent untreated onset, required to sit in the
natural-history band for the early behavioural deficit. Everything else about
the toxicity layer is what training must pin for the latency prediction to be
sharp.

## Result 1 — the κ/ρ degeneracy makes flux latency nearly free

Profiling ρ over 0.01–0.5 (a 50× range), with κ refit to hold terminal at 150:

| quantity | behaviour across the 50× ρ range |
|---|---|
| untreated onset day | **76–81 d, essentially flat** (51–54% of course) |
| κ/ρ | varies only **2.5×** |
| κ, ρ individually | vary ~**50×** |
| flux 50%-recovery latency | **6.3 → 96.1 d** (15× swing) |

The untreated rise is quasi-steady, `D ≈ (κ/ρ)·flux`, so onset fixes the
**combination κ/ρ**, not ρ. But the flux recovery latency after knockdown is a
pure decay at rate ρ (`latency ≈ ln2/ρ`), set by ρ *alone*. So the one number
the held-out test turns on is the one number training cannot pin. No feature of
the untreated course breaks this — the frank-onset→terminal interval is set by
the growth rate and the threshold *ratio*, also independent of ρ. Only a
perturbation-recovery experiment (i.e. Mallucci itself) pins ρ.

**Committed flux latency prediction interval:** [6.3, 96.1] d under a defensible
±15 d onset anchor; still [6.3, 49.3] d even under an implausibly tight ±3 d
anchor. Flux is therefore **nearly unfalsifiable on latency**: almost any
observation is consistent with flux for *some* ρ.

## Result 2 — load_neuronal overlaps flux; load_total is dead

**Committed load_neuronal latency prediction interval:** [30.2, 124.0] d
(survivors; `a` free, `beta` refit per `a`). It **overlaps** the flux interval.

**load_total:** across the whole ρ scan it can **never** satisfy the reversal
criterion — D tracks total burden, which rises after neuron-specific knockdown,
so D does not fall. This is Mallucci's observation, and it is robust to the free
parameters. load_total is falsified independent of any fit.

## The discrimination structure (pre-registered here)

Overlap of the two committed intervals = **[30, 96] d**. Therefore, *before*
seeing the observed value, the pre-registered reading is:

- observed latency **< 30 d** → consistent with flux, **inconsistent with
  load_neuronal** → evidence against neuronal-load;
- observed latency **> 96 d** → inconsistent with flux, consistent with
  load_neuronal → evidence against flux;
- observed latency **30–96 d** → consistent with both → **no discrimination**.

**Asymmetric power, stated up front:** because flux's interval is so wide, this
test can plausibly *reject load_neuronal* (a fast recovery) but can almost never
*reject flux*. A "flux is consistent" outcome is weak evidence — flux is nearly
unfalsifiable on this axis, so it must not be reported as "flux wins."

## Rescue-vs-death is not a clean second axis

Flux survived late (70%-of-course) knockdown in 25/25 profiled cases;
load_neuronal in 10/12 (the two deaths were slow-clearance `a≈0.01–0.02` with
small ρ). So "does late knockdown rescue at all" is mostly yes for both and does
not cleanly separate them either.

## Bug caught and fixed during this step

The first version calibrated κ with the single-compartment `simulate`, whose
`_rhs` hardcodes the *flux* toxicity form and ignores `toxicity_mode`. Applied
to a load mode, that fit κ against the wrong dynamics (flux ~ β·x·y vs load ~ z,
orders of magnitude apart), pushing D to ~2.1 and killing the animal at day ~25
— long before the day-105 knockdown. The reported `inf` latencies were instant
deaths, not slow recovery. Fixed by refitting κ in the two-compartment model
with the correct mode (`refit_kappa` now uses `survival_time_2c`). The flux
result was unaffected (1c and 2c agree for flux by the reduction invariant).

## Consequence for the rest of T0-1

Running the expensive Step 8 fit and then fetching Mallucci 2007 to show "flux
matches" would be **circular**: flux can match nearly any latency by choice of
ρ. The scientifically defensible T0-1 outcomes are:

1. **Falsify load_total** — robust, needs no fit.
2. **Report flux vs load_neuronal as non-identifiable on latency**, with the
   committed intervals above, unless the observed value lands in a discriminating
   zone (<30 d or >96 d) — and even then, only load_neuronal can be cleanly
   rejected.
3. Name what *would* settle it: an independent measurement of ρ (D's repair
   rate), which the survival+onset training data cannot provide.

The held-out value is fetched next (Step 6) and the pre-registered logic above
is applied — **after** this file is committed.
