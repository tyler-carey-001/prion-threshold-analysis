# T0-1 findings — is the flux toxicity assumption load-bearing?

**Headline.** The two-compartment model **falsifies standing-total-load
toxicity** robustly. It **cannot cleanly separate flux-driven from
neuronal-load-driven toxicity** on recovery *latency* — the discriminating
quantity there is non-identifiable from the available training data. It *can*
separate them on recovery *depth-dependence*, which is parameter-robust, and
that points to a concrete successor experiment tied to `x_crit`.

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

## 3. The powered discriminator: recovery vs knockdown *depth*

`t01_depth_discriminator.py` / `fig_t01_depth_discriminator.png`. Response of D
(60 d after knockdown, relative to its pre-knockdown peak) versus residual
neuronal PrP, for both modes at ρ ∈ {0.06, 0.15, 0.5}:

| mode | recovers (net) up to residual PrP of | across 50× ρ |
|---|---|---|
| flux | **~53–56%** | barely moves |
| neuronal-load | **~26–39%** | moves, ordering held |

flux recovers at ~2× shallower knockdown than neuronal-load, and the ordering is
ρ-robust. Mechanism: flux responds to the *instantaneous* drop in conversion, so
it recovers even while replication is still supercritical (residual > x_crit =
33%); neuronal-load can only recover once standing `z_n` clears, which needs
replication driven *subcritical* — hence its threshold straddles x_crit. This is
a shape/threshold difference governed by x_crit, not by ρ, which is why it has
the power latency lacks.

**Successor experiment (replaces "measure ρ", which is not actionable — ρ is a
lumped phenomenological rate).** The depth-titration survival study already
proposed in `RESEARCH_PLAN.md` §6 discriminates the hypotheses: titrate residual
neuronal PrP to ~50%, ~40%, ~30%, ~20% and look at where functional recovery
appears. Recovery at a *moderate* depth (~50% residual) is consistent with flux
and inconsistent with neuronal-load; recovery only at *deep* knockdown (<~30%
residual) favours neuronal-load. Either way the readout is knockdown depth
relative to x_crit — the number the whole project turns on.

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
