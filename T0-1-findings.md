# T0-1 findings — is the flux toxicity assumption load-bearing?

**Headline.** The two-compartment model **falsifies standing-total-load
toxicity** robustly. It **cannot separate flux-driven from
neuronal-load-driven toxicity** by any reversal experiment — not latency, not
depth. The two are *nested*: flux is the fast-clearance limit of neuronal-load,
so they differ only in the neuronal PrP-Sc clearance rate, which the training
data does not pin. The one measurement that separates them is a direct neuronal
PrP-Sc clearance timecourse — the neuronal analogue of the Mallucci 2003
observation that already killed the total-load model.

*(This file has been corrected twice under review. Earlier drafts claimed first a
latency discriminator, then a depth discriminator; both were artifacts of
leaving the clearance rate `a` pinned. The nesting result below is why no such
discriminator can exist. The correction history is kept visible rather than
rewritten away — see §3.)*

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

## 3. flux and neuronal-load are *nested*, not distinguishable by depth — a
retraction

The previous version of this file (and the message that accompanied it) claimed
recovery-vs-knockdown-*depth* was a "powered, ρ-robust discriminator": flux
recovers at ~2× shallower knockdown than neuronal-load. **That claim was wrong,
and wrong for the same reason twice.** It was computed with polymer clearance
`a` pinned at its default 0.02 — the exact parameter established two steps
earlier as one that must be free (it sets how fast standing `z_n` falls, and the
survival anchors, all at ~50% knockdown, do not pin it). The pin made
neuronal-load look incapable of fast or partial-depth recovery. Freed, the
separation collapses (`t01_depth_discriminator.py`,
`fig_t01_depth_discriminator.png`):

| neuronal-load clearance | clause A (>50% D-drop in 40 d) met at | dissociation (D-drop while supercritical) |
|---|---|---|
| a=0.02 (1/a=50 d) | never | 1% |
| a=0.12 (1/a=8 d) | never | 7% |
| a=0.20 (1/a=5 d) | ≤37% residual | **53%** |
| a=0.30 (1/a=3 d) | ≤52% residual | **69%** |
| *flux (reference)* | ≤25% residual | 38% |

At fast clearance neuronal-load **matches or exceeds** flux on both the
depth-threshold and the dissociation. The reason is structural: when clearance
is fast, `z_n` is slaved to the instantaneous conversion,
`z_n ≈ conv_n/a = β·x_n·y_n/a`, so `D_load = κ·z_n ∝ β·x_n·y_n` — the flux
driver. After κ is refit, **fast-clearance neuronal-load *is* flux. Flux is the
a→∞ limit of neuronal-load.** They are nested models, and no depth, latency, or
dissociation measurement separates nested models along the axis that nests them.

**These were never two hypotheses.** "Toxicity is driven by neuronal conversion
flux" and "toxicity is driven by neuronal standing load" are the same model
evaluated at two ends of one parameter — the clearance rate `a` — with flux at
`a→∞`. This is not a discrimination problem awaiting a cleverer experiment; it
is the discovery that **the question was mis-posed.** The entire empirical
content of the toxicity question is the *value of `a`*: how fast neuronal PrP-Sc
clears once conversion stops. Everything else in the toxicity layer is
phenomenology that cannot move the answer. That reframing — one model, one
unmeasured number — is the useful output of T0-1, more so than either mode
"winning" would have been.

So the honest status:
- **load_total: falsified** (robust; §1).
- **flux and neuronal-load: one model at two values of `a`.** No reversal-timing
  or reversal-depth experiment separates them, because those all move with `a`;
  the answer is `a` itself.

**The measurement that is the answer, and when it is decisive (`a*`).** The one
readout that pins `a` is the neuronal analogue of Mallucci 2003: after
neuron-specific knockdown, does **neuronal** PrP-Sc persist while behaviour
recovers, or does recovery track its fall? Mallucci established this dissociation
for *glial* PrP-Sc (which is what killed load_total); the intraneuronal version
is `a` measured directly.

`t01_astar.py` locates when that measurement is decisive. The flux signature
(behaviour recovered while neuronal PrP-Sc still elevated) is resolvable only if
z_n persists at least one RT-QuIC sampling interval past behavioural recovery:

> **a\* ≈ 0.048 /day — neuronal PrP-Sc half-life ≈ 14 days** (weekly sampling,
> 2-fold assay resolution, deep knockdown).
> - half-life **> ~14 d** (`a < a*`): the flux signature is resolvable; the
>   experiment discriminates. z_n still ≥50% of baseline for ~11–38 d after
>   behaviour recovers.
> - half-life **< ~14 d** (`a ≥ a*`): neuronal PrP-Sc clears as fast as behaviour
>   recovers; flux and neuronal-load are observationally identical — the
>   `a→∞`-limit regime.

PrP-Sc is protease-resistant and accumulates over the ~13-week disease, so a
half-life above 14 d (the discriminating regime) is the biologically expected
case — but this is now a stated, checkable inequality, not an assumption. The
moment anyone reports a neuronal PrP-Sc clearance rate, compare it to `a*`.

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
- **Step 4 (re-anchoring the toxicity layer) is consciously deferred, not
  skipped.** The plan called for anchoring `D_sx` and ρ to measured onset timing
  so the toxicity layer would not be free-floating. That work is *not worth
  doing for this question*: the non-identifiability lives entirely in `a`
  (§3), and no amount of constraint on `D_sx`, ρ, `D_tox`, or μ touches `a`.
  Re-anchoring would produce better-calibrated numbers on a question whose answer
  is already "unanswerable from reversal data" — motion, not progress. The
  toxicity layer remains phenomenology whose one decision-relevant parameter
  (`a`) is unmeasured; tightening the others dresses that up without changing it.
  Deferred deliberately, reason recorded. (What was checked: the *existing* `D_sx`
  already yields an emergent onset — 76 d, 51% of course — consistent with the
  early-behavioural-deficit natural history in `CITATIONS.md`, so the layer is
  not obviously mis-set; it is simply under-determined on the axis that matters.)

## What the power check bought

Catching the κ/ρ non-identifiability *before* Step 8 meant the ~2–3 h multi-start
fit + sweep was never spent on an underpowered test, and the naive
"flux-matches-Mallucci-so-flux-wins" conclusion — which the original
single-compartment figures gesture toward — was avoided. The substantive
findings are the robust death of load_total and the **nesting** of flux inside
neuronal-load, which together say precisely how much the reversibility data can
and cannot establish.

The process cost worth recording honestly: I twice proposed a discriminator
(latency, then depth) that turned out to be an artifact of leaving the clearance
rate `a` pinned, and was corrected both times under review. The nesting result
in §3 is the general reason those attempts were doomed — flux and neuronal-load
are the same model in the fast-clearance limit — and it is more useful than
either false discriminator would have been.

## Consequence for the repo's claims

Nothing here weakens the **replication-threshold** result (the ~65–90% knockdown
headline), which lives in the validated replication layer (`test_growth_rate.py`
confirms the determinant condition; `REPRODUCTION.md` confirms the calibration).
It sharpens the **toxicity** story: standing-*total*-load toxicity is falsified;
flux and neuronal-load are observationally equivalent along every axis a
reversal experiment can probe (they are nested), so the model **cannot** claim
flux is proven; and the one measurement that would separate them is a neuronal
PrP-Sc clearance timecourse, not a depth titration. The `~65–90%` knockdown
conclusion stands on the replication layer alone and does not depend on which
toxicity hypothesis is correct.
