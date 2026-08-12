# T0-1 — Try to break the toxicity model

**Estimated session length:** 1–2 hours. Run this one first; it is fully
self-contained and needs no external data.

## Objective

Determine whether the flux-driven toxicity assumption in `prion_model.py` is
*load-bearing* or merely convenient, by building a competing load-driven variant
and testing both against a held-out observation neither was fitted to.

## Why it matters

The current model assumes neuronal damage accrues in proportion to the ongoing
PrP-C → PrP-Sc conversion flux (`kappa * beta * x * y`). The obvious
alternative — the one most people intuitively assume — is that damage accrues in
proportion to standing PrP-Sc burden (`kappa * z`).

These two assumptions make identical predictions about untreated survival, so
survival data alone cannot distinguish them. They diverge sharply on one
observation: Mallucci's 2003 result that depleting neuronal PrP-C reverses
spongiosis *while PrP-Sc continues to accumulate*. Under the load-driven model
that should be impossible.

If the load-driven variant can be made to reproduce Mallucci's reversal after
fair refitting, the flux assumption is doing no work and the model's conclusions
about reversibility are weaker than they look. That is worth knowing.

## Method

1. **Pre-register first.** Write `T0-1-preregistration.md` containing:
   - the exact quantitative criteria for "reproduces Mallucci reversal"
     (suggest: dysfunction `D` must fall by >50% from its peak within 40 days of
     knockdown onset, while `z` continues to rise, and the animal must survive
     past 3× untreated lifespan);
   - which observations are used for fitting (the 50%-lowering survival anchors
     in `ANCHORS`) and which are held out (the Mallucci reversal);
   - what result would count as falsifying the flux model instead.
   Commit this before writing any model code.

2. **Implement the variant.** Add a `toxicity_mode` field to `PrionParams`
   taking `"flux"` or `"load"`, and branch in `_rhs`. Do not fork the file.

3. **Refit fairly.** Give the load-driven variant the same number of free
   parameters and refit `kappa`, `rho`, `D_tox`, `mu` to the same untreated
   course (onset ~day 77, terminal ~day 150) and the same 50%-lowering survival
   anchors. Use a proper optimiser, not hand-tuning. Log the fit quality.

4. **Test on held-out data.** Simulate the Mallucci protocol under both models:
   85–100% neuronal PrP knockdown starting at ~70% of the untreated course.
   Record whether `D` falls while `z` rises.

5. **Sweep for escape hatches.** Before concluding the load model fails, search
   its parameter space properly (Latin hypercube or similar, ≥2000 samples) for
   *any* combination that satisfies both the fitting anchors and the held-out
   criterion. Report the sampled volume and the fraction that pass. "I couldn't
   find one by hand" is not a result; "0/2000 samples passed" is.

## Deliverables

- `T0-1-preregistration.md` (committed before results exist)
- Modified `prion_model.py` with both toxicity modes
- `t01_model_comparison.py`
- `fig_t01_model_comparison.png` — side-by-side trajectories under both models
- `T0-1-findings.md` — one page, stating plainly which model survived and what
  the residual uncertainty is

## Acceptance criteria

- Pre-registration file is committed in an earlier commit than any results file
  (verifiable via `git log`).
- The parameter sweep is reproducible from a fixed random seed.
- `T0-1-findings.md` explicitly names which observations were fitted and which
  were held out.
- If the flux model *loses*, the findings file says so in the first paragraph.

## Traps

- Do not let the load-driven model cheat by making `a` (polymer clearance) large
  enough that `z` falls after knockdown. Mallucci's observation is specifically
  that `z` **rose**. Constrain the sweep so held-out simulations have `dz/dt > 0`
  throughout the post-knockdown window, or the test is vacuous.
- Do not compare models by eyeballing the figures. Use the pre-registered
  numeric criteria.
