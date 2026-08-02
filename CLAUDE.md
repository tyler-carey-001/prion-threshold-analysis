# CLAUDE.md

## What this project is

A computational investigation into whether prion disease can be halted or
reversed by lowering PrP-C, the obligate substrate for prion propagation.

The repo currently contains a lumped-parameter kinetic model
(`prion_model.py`, `run_analysis.py`) that combines:

- a **replication layer**: the Masel–Jansen–Nowak nucleated polymerization
  model, which has an analytic self-sustaining threshold `x_crit`;
- a **toxicity layer**: damage driven by *neuronal conversion flux*, not by
  standing PrP-Sc load.

Current headline output: halting replication appears to require ~65–90% PrP
knockdown, and 50% (where every published survival study clusters) is on the
wrong side of that line.

## Domain facts you must not get wrong

- PrP-C is encoded by **PRNP**, chromosome 20, and the entire open reading
  frame sits in a **single exon**. Premature stop codons therefore produce
  truncated protein rather than triggering nonsense-mediated decay. This matters
  for interpreting truncating variants.
- `Prnp-/-` mice are healthy and completely resistant to prion disease. Lowering
  PrP is a target validated by nature.
- **PrP-Sc burden and toxicity are dissociable.** Mallucci 2003 (Science
  302:871) depleted neuronal PrP-C mid-infection; spongiosis reversed and the
  animals survived *while extraneuronal PrP-Sc kept accumulating to terminal
  levels*. Any model that makes toxicity a function of standing load contradicts
  this. That is the crux of task T0-1.
- **Codon 129 (M129V, rs1799990)** is a major modifier and determines whether
  D178N presents as fatal familial insomnia (129M in cis) or CJD (129V in cis).
  Any penetrance analysis that ignores it is incomplete.
- Prion disease is ~85% sporadic, ~15% genetic. Roughly 12–18% of surveilled
  cases carry a rare PRNP variant.

## Working conventions

- Python 3, numpy/scipy/matplotlib/pandas. Add dependencies only when needed.
- Time units are **days** throughout the kinetic model.
- `x` is always PrP-C concentration normalised so untreated steady state = 1.0.
- Every figure gets a caption in the accompanying markdown that states what
  would falsify the claim it depicts.
- Commit after each logical step so diffs stay reviewable.

## Scientific integrity rules — these are hard constraints

1. **Never tune parameters until a desired result appears.** If a fit fails,
   report that it failed. A negative result is a result.
2. **Pre-register.** For any model comparison, write the acceptance criteria to
   a file and commit it *before* running the comparison.
3. **Separate fitted data from held-out data.** State explicitly which
   observations were used for fitting and which were reserved for testing.
4. **Do not invent citations.** If you are unsure a paper says something, either
   fetch and verify it or mark the claim `[UNVERIFIED]`. Never fabricate a PMID,
   DOI, effect size, or sample count.
5. **Distinguish model output from evidence.** This model is a hypothesis
   generator. Phrase conclusions as "the model predicts", never as "we show".
6. **No clinical interpretation.** This repo does not produce advice about any
   individual's genotype, risk, or care. If analysis touches on penetrance,
   frame results as population-level estimates with explicit uncertainty.

## Stop and ask the user if

- A task would require downloading individual-level genetic data.
- You cannot reproduce a published number and the discrepancy is >2-fold.
- A result looks strong enough to be worth publishing — the human should decide
  what happens next, not you.

## Out of scope

Anything involving physical handling of prions or infectious material. This is a
dry-lab repo. Prions resist autoclaving and standard disinfection, and there are
confirmed occupational fatalities. Wet work belongs in a properly equipped
institution with trained staff.
