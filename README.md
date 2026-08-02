# Prion reversal: a quantitative threshold model

Companion code for the question "can prion disease be reversed or cured?"

## Files
- `prion_model.py` — nucleated polymerization (Masel–Jansen–Nowak) replication
  layer with an analytic self-sustaining threshold, plus a toxicity layer in
  which damage is driven by *conversion flux* rather than standing PrP-Sc load.
- `run_analysis.py` — three analyses: threshold inference from published
  dose-response data, the predicted survival curve, and the depth × timing
  reversibility window.
- `RESEARCH_PLAN.md` — what's settled, what's open, and concrete entry points.
- `fig1`–`fig4`, `results.json` — outputs.

## Run
```
pip install numpy scipy matplotlib
python run_analysis.py
```

## Headline results
- Required PrP knockdown to push replication below self-sustaining: **~65–90%**,
  depending on which published anchor is used. 50% is on the wrong side.
- Calibrated to a 150-day RML course, 85% knockdown still permits **full rescue
  up to ~day 98** — past clinical onset. This was not fitted to Mallucci's data
  but lands close to it (depletion at ~9 of 13 weeks post-inoculation).
- Survival is strongly non-linear in knockdown depth, with a vertical asymptote
  at `x_crit`. This is why 50%-lowering studies all look similar and all fail.

## Key sources
- Mallucci et al., *Science* 302:871 (2003); *Neuron* 53:325 (2007) — reversal of
  spongiosis and behavioural deficits on neuronal PrP depletion.
- Masel, Jansen & Nowak, *Biophys Chem* 77:139 (1999) — nucleated polymerization kinetics.
- Sandberg et al., *Nature* 470:540 (2011) — two-phase replication/toxicity kinetics.
- Büeler et al. (1993, 1994) — Prnp-/- resistance; Prnp+/- extended incubation.
- Minikel et al., *Sci Transl Med* 8:322ra9 (2016) — PRNP penetrance; code at
  github.com/ericminikel/prnp_penetrance.
- Minikel et al. (2020) — PrP-lowering ASOs delay onset and slow progression.
- An, Davis et al., *Nat Med* (2025) — in vivo base editing, 50% PrP reduction,
  +52% lifespan in humanized PRNP mice.
- Gentile et al., *Nucleic Acids Res* 54:gkag287 (2026) — divalent siRNA; 2.7×
  survival at 49% residual PrP; 17% residual achievable with a single dose.
- Ionis PrProfile (NCT06153966) — ION717 Phase 1/2a; reopened March 2026 for a
  third, higher dosing regimen.

## Caveat
Lumped-parameter toy model. Every parameter is uncertain and several are
unidentifiable from existing data — which is precisely the argument for running
the depth-titration experiment described in `RESEARCH_PLAN.md`.
