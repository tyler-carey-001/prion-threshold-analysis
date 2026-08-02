# Reproduction of the original analysis

Run before any modification, on unmodified `prion_model.py` / `run_analysis.py`.
`results.json` is a `run_analysis.py` output and was not in the download batch,
so it is regenerated here and committed as the frozen regression baseline.

**Environment:** Python 3.10, numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9,
macOS (darwin 23.6.0).

## Calibration constants

Tolerance was fixed in advance at **~3 significant figures**. This is unmodified
original code, not a rebuild, so anything looser would let real numerical
fragility pass as "close enough" — and a miss at 3 sig figs would itself be a
finding, implying the `brentq` solve for `kappa` (`run_analysis.py:188`) is
numerically fragile and demoting the terminal-day check below even a convergence
diagnostic.

| quantity | expected | reproduced | agreement |
|---|---|---|---|
| `beta` | 6.67 | 6.665562105385422 | 4 s.f. |
| `kappa` | 2051.6 | 2051.5723240209686 | 6 s.f. |
| `x_crit` | 0.334 | 0.3339553311192629 | 4 s.f. |

**Result: PASS, comfortably.** All three exceed the 3-s.f. bar. The
root-finder is not fragile: terminal day solves to 150.000000591, i.e. `brentq`
converged to ~6e-7 d against a 150 d target.

Derived quantities: required-knockdown range **64.3–91.9%** (`README`/
`RESEARCH_PLAN` quote "roughly 65–90%"), untreated onset **76.46 d**, terminal
**150.0 d**, last day of full rescue at 85% knockdown **97.5 d**.

## Figures

Pixel-diffed against the shipped PNGs:

| figure | differing pixels | interpretation |
|---|---|---|
| fig1 threshold inference | 0 | bit-identical |
| fig2 dose response | 0 | bit-identical |
| fig3 reversibility window | 643 px (0.113%) | ~2 cells of the 16×16 outcome grid flip class |
| fig4 trajectories | 326 px (0.054%), max Δ 0.13 | antialiasing only |

fig3's differing pixels sit in two contiguous blocks in the high-knockdown
region of the heatmap, consistent with a small number of grid cells landing on
the other side of an outcome-category boundary — the categories are decided by
event detection (`survival_time` finite vs infinite) and a hard `N < 0.02`
cutoff, so cells adjacent to a boundary are expected to be sensitive to
integrator version. This does not affect any calibration constant.

## What this baseline is for

Every subsequent change is checked against `results.json`. Note that the
untreated onset of ~76.5 d is the **old model's output**, not a measured
quantity — see `T0-1-preregistration.md` on why re-running it after the
compartment split is a regression check, not a validity check.
