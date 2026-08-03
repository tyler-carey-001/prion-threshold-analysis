"""
t01_astar.py
============
The nesting result (flux = a->infinity limit of neuronal-load) says the toxicity
question reduces to a single unmeasured number: the neuronal PrP-Sc clearance
rate a. The clearance-timecourse experiment (neuronal analogue of Mallucci 2003)
discriminates flux from neuronal-load ONLY when a is slow enough that neuronal
PrP-Sc still persists after behaviour has recovered. This script locates that
boundary, a*, so "conditional on slow clearance" becomes a checkable inequality
a < a* against any future measurement.

DEFINITION (operational, tied to the assay, not to the latent variable)
-----------------------------------------------------------------------
Under the flux model, apply a deep neuron-specific knockdown (the Mallucci
protocol). Two clocks then run:
  * behaviour recovers on the repair timescale ~1/rho (a-independent);
  * neuronal PrP-Sc z_n clears on the clearance timescale ~1/a.
The flux SIGNATURE is that behaviour recovers while z_n is still elevated. It is
observable only if z_n stays above the assay's resolution for at least one
sampling interval AFTER behaviour has recovered.

Assay assumptions (state them; they set a*):
  SAMPLING_INTERVAL  RT-QuIC cadence in a mouse timecourse (default weekly).
  DETECT_FRACTION    fold-resolution: z_n counts as "still present" while above
                     this fraction of its pre-knockdown level (default 0.5,
                     i.e. a 2-fold drop is the smallest reliably-called change).
  RECOVER_FRACTION   behaviour counts as recovered when D has fallen this far
                     from its peak (default 0.5).

a* is the largest a for which z_n is still >= DETECT_FRACTION of baseline for at
least SAMPLING_INTERVAL days after behaviour has recovered. a < a* => the flux
signature is resolvable and the experiment discriminates; a >= a* => neuronal
PrP-Sc clears too fast to catch the dissociation, and flux and neuronal-load are
observationally identical.
"""

import json
import os
import numpy as np
from dataclasses import replace
from scipy.optimize import brentq

from prion_model import (PrionParams, simulate_2c, survival_time_2c,
                         solve_beta_for_target, x_crit, IX2C)

HERE = os.path.dirname(os.path.abspath(__file__))
cal = json.load(open(os.path.join(HERE, "results.json")))["calibration"]

T_TREAT = 90.0
KD_N = 0.90                      # deep neuron-specific knockdown (Mallucci-like)
SAMPLING_INTERVAL = 7.0         # weekly RT-QuIC
DETECT_FRACTION = 0.5           # 2-fold resolution
RECOVER_FRACTION = 0.5          # behaviour recovered = D down 50%


def refit_kappa(p):
    def f(k):
        st = survival_time_2c(replace(p, kappa=k), t_end=4000.0)
        return (st if np.isfinite(st) else 4000.0) - 150.0
    try:
        return brentq(f, 1e-6, 1e12, xtol=1e-9, rtol=1e-8, maxiter=300)
    except ValueError:
        return None


def flux_model(a, rho=0.06):
    """flux model at clearance a, with beta refit so x_crit is held and kappa
    refit to terminal=150."""
    beta = solve_beta_for_target(replace(PrionParams(), a=a), 0.50, 2.70,
                                 t_tox_frac=0.25)
    if beta is None:
        return None
    p = replace(PrionParams(), a=a, beta=beta, rho=rho, toxicity_mode="flux")
    k = refit_kappa(p)
    return replace(p, kappa=k) if k is not None else None


def persistence_window(a, rho=0.06):
    """Days that z_n stays >= DETECT_FRACTION of baseline AFTER behaviour has
    recovered, under the flux model. Positive => flux signature resolvable."""
    p = flux_model(a, rho)
    if p is None:
        return None
    sol = simulate_2c(p, t_end=T_TREAT + 400.0, kd_n=KD_N, kd_g=0.0,
                      t_treat=T_TREAT, ramp_days=7.0, max_step=0.5)
    t = sol.t
    D = sol.y[IX2C["D"]]
    zn = sol.y[IX2C["z_n"]]
    pre = t <= T_TREAT + 1.0
    D_peak = float(D[pre].max())
    zn_base = float(zn[np.argmin(np.abs(t - T_TREAT))])

    # time behaviour recovers (D <= RECOVER_FRACTION * peak)
    rec = np.where((t > T_TREAT) & (D <= RECOVER_FRACTION * D_peak))[0]
    if not len(rec):
        return None
    t_rec = t[rec[0]]
    # time z_n falls below DETECT_FRACTION of baseline
    clr = np.where((t > T_TREAT) & (zn <= DETECT_FRACTION * zn_base))[0]
    t_clr = t[clr[0]] if len(clr) else t[-1]
    return t_clr - t_rec, t_rec - T_TREAT, zn[rec[0]] / zn_base


def window_at(a, sampling_interval, detect_fraction, recover_fraction):
    """Persistence window (days z_n stays >= detect_fraction of baseline past
    behavioural recovery) for given assay assumptions."""
    p = flux_model(a)
    if p is None:
        return None
    sol = simulate_2c(p, t_end=T_TREAT + 400.0, kd_n=KD_N, kd_g=0.0,
                      t_treat=T_TREAT, ramp_days=7.0, max_step=0.5)
    t, D, zn = sol.t, sol.y[IX2C["D"]], sol.y[IX2C["z_n"]]
    pre = t <= T_TREAT + 1.0
    D_peak = float(D[pre].max())
    zn_base = float(zn[np.argmin(np.abs(t - T_TREAT))])
    rec = np.where((t > T_TREAT) & (D <= recover_fraction * D_peak))[0]
    if not len(rec):
        return None
    t_rec = t[rec[0]]
    clr = np.where((t > T_TREAT) & (zn <= detect_fraction * zn_base))[0]
    t_clr = t[clr[0]] if len(clr) else t[-1]
    return t_clr - t_rec


def a_star(sampling_interval, detect_fraction, recover_fraction):
    """Clearance a at which the persistence window equals one sampling interval.
    Below it the flux signature is resolvable; above it, not."""
    def g(a):
        w = window_at(a, sampling_interval, detect_fraction, recover_fraction)
        return (w if w is not None else -1e3) - sampling_interval
    try:
        return brentq(g, 0.01, 0.3, xtol=1e-4, maxiter=100)
    except ValueError:
        return None


if __name__ == "__main__":
    print("a* is NOT a single number: it depends on assay assumptions that are")
    print("placeholders for the wet-lab protocol, not measured quantities.")
    print("Reported here as a RANGE across plausible assumptions, so 0.048 is not")
    print("mistaken for a threshold (same discipline as not reporting a false-")
    print("precision latency floor).\n")

    # sensitivity sweep over the three unpinned assay assumptions
    sampling = [3.5, 7.0, 14.0]          # twice-weekly, weekly, biweekly
    detect = [0.33, 0.5, 0.7]            # 3-fold, 2-fold, 1.4-fold resolution
    recover = [0.5]                      # "recovered" = D down 50% (definitional)

    print(f"{'sampling(d)':>12} {'detect':>7} {'a*':>8} {'half-life(d)':>13}")
    vals = []
    for si in sampling:
        for df in detect:
            ast = a_star(si, df, recover[0])
            if ast is None:
                print(f"{si:>12.1f} {df:>7.2f}   (none in range)")
                continue
            vals.append(ast)
            print(f"{si:>12.1f} {df:>7.2f} {ast:>8.3f} {np.log(2)/ast:>13.0f}")

    lo, hi = min(vals), max(vals)
    print()
    print(f"a* RANGE = [{lo:.3f}, {hi:.3f}] /day  "
          f"(neuronal PrP-Sc half-life ~{np.log(2)/hi:.0f}-{np.log(2)/lo:.0f} d)")
    print(f"  nominal (weekly, 2-fold): a* ~ {a_star(7.0,0.5,0.5):.3f}")
    print("  Interpretation: the clearance-timecourse experiment discriminates")
    print("  flux from neuronal-load when neuronal PrP-Sc clears slower than this")
    print("  band; faster, and the two are observationally identical. The band")
    print("  itself must be rebuilt on the actual assay before it is a threshold.")
    print("\nCheckable against literature: find any measurement of neuronal PrP-Sc")
    print("clearance after conversion arrest; compare to the band, not to 0.048.")
