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


if __name__ == "__main__":
    print("Locating a*: clearance below which the flux signature (behaviour "
          "recovered\nwhile neuronal PrP-Sc persists) is resolvable by RT-QuIC.")
    print(f"assumptions: weekly sampling, 2-fold detection, D-50% = recovered, "
          f"deep KD {KD_N:.0%}\n")
    print(f"{'a':>6} {'1/a (d)':>8} {'behav.recovers(d)':>18} "
          f"{'z_n at recovery':>16} {'persistence window(d)':>22}")
    aa = [0.02, 0.04, 0.06, 0.10, 0.15, 0.20, 0.30, 0.50]
    windows = {}
    for a in aa:
        r = persistence_window(a)
        if r is None:
            print(f"{a:>6.2f} {1/a:>8.1f}   (no calibration)")
            continue
        W, trec, zfrac = r
        windows[a] = W
        print(f"{a:>6.2f} {1/a:>8.1f} {trec:>18.1f} {100*zfrac:>15.0f}% "
              f"{W:>22.1f}")

    # a* where persistence window == one sampling interval
    xs = sorted(windows)
    a_star = None
    for i in range(len(xs) - 1):
        w0, w1 = windows[xs[i]], windows[xs[i + 1]]
        if (w0 - SAMPLING_INTERVAL) * (w1 - SAMPLING_INTERVAL) <= 0:
            # linear interpolation in log-a
            la0, la1 = np.log(xs[i]), np.log(xs[i + 1])
            a_star = np.exp(la0 + (SAMPLING_INTERVAL - w0) * (la1 - la0) / (w1 - w0))
            break

    print()
    if a_star:
        print(f"a* ~ {a_star:.3f} /day  (neuronal PrP-Sc half-life "
              f"~{np.log(2)/a_star:.0f} d)")
        print(f"  a < {a_star:.3f}  -> flux signature resolvable, experiment "
              f"discriminates flux from neuronal-load.")
        print(f"  a >= {a_star:.3f} -> neuronal PrP-Sc clears too fast to catch; "
              f"the two modes are observationally identical.")
    else:
        print("a* outside the scanned range.")
    print("\nCheckable against literature: find any measurement of neuronal "
          "PrP-Sc\nclearance kinetics after conversion arrest; compare to a*.")
