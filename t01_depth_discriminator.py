"""
t01_depth_discriminator.py
==========================
Does a knockdown-DEPTH titration separate flux from neuronal-load? Earlier
drafts of this file claimed yes. With the polymer-clearance rate a held at its
default 0.02 the two look cleanly separable -- but a is exactly the parameter
that must be free (it sets how fast standing z_n falls), and training does not
pin it (every survival anchor sits at ~50% knockdown). Freed, the claimed
discriminator collapses. This script shows the collapse honestly.

THE STRUCTURE
-------------
Under neuronal-load, D tracks z_n. When clearance is fast, z_n is slaved to the
instantaneous conversion: dz_n/dt ~ 0 gives z_n ~ conv_n/a = beta*x_n*y_n/a, so
D_load = kappa*z_n is proportional to the flux driver beta*x_n*y_n. After kappa
is refit, fast-clearance neuronal-load IS flux. flux is the a->infinity limit of
neuronal-load; the two are nested, not distinct.

So the depth-response of neuronal-load migrates onto the flux response as a
grows. At small a (slow clearance) they differ -- neuronal-load recovers slowly
and shows no dissociation; at large a they coincide. The modes are therefore
separable ONLY if neuronal PrP-Sc clearance is independently known to be slow.

WHAT THIS LEAVES AS THE ACTUAL DISCRIMINATOR
--------------------------------------------
Not depth, not latency, not the dissociation alone -- all of these move with a.
The one thing that pins a is a direct measurement: after neuronal conversion is
switched off, does NEURONAL PrP-Sc persist while behaviour recovers (flux-like),
or does behavioural recovery track the fall of neuronal PrP-Sc (load-like)?
This is the neuronal analogue of Mallucci 2003's extraneuronal observation, and
it is what RESEARCH_PLAN section 6 now specifies.

Output: fig_t01_depth_discriminator.png + printed table.
"""

import json
import os
import numpy as np
from dataclasses import replace
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from prion_model import (PrionParams, simulate_2c, survival_time_2c,
                         solve_beta_for_target, x_crit, IX2C)

HERE = os.path.dirname(os.path.abspath(__file__))
cal = json.load(open(os.path.join(HERE, "results.json")))["calibration"]

T_TREAT = 90.0
WIN = 40.0
RESIDUALS = np.linspace(0.05, 0.60, 23)


def refit_kappa(p):
    def f(k):
        st = survival_time_2c(replace(p, kappa=k), t_end=4000.0)
        return (st if np.isfinite(st) else 4000.0) - 150.0
    try:
        return brentq(f, 1e-6, 1e12, xtol=1e-9, rtol=1e-8, maxiter=300)
    except ValueError:
        return None


def d_drop_curve(p):
    """Max fractional D-drop within WIN days of knockdown, vs residual PrP."""
    out = []
    for xr in RESIDUALS:
        sol = simulate_2c(p, t_end=T_TREAT + WIN + 5, kd_n=1.0 - xr, kd_g=0.0,
                          t_treat=T_TREAT, ramp_days=7.0, max_step=1.0)
        t, D = sol.t, sol.y[IX2C["D"]]
        pre = t <= T_TREAT + 1.0
        Dpk = float(D[pre].max())
        win = (t > T_TREAT) & (t <= T_TREAT + WIN)
        out.append(100 * (Dpk - float(D[win].min())) / Dpk if win.any() else 0.0)
    return np.array(out)


if __name__ == "__main__":
    # flux reference (a is irrelevant to flux; use default beta)
    pf = replace(PrionParams(), beta=cal["beta"], rho=0.06, toxicity_mode="flux")
    pf = replace(pf, kappa=refit_kappa(pf))
    xc = x_crit(pf)
    flux_curve = d_drop_curve(pf)

    # neuronal-load for a family of clearance rates a (beta refit per a so
    # x_crit is held; kappa refit to terminal)
    a_values = [0.02, 0.05, 0.12, 0.20, 0.30]
    load_curves = {}
    print(f"x_crit = {xc:.3f} residual PrP. Knockdown at day {T_TREAT:.0f}.")
    print("Behavioural recovery (% D-drop in 40d) and dissociation "
          "(max D-drop while residual>x_crit):\n")
    print(f"{'mode/a':>16} {'clauseA @resid<=':>16} {'dissociation':>13}")
    A = RESIDUALS[flux_curve > 50]
    sup = RESIDUALS > xc
    print(f"{'flux':>16} {(str(int(100*A.max()))+'%') if len(A) else 'never':>16} "
          f"{100 if False else int(flux_curve[sup].max()):>12}%")
    for a in a_values:
        beta = solve_beta_for_target(replace(PrionParams(), a=a), 0.50, 2.70,
                                     t_tox_frac=0.25)
        p = replace(PrionParams(), a=a, beta=beta, rho=0.06,
                    toxicity_mode="load_neuronal")
        p = replace(p, kappa=refit_kappa(p))
        load_curves[a] = d_drop_curve(p)
        Al = RESIDUALS[load_curves[a] > 50]
        print(f"{'load a=' + str(a):>16} "
              f"{(str(int(100*Al.max()))+'%') if len(Al) else 'never':>16} "
              f"{int(load_curves[a][sup].max()):>12}%")

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.plot(100 * RESIDUALS, flux_curve, "k-", lw=2.5, label="flux (reference)")
    cmap = plt.cm.viridis(np.linspace(0.15, 0.9, len(a_values)))
    for a, c in zip(a_values, cmap):
        ax.plot(100 * RESIDUALS, load_curves[a], "o-", ms=3, color=c,
                label=f"neuronal-load, a={a} (1/a={1/a:.0f}d)")
    ax.axvspan(100 * xc, 100 * RESIDUALS.max(), color="gray", alpha=0.10)
    ax.axvline(100 * xc, color="k", ls="--", lw=1)
    ax.text(100 * xc + 1, 30, "x_crit; right of here\nreplication supercritical\n(seeding still rising)",
            fontsize=7, va="top")
    ax.axhline(50, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("residual neuronal PrP after knockdown (%)")
    ax.set_ylabel("behavioural recovery (% D-drop in 40 d)")
    ax.set_title("Neuronal-load migrates onto flux as clearance a increases:\n"
                 "depth does NOT separate the modes unless a is known to be slow",
                 fontsize=10, weight="bold")
    ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_t01_depth_discriminator.png"), dpi=140)
    print("\nWrote fig_t01_depth_discriminator.png")
    print("-> at a=0.02 load looks nothing like flux; by a=0.30 it reproduces "
          "flux's depth-response. Not separable without pinning a.")
