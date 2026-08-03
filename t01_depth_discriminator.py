"""
t01_depth_discriminator.py
==========================
How the toxicity modes separate under a knockdown-DEPTH titration, resolved
against the actual pre-registered reversal clauses rather than a single latent
D metric. (An earlier version reported a "~55% residual recovery threshold"
using only whether D nets a drop at +60 d -- the weakest possible notion. That
conflated transient dysfunction relief with rescue; this version reports each
clause separately.)

Pre-registered reversal clauses:
  A  D falls >50% from its pre-knockdown peak within 40 d
  B  total load z_n+z_g still rising across that window
  C  survival past 3x untreated (>450 d)

Knockdown applied at day 90 (60% of the 150 d course), depth swept.

WHAT SEPARATES THE MODES, AND WHAT DOES NOT
-------------------------------------------
* Survival rescue (clause C) requires driving replication subcritical -- residual
  PrP below x_crit -- for BOTH modes. A survival titration therefore locates
  x_crit; it does NOT by itself reveal the toxicity mechanism.
* Fast behavioural recovery (clause A, D down >50% within 40 d) separates them,
  but load_neuronal's ability to show it depends on the polymer clearance rate
  a, which training does not pin -- the same non-identifiability that sank the
  latency test. So clause A alone is not a robust discriminator either.
* The robust, x_crit-anchored discriminator is the DISSOCIATION: does behaviour
  improve at a knockdown depth where replication is still SUPERCRITICAL
  (residual > x_crit, seeding still rising)? Flux permits this transient
  dissociation; neuronal-load cannot -- its D can only fall once z_n clears,
  which needs subcritical replication.

THE TWO-UNKNOWNS PROBLEM AND ITS FIX
------------------------------------
An observed behavioural-recovery threshold at, say, 30% residual is equally
consistent with flux at x_crit=15% and neuronal-load at x_crit=30%: one
measurement, two unknowns (recovery depth AND x_crit, which the original
analysis pins only loosely to 10-35% residual). The fix, already half-present in
RESEARCH_PLAN section 6: add an RT-QuIC seeding readout. Seeding-decline depth
gives x_crit from REPLICATION alone, independent of any toxicity assumption;
behavioural-recovery depth gives the toxicity threshold. Recovery at a shallower
depth than seeding-decline => flux; the two coincide => neuronal-load.

Output: fig_t01_depth_discriminator.png + printed clause table.
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
                         x_crit, IX2C, z_total_2c)

HERE = os.path.dirname(os.path.abspath(__file__))
cal = json.load(open(os.path.join(HERE, "results.json")))["calibration"]

T_TREAT = 90.0
WIN = 40.0
RESCUE_DAY = 450.0
RESIDUALS = np.linspace(0.05, 0.65, 25)


def refit_kappa(p):
    def f(k):
        st = survival_time_2c(replace(p, kappa=k), t_end=4000.0)
        return (st if np.isfinite(st) else 4000.0) - 150.0
    try:
        return brentq(f, 1e-6, 1e12, xtol=1e-9, rtol=1e-8, maxiter=300)
    except ValueError:
        return None


def calibrate(mode, rho=0.06):
    p = replace(PrionParams(), beta=cal["beta"], rho=rho, toxicity_mode=mode)
    k = refit_kappa(p)
    return replace(p, kappa=k) if k is not None else None


def clause_metrics(p, xr):
    """Return (D_drop_frac, load_rising, survival_day) at residual PrP xr."""
    kd = 1.0 - xr
    sol = simulate_2c(p, t_end=T_TREAT + WIN + 5, kd_n=kd, kd_g=0.0,
                      t_treat=T_TREAT, ramp_days=7.0, max_step=1.0)
    t, D = sol.t, sol.y[IX2C["D"]]
    zt = z_total_2c(sol)
    pre = t <= T_TREAT + 1.0
    D_peak = float(D[pre].max())
    win = (t > T_TREAT) & (t <= T_TREAT + WIN)
    D_drop = (D_peak - float(D[win].min())) / D_peak if win.any() else 0.0
    i_min = np.where(win)[0][int(np.argmin(D[win]))] if win.any() else 0
    load_rising = bool(zt[i_min] > zt[int(np.argmin(np.abs(t - T_TREAT)))])
    surv = survival_time_2c(p, t_end=2000.0, kd_n=kd, kd_g=0.0, t_treat=T_TREAT,
                            ramp_days=7.0, max_step=2.0)
    return D_drop, load_rising, surv


if __name__ == "__main__":
    xc = x_crit(replace(PrionParams(), beta=cal["beta"]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)

    print(f"x_crit = {xc:.3f} residual PrP ({100*(1-xc):.0f}% knockdown); "
          f"knockdown at day {T_TREAT:.0f}. Seeding declines below x_crit.\n")
    for mode, ax in zip(["flux", "load_neuronal"], axes):
        p = calibrate(mode)
        drops, survs = [], []
        for xr in RESIDUALS:
            d, _, s = clause_metrics(p, xr)
            drops.append(100 * d)
            survs.append(min(s, 600.0) if np.isfinite(s) else 600.0)
        drops = np.array(drops); survs = np.array(survs)

        ax.plot(100 * RESIDUALS, drops, "o-", ms=3, color="#1d4ed8",
                label="behavioural recovery\n(% D-drop in 40 d)")
        ax.axhline(50, color="#1d4ed8", ls=":", lw=0.8)
        ax2 = ax.twinx()
        ax2.plot(100 * RESIDUALS, survs, "s-", ms=3, color="#c2410c",
                 label="survival (d, capped 600)")
        ax2.axhline(RESCUE_DAY, color="#c2410c", ls=":", lw=0.8)
        ax2.set_ylim(100, 640)
        ax.axvspan(100 * xc, 100 * RESIDUALS.max(), color="gray", alpha=0.10)
        ax.axvline(100 * xc, color="k", ls="--", lw=1)
        ax.text(100 * xc + 1, 90, "supercritical\n(seeding rising)", fontsize=7,
                va="top")
        ax.set_title(mode, weight="bold")
        ax.set_xlabel("residual neuronal PrP after knockdown (%)")
        ax.set_ylim(0, 100)
        if mode == "load_neuronal":
            ax2.set_ylabel("survival (days)", color="#c2410c")
        ax.legend(loc="upper left", fontsize=7, frameon=False)
        ax2.legend(loc="lower right", fontsize=7, frameon=False)

        # report where each clause turns on
        A = 100 * RESIDUALS[drops > 50]
        C = 100 * RESIDUALS[survs > RESCUE_DAY]
        print(f"{mode}: clause A (>50% D-drop) met at residual <= "
              f"{A.max():.0f}%" if len(A) else f"{mode}: clause A NEVER met")
        print(f"        clause C (survival>450d) met at residual <= "
              f"{C.max():.0f}%  (x_crit={100*xc:.0f}%)" if len(C)
              else f"        clause C never met")
        # dissociation: behaviour improves while supercritical (residual>xc)
        sup = RESIDUALS > xc
        diss = drops[sup].max() if sup.any() else 0.0
        print(f"        max behavioural improvement while SUPERCRITICAL: "
              f"{diss:.0f}%  -> {'dissociation present (flux-like)' if diss>15 else 'no dissociation (load-like)'}\n")

    axes[0].set_ylabel("behavioural recovery (% D-drop in 40 d)", color="#1d4ed8")
    fig.suptitle("Depth titration: survival rescue locates x_crit for both modes; "
                 "dissociation above x_crit separates them", weight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_t01_depth_discriminator.png"), dpi=140)
    print("Wrote fig_t01_depth_discriminator.png")
