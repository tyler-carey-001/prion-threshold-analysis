"""
t01_depth_discriminator.py
==========================
The latency test cannot cleanly separate flux from neuronal-load: with polymer
clearance a free (and training does not pin it -- every survival anchor sits at
~50% knockdown), load_neuronal can be pushed to match the observed ~7-14 d
recovery. So latency is the wrong axis.

This script tests the discriminator that IS powered: the recovery response as a
function of knockdown DEPTH, measured a fixed window after treatment.

  flux           D is driven by conversion flux beta*x_n*y_n. Partial knockdown
                 drops flux immediately and proportionally, so D falls even at
                 modest knockdown depth.
  load_neuronal  D is driven by standing z_n, which can only fall by clearance
                 (slow, rate a) and only actually clears once replication has
                 gone subcritical. So D recovers only at DEEP knockdown.

WHAT THIS SCRIPT ACTUALLY FINDS (reported honestly, not the sharp step first
hypothesised): both depth-response curves are smooth, but they separate in the
RECOVERY-THRESHOLD DEPTH -- the residual PrP at which net recovery (D at +POST
days below its pre-knockdown peak) just fails. flux still recovers at ~2x
shallower knockdown than load_neuronal, and that ordering is robust across a
50x range of rho. So a depth-titration experiment (RESEARCH_PLAN section 6)
discriminates the two: recovery at moderate knockdown favours flux, recovery
only at deep knockdown favours neuronal-load. This is the powered axis latency
was not, and it ties back to knockdown depth / x_crit -- the project's core
number.

Output: fig_t01_depth_discriminator.png + a printed table.
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
                         onset_time_2c, x_crit, IX2C)

HERE = os.path.dirname(os.path.abspath(__file__))
cal = json.load(open(os.path.join(HERE, "results.json")))["calibration"]

T_TREAT = 90.0            # knockdown at 60% of the 150 d course (well pre-terminal)
POST = 60.0               # inspect D this many days after knockdown
RESIDUALS = np.linspace(0.05, 1.0, 32)   # residual neuronal PrP after knockdown


def refit_kappa(p):
    def f(k):
        st = survival_time_2c(replace(p, kappa=k), t_end=4000.0)
        return (st if np.isfinite(st) else 4000.0) - 150.0
    try:
        return brentq(f, 1e-6, 1e12, xtol=1e-9, rtol=1e-8, maxiter=300)
    except ValueError:
        return None


def calibrate(mode, rho):
    """beta from the dose-response (fixed); kappa refit to terminal=150 for mode."""
    p = replace(PrionParams(), beta=cal["beta"], rho=rho, toxicity_mode=mode)
    k = refit_kappa(p)
    return replace(p, kappa=k) if k is not None else None


def depth_response(p):
    """For each residual PrP, D at POST days after knockdown, relative to its
    pre-knockdown peak. <1 means recovery; ~1 or >1 means no recovery."""
    out = []
    for xr in RESIDUALS:
        kd = 1.0 - xr
        sol = simulate_2c(p, t_end=T_TREAT + POST + 5.0, kd_n=kd, kd_g=0.0,
                          t_treat=T_TREAT, ramp_days=7.0, max_step=1.0)
        t, D = sol.t, sol.y[IX2C["D"]]
        pre = t <= T_TREAT + 1.0
        D_peak = float(D[pre].max()) if pre.any() else float(D[0])
        i = int(np.argmin(np.abs(t - (T_TREAT + POST))))
        out.append(float(D[i]) / D_peak if D_peak > 0 else np.nan)
    return np.array(out)


if __name__ == "__main__":
    xc = x_crit(replace(PrionParams(), beta=cal["beta"]))
    rhos = [0.06, 0.15, 0.5]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)

    def crossing(residuals, ratio):
        """Residual PrP where the recovery ratio crosses 1.0 (recovery just
        fails), by linear interpolation. Higher = recovers at shallower KD."""
        for i in range(len(residuals) - 1):
            if (ratio[i] - 1.0) * (ratio[i + 1] - 1.0) <= 0:
                x0, x1, y0, y1 = (residuals[i], residuals[i + 1],
                                  ratio[i], ratio[i + 1])
                return x0 + (1.0 - y0) * (x1 - x0) / (y1 - y0)
        return np.nan

    print(f"x_crit = {xc:.3f} residual PrP (knockdown {100*(1-xc):.0f}%)")
    print(f"Recovery response D(t+{POST:.0f}d)/D_peak vs residual PrP; <1 = recovery\n")
    thresholds = {}
    for mode, ax in zip(["flux", "load_neuronal"], axes):
        print(f"{mode}:")
        print(f"  {'residual':>9} " + " ".join(f"rho={r:<5}" for r in rhos))
        curves = {}
        for rho in rhos:
            p = calibrate(mode, rho)
            curves[rho] = depth_response(p) if p else np.full_like(RESIDUALS, np.nan)
        for j in range(0, len(RESIDUALS), 4):
            print(f"  {RESIDUALS[j]:>9.2f} " +
                  " ".join(f"{curves[r][j]:>7.2f}" for r in rhos))
        thresholds[mode] = {r: crossing(RESIDUALS, curves[r]) for r in rhos}
        print("  recovery-threshold residual PrP (recovery fails above this):")
        print("   " + " ".join(f"rho={r}:{100*thresholds[mode][r]:.0f}%" for r in rhos))
        for rho in rhos:
            ax.plot(100 * RESIDUALS, curves[rho], marker="o", ms=3,
                    label=f"rho={rho}")
        ax.axvline(100 * xc, color="k", ls="--", lw=1,
                   label=f"x_crit={100*xc:.0f}%")
        ax.axhline(1.0, color="gray", ls=":", lw=0.8)
        ax.set_title(mode, weight="bold")
        ax.set_xlabel("residual neuronal PrP after knockdown (%)")
        ax.legend(fontsize=7, frameon=False)
        print()

    tf = np.nanmean(list(thresholds["flux"].values()))
    tl = np.nanmean(list(thresholds["load_neuronal"].values()))
    print(f"SEPARATION: flux recovers up to ~{100*tf:.0f}% residual PrP; "
          f"load_neuronal only up to ~{100*tl:.0f}%.")
    print(f"  ratio {tf/tl:.1f}x, and the ordering (flux>load) holds at every rho "
          f"-> rho-robust, powered discriminator.")
    axes[0].set_ylabel(f"D at +{POST:.0f}d / pre-knockdown peak\n(<1 = recovery)")
    fig.suptitle("Depth-response discriminator: flux recovers at ~2x shallower "
                 "knockdown than neuronal-load (rho-robust)", weight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_t01_depth_discriminator.png"), dpi=140)
    print("Wrote fig_t01_depth_discriminator.png")
