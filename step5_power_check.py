"""
step5_power_check.py
====================
Does the held-out latency test actually have power to discriminate the toxicity
modes? Committed BEFORE the Mallucci 2007 recovery timecourse is fetched, so the
answer cannot be tuned to the observed value.

A held-out test discriminates only if the training data pinned the prediction
before the held-out value was seen. The prediction here is the 50%-recovery
latency of dysfunction D after neuron-specific PrP knockdown. We ask, per mode:
given everything training can constrain, how wide is the latency prediction
interval? A wide interval means the mode can accommodate almost any observation
post-hoc, and "mode X matches Mallucci 2007" is a second fitting exercise
wearing a prediction costume.

TRAINING DATA (same as the original calibration):
  * beta pinned by the ANCHORS survival dose-response (independent of the
    toxicity layer -- relative_survival_analytic uses only growth_rate);
  * kappa refit by root-finding so the untreated animal reaches terminal at
    150 d, exactly as run_analysis.py does;
  * the emergent untreated onset day, which must land in the natural-history
    band for the early behavioural (burrowing) deficit.

Everything else about the toxicity layer -- crucially rho (D's repair rate) and,
for the load modes, a (polymer clearance) -- is what training must pin if the
latency prediction is to be sharp.
"""

import json
import os
import numpy as np
from dataclasses import replace
from scipy.optimize import brentq

from prion_model import (PrionParams, simulate, survival_time, onset_time,
                         simulate_2c, survival_time_2c, onset_time_2c,
                         solve_beta_for_target, x_crit, IX2C, z_total_2c)

HERE = os.path.dirname(os.path.abspath(__file__))
cal = json.load(open(os.path.join(HERE, "results.json")))["calibration"]

# knockdown protocol used to elicit recovery (fixed; not a fitted quantity)
KD_N = 0.95
T_TREAT = 105.0          # ~70% of the 150 d untreated course
TERMINAL_TARGET = 150.0

# natural-history onset band for the early behavioural deficit.
#   LOOSE: what a lumped toy model calibrated to a different mouse line can
#          honestly claim (+/-15 d around the model's emergent 76 d).
#   TIGHT: an optimistic +/-3 d, to show the interval even in the best case.
ONSET_LOOSE = (76.0 - 15.0, 76.0 + 15.0)
ONSET_TIGHT = (76.0 - 3.0, 76.0 + 3.0)


def refit_kappa(p):
    """kappa such that the untreated animal terminates at 150 d.

    Uses the TWO-COMPARTMENT untreated course with p's own toxicity_mode. This
    matters: the single-compartment _rhs hardcodes the flux toxicity form and
    ignores toxicity_mode, so calibrating a load mode against it fits kappa to
    the wrong dynamics (flux ~ beta*x*y vs load ~ z, orders of magnitude apart)
    and the animal dies almost immediately. For flux this agrees with the 1c
    calibration by the reduction invariant.
    """
    def f(k):
        st = survival_time_2c(replace(p, kappa=k), t_end=4000.0)
        return (st if np.isfinite(st) else 4000.0) - TERMINAL_TARGET
    try:
        return brentq(f, 1e-6, 1e12, xtol=1e-9, rtol=1e-8, maxiter=300)
    except ValueError:
        return None


def recovery_latency(p, mode, sigma=0.0):
    """Response to neuron-specific knockdown at T_TREAT.

    Returns (latency_d, died, D_drop_frac):
      died         animal reached terminal (N=N_death) despite knockdown
      latency_d    days to 50% recovery of D (inf if it never halves)
      D_drop_frac  fractional fall of D from its pre-knockdown peak by the end
    'inf latency' therefore has two distinct meanings, now separated: the
    animal DIED (no rescue) vs it survived but D never halved (weak recovery).
    """
    q = replace(p, toxicity_mode=mode, sigma=sigma)
    sol = simulate_2c(q, t_end=T_TREAT + 600.0, kd_n=KD_N, kd_g=0.0,
                      t_treat=T_TREAT, ramp_days=7.0, max_step=1.0)
    t, D = sol.t, sol.y[IX2C["D"]]
    died = bool(sol.t_events and len(sol.t_events[0]) > 0)
    pre = t <= T_TREAT + 1.0
    D_peak = float(D[pre].max()) if pre.any() else float(D.max())
    D_end = float(D[-1])
    drop = (D_peak - D_end) / D_peak if D_peak > 0 else np.nan
    half = 0.5 * D_peak
    after = np.where((t > T_TREAT) & (D <= half))[0]
    lat = float(t[after[0]] - T_TREAT) if len(after) else np.inf
    return lat, died, drop


def flux_profile():
    """Profile rho (kappa refit each time); record onset and flux latency."""
    base = replace(PrionParams(), beta=cal["beta"], toxicity_mode="flux")
    rows = []
    for rho in np.geomspace(0.01, 0.5, 25):
        p = replace(base, rho=rho)
        k = refit_kappa(p)
        if k is None:
            continue
        p = replace(p, kappa=k)
        onset = onset_time_2c(simulate_2c(p, t_end=400.0), p)
        lat, died, drop = recovery_latency(p, "flux")
        rows.append((rho, k, k / rho, onset, lat, died, drop))
    return rows


def load_neuronal_profile():
    """Profile a (polymer clearance) and rho for the neuronal-load mode.

    a is FREE (pinning it would rig this mode's latency -- see plan item 1).
    Freeing a moves x_crit, so beta is refit to the dose-response at each a,
    then kappa to terminal. Latency measured with sigma=0 (best case for
    recovery; nonzero sigma props z_n up via glial reseeding and only lengthens
    or abolishes recovery)."""
    rows = []
    for a in [0.01, 0.02, 0.04, 0.08]:
        beta = solve_beta_for_target(replace(PrionParams(), a=a), 0.50, 2.70,
                                     t_tox_frac=0.25)
        if beta is None:
            continue
        for rho in [0.02, 0.06, 0.15]:
            p = replace(PrionParams(), a=a, beta=beta, rho=rho,
                        toxicity_mode="load_neuronal")
            k = refit_kappa(p)
            if k is None:
                continue
            p = replace(p, kappa=k)
            onset = onset_time_2c(simulate_2c(p, t_end=600.0), p)
            lat, died, drop = recovery_latency(p, "load_neuronal")
            rows.append((a, rho, x_crit(p), onset, lat, died, drop))
    return rows


def load_total_reverses():
    """Can load_total ever satisfy the reversal criterion (D falls >50% while
    total burden rises)? Scan rho with kappa refit."""
    base = replace(PrionParams(), beta=cal["beta"])
    any_reversal = False
    details = []
    for rho in np.geomspace(0.01, 1.0, 15):
        p = replace(base, rho=rho, toxicity_mode="load_total")
        k = refit_kappa(p)
        if k is None:
            continue
        p = replace(p, kappa=k)
        sol = simulate_2c(p, t_end=T_TREAT + 200.0, kd_n=KD_N, kd_g=0.0,
                          t_treat=T_TREAT, max_step=1.0)
        t, D = sol.t, sol.y[IX2C["D"]]
        zt = z_total_2c(sol)
        pre = t <= T_TREAT + 1.0
        i_end = np.argmin(np.abs(t - (T_TREAT + 40.0)))
        i_tr = np.argmin(np.abs(t - T_TREAT))
        D_drop = (D[pre].max() - D[i_end]) / D[pre].max()
        zt_rises = zt[i_end] > zt[i_tr]
        reversed_ = D_drop > 0.5 and zt_rises
        any_reversal = any_reversal or reversed_
        details.append((rho, D_drop, zt_rises, reversed_))
    return any_reversal, details


def interval(latencies):
    lat = [l for l in latencies if np.isfinite(l)]
    return (min(lat), max(lat)) if lat else (np.inf, np.inf)


if __name__ == "__main__":
    print("=" * 92)
    print("STEP 5 POWER CHECK -- committed before the Mallucci 2007 latency is fetched")
    print("=" * 92)

    # ---- flux ----
    fr = flux_profile()
    print("\nFLUX MODE  (profile rho; kappa refit to terminal=150 in 2c)")
    print(f"{'rho':>7} {'kappa':>11} {'kappa/rho':>11} {'onset_d':>8} {'latency_d':>10} {'died':>5}")
    for rho, k, kr, onset, lat, died, drop in fr:
        print(f"{rho:>7.3f} {k:>11.1f} {kr:>11.1f} {onset:>8.2f} {lat:>10.2f} {str(died):>5}")

    kr_vals = [r[2] for r in fr]
    print(f"\n  kappa/rho spread: {min(kr_vals):.1f}-{max(kr_vals):.1f} "
          f"(x{max(kr_vals)/min(kr_vals):.2f}) -- the degeneracy: onset fixes "
          f"the COMBINATION kappa/rho, not rho.")

    for name, band in (("loose +/-15d", ONSET_LOOSE), ("tight +/-3d", ONSET_TIGHT)):
        keep = [r for r in fr if band[0] <= r[3] <= band[1] and not r[5]]
        lo, hi = interval([r[4] for r in keep])
        sp = f"x{hi/lo:.1f}" if np.isfinite(hi) and lo > 0 else "n/a"
        print(f"  onset anchor {name}: {len(keep)} rho values admitted (survived), "
              f"flux latency interval = [{lo:.1f}, {hi:.1f}] d ({sp} spread)")

    # ---- load_neuronal ----
    ln = load_neuronal_profile()
    print("\nLOAD_NEURONAL MODE  (a FREE; beta refit per a, kappa refit to terminal in 2c)")
    print(f"{'a':>6} {'rho':>6} {'x_crit':>8} {'onset_d':>8} {'latency_d':>10} {'died':>5} {'Ddrop':>7}")
    for a, rho, xc, onset, lat, died, drop in ln:
        print(f"{a:>6.3f} {rho:>6.3f} {xc:>8.3f} {onset:>8.1f} {lat:>10.2f} "
              f"{str(died):>5} {100*drop:>6.0f}%")
    surv = [r for r in ln if not r[5]]
    lo, hi = interval([r[4] for r in surv])
    n_died = sum(1 for r in ln if r[5])
    print(f"  of {len(ln)} (a,rho) combos: {n_died} DIED (no late rescue), "
          f"{len(surv)} survived; survivors' latency interval [{lo:.1f}, {hi:.1f}] d")

    # ---- load_total ----
    rev, det = load_total_reverses()
    print("\nLOAD_TOTAL MODE  (can it reverse at all?)")
    print(f"  any rho gives D-drop>50% while total burden rises: {rev}")

    # ---- verdict ----
    flux_surv = [r for r in fr if ONSET_LOOSE[0] <= r[3] <= ONSET_LOOSE[1]
                 and not r[5]]
    lo_f, hi_f = interval([r[4] for r in flux_surv])
    ln_surv = [r for r in ln if not r[5]]
    lo_n, hi_n = interval([r[4] for r in ln_surv])
    overlap = (np.isfinite(hi_f) and np.isfinite(hi_n)
               and max(lo_f, lo_n) <= min(hi_f, hi_n))
    ln_all_die = len(ln_surv) == 0
    print("\n" + "=" * 92)
    print("VERDICT")
    print(f"  flux (survivors):          latency [{lo_f:.1f}, {hi_f:.1f}] d, "
          f"{len(flux_surv)}/{len(fr)} survive late knockdown")
    print(f"  load_neuronal (survivors): latency [{lo_n:.1f}, {hi_n:.1f}] d, "
          f"{len(ln_surv)}/{len(ln)} survive late knockdown")
    print(f"  latency intervals overlap: {overlap}")
    print("  Two separate discrimination axes exist -- read both, do not collapse them:")
    print("   (1) LATENCY: flux's interval alone spans an order of magnitude because")
    print("       training fixes only kappa/rho, not rho. On latency alone the test is")
    print("       UNDERPOWERED regardless of what load_neuronal does.")
    print("   (2) RESCUE-AT-LATE-KNOCKDOWN: whether a mode rescues at 70% of course is")
    print("       itself mode- and parameter-dependent. If load_neuronal universally")
    print("       fails to rescue late while flux can, THAT is the real discriminator,")
    print("       not latency -- and it is a survival readout, not a timecourse.")
    print("=" * 92)
