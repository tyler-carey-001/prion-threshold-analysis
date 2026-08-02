"""
test_compartment_reduction.py
=============================
The design invariant for the two-compartment extension: in the symmetric limit
it must reproduce the validated single-compartment model *exactly*.

Symmetric limit = identical compartments, identical (global) knockdown,
toxicity_mode="flux", inoculum seeded into both compartments. Then:

  * the two compartments evolve identically, so the diffusive coupling
    sigma*(other - self) is identically zero and the result is independent of
    sigma (we deliberately use sigma != 0 to prove this);
  * neuronal conversion flux equals the whole-brain flux of the 1c model, so
    the flux-driven D and the neuron loss N evolve identically;
  * onset day and terminal day therefore match to integrator tolerance.

If this passes, any later divergence between the models is attributable to the
compartment structure or the toxicity mode, never to a transcription error in
the replication equations that both are supposed to share.

Run:  python test_compartment_reduction.py   (or: pytest)
"""

import json
import os

import numpy as np
from dataclasses import replace

from prion_model import (PrionParams, simulate, simulate_2c,
                         survival_time, survival_time_2c,
                         onset_time, onset_time_2c, IX2C)

HERE = os.path.dirname(os.path.abspath(__file__))
RTOL_TRAJ = 1e-6
ATOL_TRAJ = 1e-9
RTOL_DAY = 1e-4      # onset/terminal day agreement (fraction)


def _calibrated_params():
    """The calibrated (beta, kappa) from the frozen regression baseline."""
    with open(os.path.join(HERE, "results.json")) as fh:
        cal = json.load(fh)["calibration"]
    # sigma != 0 on purpose: it must not matter in the symmetric limit
    return replace(PrionParams(), beta=cal["beta"], kappa=cal["kappa"],
                   toxicity_mode="flux", sigma=0.037)


def _compare_case(p, kd, t_treat, t_end, verbose):
    """Run 1c and 2c under a matched (global) knockdown and compare."""
    sol1 = simulate(p, t_end=t_end, knockdown=kd, t_treat=t_treat)
    sol2 = simulate_2c(p, t_end=t_end, kd_n=kd, kd_g=kd, t_treat=t_treat,
                       seed_both=True)

    # common grid inside both integration spans
    t_lo = 0.0
    t_hi = min(sol1.t[-1], sol2.t[-1])
    tg = np.linspace(t_lo, t_hi, 500)
    a = sol1.sol(tg)            # [x, y, z, D, N]
    b = sol2.sol(tg)            # [x_n, x_g, y_n, z_n, y_g, z_g, D, N]

    # 1. compartments stayed symmetric (coupling vanished despite sigma != 0)
    sym_x = np.max(np.abs(b[IX2C["x_n"]] - b[IX2C["x_g"]]))
    sym_y = np.max(np.abs(b[IX2C["y_n"]] - b[IX2C["y_g"]]))
    sym_z = np.max(np.abs(b[IX2C["z_n"]] - b[IX2C["z_g"]]))
    sym = max(sym_x, sym_y, sym_z)

    # 2. per-variable agreement 1c vs neuronal compartment of 2c
    pairs = {
        "x": (a[0], b[IX2C["x_n"]]),
        "y": (a[1], b[IX2C["y_n"]]),
        "z": (a[2], b[IX2C["z_n"]]),
        "D": (a[3], b[IX2C["D"]]),
        "N": (a[4], b[IX2C["N"]]),
    }
    worst = {}
    ok = True
    for name, (u, v) in pairs.items():
        err = np.abs(u - v)
        tol = ATOL_TRAJ + RTOL_TRAJ * np.abs(u)
        worst[name] = float(np.max(err / (tol)))    # >1 means out of tolerance
        ok = ok and np.all(err <= tol)

    # 3. scalar onset / terminal days
    o1, o2 = onset_time(sol1, p), onset_time_2c(sol2, p)
    s1 = survival_time(p, t_end=t_end, knockdown=kd, t_treat=t_treat)
    s2 = survival_time_2c(p, t_end=t_end, kd_n=kd, kd_g=kd, t_treat=t_treat,
                          seed_both=True)
    day_ok = True
    for u, v in ((o1, o2), (s1, s2)):
        if np.isfinite(u) and np.isfinite(v):
            day_ok = day_ok and abs(u - v) <= RTOL_DAY * max(u, 1.0)
        else:
            day_ok = day_ok and (np.isfinite(u) == np.isfinite(v))

    sym_ok = sym <= 1e-9
    passed = ok and day_ok and sym_ok

    if verbose:
        label = f"kd={kd:.2f}@d{t_treat:.0f}"
        print(f"  {label:<16} symmetry={sym:.1e}  "
              f"traj worst-ratio={max(worst.values()):.2e}  "
              f"onset {o1:.3f}/{o2:.3f}  terminal {s1:.3f}/{s2:.3f}  "
              f"{'PASS' if passed else 'FAIL'}")
    return passed


def check(verbose=True):
    p = _calibrated_params()
    if verbose:
        print("=" * 88)
        print("Reduction invariant: 2c symmetric limit == validated 1c model")
        print("=" * 88)
        print(f"  beta={p.beta:.6f} kappa={p.kappa:.3f} sigma={p.sigma} "
              f"(sigma must not matter here)")
    cases = [
        (0.00, 0.0, 400.0),     # untreated course (onset ~76, terminal ~150)
        (0.85, 100.0, 400.0),   # late global knockdown (exercises the ramp)
        (0.50, 40.0, 400.0),    # partial early global knockdown
    ]
    results = [_compare_case(p, kd, tt, te, verbose) for kd, tt, te in cases]
    passed = all(results)
    if verbose:
        print()
        print("  RESULT:", "PASS -- the 2c model reduces exactly to the 1c model."
              if passed else "FAIL -- 2c does not reduce to 1c; replication code diverges.")
        print("=" * 88)
    return passed


def test_reduction():
    assert check(verbose=False), "2c model does not reduce to 1c in symmetric limit"


if __name__ == "__main__":
    import sys
    sys.exit(0 if check() else 1)
