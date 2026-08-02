"""
expressibility_check.py
=======================
QUALITATIVE check, Step 3 of the plan. NO acceptance threshold is derived from
this script -- it exists only to confirm that the two-compartment structure can
*express* Mallucci's observation at all, which the single-compartment model
provably cannot (there, knockdown drops conversion flux everywhere and standing
load z always falls -- every treated trajectory in fig4 bends downward).

The claim under test here is narrow and about ONE mode:

    Under the flux model, neuron-specific knockdown (kd_n high, kd_g = 0) should
    make dysfunction D fall while total PrP-Sc load z_total = z_n + z_g RISES --
    exactly Mallucci's dissociation of toxicity from burden.

The three-way discrimination (flux vs load_total vs load_neuronal) is NOT done
here. It requires each mode to carry its own kappa, refit so that the untreated
animal still dies at ~150 d; with the flux-calibrated kappa the load models put
D far above threshold and the animal is dead long before knockdown. That fair
refit and comparison happen in Step 8 (t01_model_comparison.py), after the
numeric criteria are pre-registered and committed.
"""

import json
import os
import numpy as np
from dataclasses import replace

from prion_model import PrionParams, simulate_2c, z_total_2c, IX2C

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "results.json")) as fh:
    cal = json.load(fh)["calibration"]

TERMINAL = cal["untreated_terminal_d"]           # ~150 d
T_TREAT = 0.70 * TERMINAL                         # neuron-specific knockdown onset
WINDOW = 40.0                                     # days after onset to inspect


def flux_case(kd_n=0.9, sigma=0.0):
    p = replace(PrionParams(), beta=cal["beta"], kappa=cal["kappa"],
                toxicity_mode="flux", sigma=sigma)
    sol = simulate_2c(p, t_end=T_TREAT + WINDOW + 5.0, kd_n=kd_n, kd_g=0.0,
                      t_treat=T_TREAT, seed_both=True, max_step=1.0)
    t = sol.t
    D = sol.y[IX2C["D"]]
    zt = z_total_2c(sol)
    zn = sol.y[IX2C["z_n"]]
    zg = sol.y[IX2C["z_g"]]

    i_treat = int(np.argmin(np.abs(t - T_TREAT)))
    i_after = int(np.argmin(np.abs(t - (T_TREAT + WINDOW))))
    D_peak = float(np.max(D[t <= T_TREAT + 1.0]))
    return {
        "sigma": sigma,
        "D_peak": D_peak, "D_after": float(D[i_after]),
        "zt": (float(zt[i_treat]), float(zt[i_after])),
        "zn": (float(zn[i_treat]), float(zn[i_after])),
        "zg": (float(zg[i_treat]), float(zg[i_after])),
        "terminated": bool(sol.t_events and len(sol.t_events[0]) > 0),
    }


if __name__ == "__main__":
    print("=" * 84)
    print("Expressibility (QUALITATIVE -- flux mode only, no thresholds frozen)")
    print(f"neuron-specific knockdown kd_n=0.9, kd_g=0, onset day {T_TREAT:.0f} "
          f"(70% of {TERMINAL:.0f} d), inspected {WINDOW:.0f} d later")
    print("=" * 84)
    all_ok = True
    for sigma in (0.0, 0.02, 0.1):
        r = flux_case(sigma=sigma)
        D_drop = 100 * (r["D_peak"] - r["D_after"]) / r["D_peak"]
        zt_rises = r["zt"][1] > r["zt"][0]
        D_falls = r["D_after"] < r["D_peak"]
        ok = D_falls and zt_rises and not r["terminated"]
        all_ok = all_ok and ok
        print(f"\nsigma = {sigma}")
        print(f"  D        {r['D_peak']:.4f} -> {r['D_after']:.4f}   "
              f"(drop {D_drop:+.0f}% over {WINDOW:.0f} d)   "
              f"{'FALLS' if D_falls else 'does not fall'}")
        print(f"  z_total  {r['zt'][0]:.4e} -> {r['zt'][1]:.4e}   "
              f"{'RISES' if zt_rises else 'falls'}")
        print(f"    z_n    {r['zn'][0]:.4e} -> {r['zn'][1]:.4e}   (neuronal, knocked down)")
        print(f"    z_g    {r['zg'][0]:.4e} -> {r['zg'][1]:.4e}   (glial, untreated)")
        print(f"  -> {'EXPRESSIBLE' if ok else 'NOT expressible'}: "
              f"D falls while total burden rises")
    print()
    print("RESULT:", "PASS -- the two-compartment flux model expresses Mallucci's"
          " dissociation." if all_ok else "FAIL -- structure still cannot express it.")
