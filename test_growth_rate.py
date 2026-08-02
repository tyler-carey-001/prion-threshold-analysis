"""
test_growth_rate.py  --  V1
===========================
Independent validation of the hand-derived determinant condition behind
``growth_rate()`` and ``x_crit()``.

WHY THIS EXISTS
---------------
``growth_rate()`` and ``x_crit()`` are both closed form and both descend from
the same hand-derived 2x2 determinant condition. Comparing them to each other
therefore proves nothing: an algebra slip in that derivation is invisible to
any check built from the same two functions. And the whole "65-90% knockdown"
headline rests on it.

The only independent reference is the ODE itself. So: integrate the actual
``_rhs`` used by the model, measure the empirical exponential growth rate of
``z``, and compare against the analytic ``r(x)``.

THREE THINGS THAT WOULD OTHERWISE PRODUCE FALSE ALARMS
------------------------------------------------------
1. ``x`` must be CLAMPED. "Integrate the full nonlinear system" and "hold x
   fixed" are contradictory instructions: left free, ``x`` depletes as prion
   load grows, the empirical growth rate falls below analytic ``r(x)``, and
   substrate depletion gets misread as an algebra error. We call the real
   ``_rhs`` and zero its dx component, so the (y, z) code path under test is
   the production one.

2. The early transient must be DISCARDED. The inoculum
   ``[y, z] = [y0, n*y0]`` is not the dominant eigenvector, so the trajectory
   is a mix of both eigenmodes. The subdominant mode decays relative to the
   dominant one at rate ``sqrt(disc) = r1 - r2``; we wait 25 e-foldings of
   that separation before fitting, and we fit a regression over a window and
   report its residual rather than taking a two-point slope.

3. Tolerance must be REGIME-DEPENDENT. Near ``x_crit`` the true rate passes
   through zero, so a relative tolerance is arithmetically impossible to
   satisfy there. Absolute tolerance near threshold, relative away from it.

Run:  python test_growth_rate.py     (or: pytest test_growth_rate.py)
"""

import numpy as np
from scipy.integrate import solve_ivp

from prion_model import PrionParams, growth_rate, x_crit, _rhs

# --- tolerances -------------------------------------------------------------
RTOL_RATE = 2e-6      # relative, away from threshold
ATOL_RATE = 2e-8      # absolute, near threshold (per day)
NEAR_THRESHOLD = 1e-3  # |r| below this -> use ATOL_RATE
MAX_LOG_RESID = 1e-6  # RMSE of log(z) about the linear fit


def _clamped_rhs(t, s, p, kd_fn):
    """The production `_rhs`, with dx/dt (and the toxicity layer) forced to zero.

    This is the point of the test: exercise the same (y, z) equations the
    model actually integrates, with substrate depletion removed so that the
    linear-regime growth rate is well defined.

    D and N are frozen as well. They cannot influence (y, z) -- the coupling is
    one-way -- but leaving them live makes the integrator fail spuriously: at
    deeply subcritical x, z decays to ~1e-13, so dD = kappa*conv_flux is a tiny
    nonzero derivative on a component whose value is exactly 0. The error scale
    for that component is then atol alone, the local error ratio explodes, and
    the step size collapses. Freezing them keeps the test about the algebra.
    """
    ds = _rhs(t, s, p, kd_fn)
    ds[0] = 0.0
    ds[3] = 0.0
    ds[4] = 0.0
    return ds


def empirical_growth_rate(x, p, n_eft_transient=25.0, n_eft_fit=6.0,
                          max_fit_days=4000.0, n_samples=400):
    """Measure dz/dt's exponential rate by integrating the ODE at fixed x.

    Returns (rate, log_residual_rmse, window).
    """
    # eigenvalue separation sets how fast the subdominant mode dies out
    A = p.a + p.b * (2 * p.n - 1)
    C = p.beta * x - p.b * p.n * (p.n - 1)
    disc = (A - p.a) ** 2 + 4.0 * p.b * C
    sep = np.sqrt(disc) if disc > 0 else np.nan
    if not np.isfinite(sep) or sep <= 0:
        raise ValueError(f"complex/degenerate eigenvalues at x={x}")

    t_transient = n_eft_transient / sep

    r_analytic = growth_rate(x, p)
    if abs(r_analytic) > NEAR_THRESHOLD:
        fit_len = min(n_eft_fit / abs(r_analytic), max_fit_days)
    else:
        fit_len = max_fit_days          # r ~ 0: no e-folding to wait for

    t0, t1 = t_transient, t_transient + fit_len
    t_eval = np.linspace(t0, t1, n_samples)

    s0 = [x, p.y_inoculum, p.y_inoculum * p.n, 0.0, 0.0]
    # per-component atol: y and z span many orders of magnitude and must be
    # governed by rtol, while the frozen components keep a sane error scale
    atol = np.array([1e-14, 1e-40, 1e-40, 1e-14, 1e-14])
    sol = solve_ivp(_clamped_rhs, (0.0, t1), s0, args=(p, lambda t: 0.0),
                    method="DOP853", rtol=1e-12, atol=atol,
                    t_eval=t_eval, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"integration failed at x={x}: {sol.message}")

    z = sol.y[2]
    if np.any(z <= 0):
        raise RuntimeError(f"non-positive z at x={x}; cannot take log")

    # regression, not a two-point slope; residual is reported
    logz = np.log(z)
    slope, intercept = np.polyfit(sol.t, logz, 1)
    resid = logz - (slope * sol.t + intercept)
    rmse = float(np.sqrt(np.mean(resid ** 2)))

    # x must actually have stayed put
    assert np.allclose(sol.y[0], x, rtol=0, atol=1e-12), "x drifted despite clamp"

    return float(slope), rmse, (t0, t1)


def check(p=None, verbose=True):
    p = p or PrionParams()
    xc = x_crit(p)

    # bracket x_crit on both sides, including very close to it
    fracs = [0.30, 0.70, 0.95, 0.99, 1.00, 1.01, 1.05, 1.40, 2.00, 2.70]
    xs = [f * xc for f in fracs]

    failures = []
    rows = []
    for x in xs:
        r_an = growth_rate(x, p)
        r_em, rmse, win = empirical_growth_rate(x, p)

        near = abs(r_an) < NEAR_THRESHOLD
        if near:
            ok_rate = abs(r_em - r_an) < ATOL_RATE
            tol_desc = f"abs<{ATOL_RATE:g}"
        else:
            ok_rate = abs(r_em - r_an) <= RTOL_RATE * abs(r_an)
            tol_desc = f"rel<{RTOL_RATE:g}"

        ok_resid = rmse < MAX_LOG_RESID
        # minimum bar: the sign must agree either side of the threshold
        ok_sign = (np.sign(r_em) == np.sign(r_an)) or near

        rows.append((x / xc, x, r_an, r_em, r_em - r_an, rmse,
                     tol_desc, ok_rate and ok_resid and ok_sign))
        if not (ok_rate and ok_resid and ok_sign):
            failures.append((x, r_an, r_em, rmse, ok_rate, ok_resid, ok_sign))

    if verbose:
        print("=" * 96)
        print("V1  determinant condition vs the ODE  (x clamped, transient discarded)")
        print("=" * 96)
        print(f"  params: beta={p.beta} b={p.b} a={p.a} n={p.n}")
        print(f"  x_crit = {xc:.10f}   (required knockdown {100*(1-xc/p.x0):.2f}%)")
        print()
        print(f"{'x/x_crit':>9} {'x':>10} {'r_analytic':>14} {'r_empirical':>14} "
              f"{'diff':>12} {'log resid':>11} {'tol':>13}  ok")
        for f, x, ra, re_, d, rm, td, ok in rows:
            print(f"{f:>9.3f} {x:>10.5f} {ra:>14.9f} {re_:>14.9f} "
                  f"{d:>12.2e} {rm:>11.2e} {td:>13}  {'PASS' if ok else 'FAIL'}")
        print()

    # r(x_crit) must itself be ~0 -- the definition of the threshold
    r_at_xc = growth_rate(xc, p)
    if abs(r_at_xc) > 1e-12:
        failures.append(("x_crit self-consistency", r_at_xc, 0.0, 0.0,
                         False, True, True))
        if verbose:
            print(f"  FAIL: growth_rate(x_crit) = {r_at_xc:.3e}, expected ~0")
    elif verbose:
        print(f"  growth_rate(x_crit) = {r_at_xc:.3e}  (expected ~0)  PASS")

    if verbose:
        print()
        if failures:
            print(f"  RESULT: {len(failures)} FAILURE(S). The determinant condition "
                  f"in growth_rate()/x_crit() does not match the ODE.")
            print("  Everything downstream of x_crit is suspect until this is fixed.")
        else:
            print("  RESULT: PASS. The closed-form growth rate matches the integrated")
            print("  ODE across the threshold, so the hand-derived determinant holds.")
        print("=" * 96)

    return failures


def test_growth_rate_matches_ode():
    """pytest entry point."""
    failures = check(verbose=False)
    assert not failures, f"{len(failures)} growth-rate mismatches: {failures}"


if __name__ == "__main__":
    import sys
    sys.exit(1 if check() else 0)
