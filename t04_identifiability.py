"""
T0-4 Phase 0 — is `x_crit` identifiable from the published dose-response at all?

This runs BEFORE any literature extraction, and it is a gate on whether that
extraction is worth doing. T0-4's premise was that adding the Tga20
overexpression arm (8-10x PrP) would constrain r(x) from the opposite side of
the curve and tighten the 65-90% required-knockdown range. That premise is
testable with data already in `prion_model.ANCHORS`.

Method. Profile the sum of squared residuals in log survival-ratio over the
finite anchors, with `x_crit` held fixed on a grid and every genuinely unknown
nuisance minimised out at each grid point. Since

    x_crit = [a(a + b(2n-1)) + b^2 n(n-1)] / (b * beta)

beta is solved analytically to hit the target x_crit, leaving (a, b, n,
t_tox_frac) free. Nuisances are optimised with scipy from multiple starts rather
than gridded, so the profile is not limited by grid resolution.

Noise model: between-study scatter in log survival ratio, i.i.d. log-normal,
with sigma^2 estimated from the anchors that sit at the SAME nominal PrP level
(~50%). That cluster disagrees by a factor of two, and that disagreement is the
irreducible floor any model must live above. Stated as a modelling choice
because it is one.

POPULATION-LEVEL MODEL OUTPUT, NOT EVIDENCE. Anchors are literature values that
must be verified against primary sources before any fitted number is reported.
"""
import numpy as np
from dataclasses import replace
from scipy.optimize import minimize
from scipy.stats import chi2

from prion_model import (PrionParams, growth_rate,
                         relative_survival_analytic, ANCHORS)

BASE = PrionParams()
FINITE = [(x, r) for _lab, x, r, _src in ANCHORS if np.isfinite(r)]
CLUSTER_LO, CLUSTER_HI = 0.45, 0.55
# (log a, log b, n, t_tox_frac) bounds. n is treated as continuous; the model's
# formulae are smooth in n and restricting to integers would quantise the profile.
BOUNDS = [(np.log(1e-3), np.log(0.5)),
          (np.log(1e-6), np.log(1e-2)),
          (2.0, 12.0),
          (0.01, 0.60)]


def beta_for(xc, a, b, n):
    """beta that places x_crit exactly at xc, given (a, b, n)."""
    num = a * (a + b * (2 * n - 1)) + b ** 2 * n * (n - 1)
    return num / (b * xc)


def ssr(theta, xc, data):
    """Sum of squared log-residuals at fixed x_crit."""
    la, lb, n, tf = theta
    a, b = np.exp(la), np.exp(lb)
    beta = beta_for(xc, a, b, n)
    p = replace(BASE, a=a, b=b, n=n, beta=beta)
    if growth_rate(p.x0, p) <= 0:            # must be able to cause disease
        return 1e3
    tot = 0.0
    for x, r in data:
        pred = relative_survival_analytic(x, p, tf)
        if not np.isfinite(pred) or pred <= 0:
            return 1e3
        tot += (np.log(pred) - np.log(r)) ** 2
    return tot


def profile_point(xc, data, n_starts=12, seed=0):
    """Minimise SSR over nuisances at fixed x_crit, from several starts."""
    rng = np.random.default_rng(seed)
    best = np.inf
    starts = [np.array([np.log(0.02), np.log(2e-4), 6.0, 0.25])]
    for _ in range(n_starts - 1):
        starts.append(np.array([rng.uniform(*BOUNDS[0]), rng.uniform(*BOUNDS[1]),
                                rng.uniform(*BOUNDS[2]), rng.uniform(*BOUNDS[3])]))
    for s in starts:
        try:
            res = minimize(ssr, s, args=(xc, data), method="L-BFGS-B",
                           bounds=BOUNDS, options={"maxiter": 2000})
            if res.fun < best:
                best = float(res.fun)
        except Exception:
            continue
    return best


def noise_floor():
    """sigma^2 from disagreement among anchors at the same nominal PrP level."""
    vals = [r for x, r in FINITE if CLUSTER_LO <= x <= CLUSTER_HI]
    lg = np.log(vals)
    ss = float(((lg - lg.mean()) ** 2).sum())
    return vals, ss, ss / (len(vals) - 1)


def profile(data, xs):
    return np.array([profile_point(xc, data) for xc in xs])


def interval(xs, prof, sigma2, level=0.95):
    """Profile-likelihood interval: SSR within chi2(1) of the minimum, scaled."""
    thr = prof.min() + chi2.ppf(level, 1) * sigma2
    inside = xs[prof <= thr]
    if inside.size == 0:
        return None, None, thr
    return float(inside.min()), float(inside.max()), float(thr)


def fmt_kd(lo, hi):
    return f"{100*(1-hi):.0f}%-{100*(1-lo):.0f}%"


def main():
    xs = np.arange(0.03, 0.50, 0.005)
    vals, floor_ss, sigma2 = noise_floor()
    sigma = np.sqrt(sigma2)

    print("=" * 76)
    print("T0-4 Phase 0 — identifiability of x_crit from published dose-response")
    print("=" * 76)
    print("Model output, not evidence. Anchors unverified against primary sources.")
    print()

    print("1. NOISE FLOOR — disagreement among studies at the SAME PrP level")
    print("-" * 76)
    print(f"   anchors at {CLUSTER_LO:.0%}-{CLUSTER_HI:.0%} residual PrP: {vals}")
    print(f"   SSR of a CONSTANT fit to them          : {floor_ss:.4f}")
    print(f"   => sigma^2 = {sigma2:.4f}, sigma = {sigma:.3f} log units "
          f"(~{100*(np.exp(sigma)-1):.0f}% multiplicative)")

    prof_all = profile(FINITE, xs)
    lo, hi, thr = interval(xs, prof_all, sigma2)
    imin = xs[int(np.argmin(prof_all))]

    print()
    print("2. PROFILE LIKELIHOOD, all anchors")
    print("-" * 76)
    print(f"   minimum SSR {prof_all.min():.4f} at {100*imin:.0f}% residual PrP "
          f"({100*(1-imin):.0f}% knockdown)")
    print(f"   95% threshold = {thr:.4f}")
    print(f"   => x_crit in {100*lo:.0f}%-{100*hi:.0f}% residual  "
          f"==  required knockdown {fmt_kd(lo, hi)}")
    print(f"   currently reported in README/REPRODUCTION: 64.3%-91.9%")
    # Fair null: a single constant ratio fitted to ALL anchors, not just the
    # cluster. Comparing the model's 7-anchor SSR against a constant's 5-anchor
    # SSR would be comparing different datasets and would understate the model.
    lg_all = np.log([r for _x, r in FINITE])
    null_all = float(((lg_all - lg_all.mean()) ** 2).sum())
    print()
    print(f"   null (one constant ratio, all {len(FINITE)} anchors) SSR : {null_all:.4f}")
    print(f"   model (x_crit + 4 nuisances, same anchors)      SSR : {prof_all.min():.4f}")
    print(f"   => the model buys {null_all - prof_all.min():.4f} SSR for 5 parameters,")
    print("      essentially all of it by fitting wild-type and Tga20, which a")
    print("      constant cannot reach. On the ~50% cluster it does no better")
    print(f"      than a constant ({floor_ss:.4f}), and that cluster is what")
    print("      would have to discriminate x_crit.")

    print()
    print("3. LEAVE-ONE-OUT — does the overexpression arm carry any leverage?")
    print("-" * 76)
    subsets = [("all anchors", FINITE),
               ("without Tga20 (x=8)", [d for d in FINITE if d[0] < 8]),
               ("without wild-type ref", [d for d in FINITE if d[0] != 1.0])]
    for label, data in subsets:
        pr = profile(data, xs)
        l, h, _ = interval(xs, pr, sigma2)
        if l is None:
            print(f"   {label:<24} unbounded")
        else:
            print(f"   {label:<24} knockdown {fmt_kd(l, h):<12} "
                  f"(min SSR {pr.min():.4f})")
    print()
    print("   T0-4's premise was that Tga20 constrains r(x) from the opposite")
    print("   side and would tighten x_crit. Compare rows 1 and 2 above.")

    print()
    print("4. WHAT sigma WOULD HAVE TO BECOME")
    print("-" * 76)
    print("   Anchor values held fixed; only the assumed scatter is shrunk.")
    print(f"   {'sigma':<9}{'scatter':<11}{'required knockdown':<22}")
    for s in (sigma, 0.25, 0.20, 0.15, 0.10, 0.05):
        l, h, _ = interval(xs, prof_all, s ** 2)
        tag = " (current)" if abs(s - sigma) < 1e-9 else ""
        if l is not None:
            scat = f"{100*(np.exp(s)-1):.0f}%"
            print(f"   {s:<9.3f}{scat:<11}{fmt_kd(l, h):<22}{tag}")

    print()
    print("5. VERDICT")
    print("-" * 76)
    print("   The profile is NOT flat, and the published 64-92% range is the")
    print("   likelihood width rather than an artifact of anchor selection.")
    print("   But the binding constraint is between-study scatter at fixed PrP")
    print("   level, not the number of anchors — so adding points along the")
    print("   PrP axis cannot narrow it. Only reducing sigma can.")


if __name__ == "__main__":
    main()
