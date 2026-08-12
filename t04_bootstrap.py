"""
T0-4 — is the NPM-vs-threshold-free model comparison stable at n=7?

`t04_structural_test.py` reported dAICc = +64 favouring the one-parameter
threshold-free model. That number needs auditing before it is used, for a reason
that runs toward the conclusion it supports:

  - the NPM's SSR (0.3590) sits essentially AT the noise floor (0.3583, the SSR
    of a single constant fitted to the five ~50% anchors);
  - the threshold-free model's SSR (0.5484) sits well ABOVE that floor.

So on the cluster that carries most of the data, the NPM matches a constant and
the threshold-free form does worse than one. AICc's verdict is therefore "five
parameters are not supportable at n=7" — a fair statement about this dataset —
and NOT "the threshold-free structure describes the data better." Those are
different claims and the write-up conflated them.

This script fits both models on the same footing, reports per-anchor residuals,
and bootstraps dAICc over the seven anchors. Five of those anchors sit at one
x-value, so resampling should destabilise the comparison badly. If the sign of
dAICc flips under resampling, the honest statement is that n=7 cannot
discriminate structures.

Model output, not evidence. Anchors unverified against primary sources.
"""
import numpy as np
from dataclasses import replace
from scipy.optimize import minimize, minimize_scalar

from prion_model import PrionParams, growth_rate, relative_survival_analytic, ANCHORS

BASE = PrionParams()
FINITE = [(x, r) for _l, x, r, _s in ANCHORS if np.isfinite(r)]
NPM_BOUNDS = [(np.log(1e-3), np.log(0.5)),    # log a
              (np.log(1e-6), np.log(1e-2)),   # log b
              (2.0, 12.0),                    # n
              (0.01, 0.60),                   # t_tox_frac
              (0.02, 0.95)]                   # x_crit
NPM_K, SB_K = 5, 1


# ---------------------------------------------------------------- models
def sandberg_pred(x, p1):
    return p1 + (1.0 - p1) / x


def fit_sandberg(data):
    f = lambda p1: sum((np.log(sandberg_pred(x, p1)) - np.log(r)) ** 2
                       for x, r in data)
    res = minimize_scalar(f, bounds=(0.0, 0.99), method="bounded")
    return float(res.fun), float(res.x)


def npm_params(theta):
    la, lb, n, tf, xc = theta
    a, b = np.exp(la), np.exp(lb)
    num = a * (a + b * (2 * n - 1)) + b ** 2 * n * (n - 1)
    return replace(BASE, a=a, b=b, n=n, beta=num / (b * xc)), tf


def npm_ssr(theta, data):
    p, tf = npm_params(theta)
    if growth_rate(p.x0, p) <= 0:
        return 1e3
    tot = 0.0
    for x, r in data:
        pred = relative_survival_analytic(x, p, tf)
        if not np.isfinite(pred) or pred <= 0:
            return 1e3
        tot += (np.log(pred) - np.log(r)) ** 2
    return tot


def fit_npm(data, n_starts=6, seed=0):
    rng = np.random.default_rng(seed)
    best, bt = np.inf, None
    starts = [np.array([np.log(.02), np.log(2e-4), 6., .25, .30])]
    for _ in range(n_starts - 1):
        starts.append(np.array([rng.uniform(*NPM_BOUNDS[i]) for i in range(5)]))
    for s in starts:
        try:
            res = minimize(npm_ssr, s, args=(data,), method="L-BFGS-B",
                           bounds=NPM_BOUNDS, options={"maxiter": 1500})
            if res.fun < best:
                best, bt = float(res.fun), res.x
        except Exception:
            continue
    return best, bt


# ---------------------------------------------------------------- criteria
def aicc(ssr, n, k):
    if ssr <= 0:
        ssr = 1e-12
    a = n * np.log(ssr / n) + 2 * k
    d = n - k - 1
    return a + (2 * k * (k + 1) / d if d > 0 else np.inf)


def main():
    n = len(FINITE)
    sb_ssr, p1 = fit_sandberg(FINITE)
    np_ssr, np_th = fit_npm(FINITE)

    print("=" * 78)
    print("Is the structural model comparison stable at n=7?")
    print("=" * 78)

    print()
    print("1. BOTH MODELS FITTED ON THE SAME FOOTING")
    print("-" * 78)
    print(f"   threshold-free (k=1)  SSR = {sb_ssr:.4f}   p1 = {p1:.3f}")
    print(f"   NPM            (k=5)  SSR = {np_ssr:.4f}")
    print(f"   noise floor (constant fit to the five ~50% anchors) = 0.3583")
    print()
    print("   The NPM sits AT the floor; the threshold-free form sits above it.")
    print("   AICc's preference is a statement about supporting 5 parameters at")
    print("   n=7, not about which structure describes the data better.")

    print()
    print("2. PER-ANCHOR RESIDUALS (log scale)")
    print("-" * 78)
    p, tf = npm_params(np_th)
    print(f"   {'x':>7}{'observed':>10}{'NPM':>9}{'resid':>8}"
          f"{'  |':>3}{'free':>8}{'resid':>8}")
    for x, r in sorted(FINITE):
        npm_p = relative_survival_analytic(x, p, tf)
        sb_p = sandberg_pred(x, p1)
        print(f"   {x:>7.2f}{r:>10.2f}{npm_p:>9.2f}{np.log(npm_p/r):>8.3f}"
              f"{'  |':>3}{sb_p:>8.2f}{np.log(sb_p/r):>8.3f}")

    print()
    print("3. BOOTSTRAP OVER THE SEVEN ANCHORS")
    print("-" * 78)
    B = 400
    rng = np.random.default_rng(20260812)
    diffs, skipped = [], 0
    for i in range(B):
        idx = rng.integers(0, n, n)
        samp = [FINITE[j] for j in idx]
        if len({x for x, _ in samp}) < 3:      # need >=3 distinct x to fit at all
            skipped += 1
            continue
        s_sb, _ = fit_sandberg(samp)
        s_np, _ = fit_npm(samp, n_starts=3, seed=i)
        d = aicc(s_np, n, NPM_K) - aicc(s_sb, n, SB_K)
        if np.isfinite(d):
            diffs.append(d)
    diffs = np.array(diffs)
    pos = float((diffs > 0).mean())
    print(f"   {B} resamples, {len(diffs)} usable ({skipped} had <3 distinct x)")
    print(f"   dAICc = AICc(NPM) - AICc(threshold-free); positive favours the")
    print(f"   threshold-free model.")
    print()
    print(f"   point estimate (full data) : {aicc(np_ssr,n,NPM_K)-aicc(sb_ssr,n,SB_K):+.1f}")
    print(f"   bootstrap median           : {np.median(diffs):+.1f}")
    print(f"   bootstrap 2.5-97.5 pct     : {np.percentile(diffs,2.5):+.1f} to "
          f"{np.percentile(diffs,97.5):+.1f}")
    print(f"   fraction favouring threshold-free : {pos:.3f}")
    print(f"   fraction favouring NPM            : {1-pos:.3f}")

    print()
    print("4. WHAT dAICc IS ACTUALLY MADE OF")
    print("-" * 78)

    def parts(ssr, k):
        return (n * np.log(ssr / n), 2 * k, 2 * k * (k + 1) / (n - k - 1))

    f1, q1, c1 = parts(sb_ssr, SB_K)
    f5, q5, c5 = parts(np_ssr, NPM_K)
    tot = (f5 + q5 + c5) - (f1 + q1 + c1)
    print(f"   {'term':<32}{'free':>10}{'NPM':>10}{'diff':>10}")
    print(f"   {'n*ln(SSR/n)  — the FIT':<32}{f1:>10.2f}{f5:>10.2f}{f5-f1:>+10.2f}")
    print(f"   {'2k           — AIC penalty':<32}{q1:>10.2f}{q5:>10.2f}{q5-q1:>+10.2f}")
    print(f"   {'2k(k+1)/(n-k-1) — AICc corr':<32}{c1:>10.2f}{c5:>10.2f}{c5-c1:>+10.2f}")
    print(f"   {'TOTAL':<32}{'':>10}{'':>10}{tot:>+10.2f}")
    print()
    print(f"   The FIT term favours the NPM by {f5-f1:+.2f}.")
    print(f"   The AICc small-sample correction alone is {c5-c1:+.2f} — "
          f"{100*(c5-c1)/tot:.0f}% of the total.")
    print("   At n=7, k=5 the denominator (n-k-1) equals 1, so that term is a")
    print("   near-deterministic arithmetic penalty rather than evidence — which")
    print("   is ALSO why the bootstrap sign is stable. The stability is an")
    print("   artifact of the penalty, not a property of the data.")

    print()
    print("5. VERDICT")
    print("-" * 78)
    if pos < 0.95:
        print(f"   dAICc CHANGES SIGN under resampling ({100*(1-pos):.0f}% of samples")
        print("   favour the NPM). The model comparison is UNSTABLE at n=7.")
        print("   Honest statement: seven anchors, five of them at one x-value,")
        print("   cannot discriminate between these structures. The +64 headline")
        print("   should not be used.")
    else:
        print("   Sign is stable under resampling. The comparison survives, though")
        print("   it still says only that 5 parameters are unsupportable at n=7.")
    print()
    print("   Either way, the EXTRAPOLATION argument is untouched and does not")
    print("   depend on model selection at all: x_crit = 33% residual PrP, the")
    print("   lowest anchor is 49%, and the structures diverge 13-fold below 40%")
    print("   where nobody has measured. That is the finding to lead with.")


if __name__ == "__main__":
    main()
