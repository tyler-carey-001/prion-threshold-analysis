"""
T0-4 — structural test: does the dose-response support a THRESHOLD at all?

Found during the T0-4 literature check, which was itself added after T0-2's
process failure. Sandberg et al. 2011 (Nature 470:540) — already cited in this
repo's README — measured the two-phase kinetics directly, and its structure is
the reverse of this model's.

Verified quotes (Sandberg 2014, PMC4104459, restating the 2011 result):

  "In phase 1, prions propagate exponentially, not rate limited by cellular
   prion protein (PrP^C) concentration, rapidly rising by ~10^6-fold to reach a
   maximal prion titre, which is also independent of PrP^C concentration over
   the range we studied."

  "This is followed by a plateau phase (phase 2), which determines time to
   clinical onset of disease, the duration of which is inversely proportional to
   PrP^C concentration."

  "Remarkably, total PrP levels, although starting at different levels in the
   three mouse lines due to their different PrP^C expression levels, were
   essentially unchanged throughout phase 1 in wild-type, Prnp+/o and Tg20 mice,
   a period during which titres rose in each case by ~10^6-fold."

This model (`prion_model.relative_survival_analytic`) assumes the OPPOSITE
assignment: the replication phase carries the PrP dependence via r(x), and the
toxic phase is a fixed, PrP-independent fraction. Sandberg puts the PrP
dependence in the plateau and finds the exponential phase PrP-independent.

Taking Sandberg's mechanism literally gives a one-parameter model with NO
threshold:

    survival_ratio(x) = p1 + (1 - p1) / x

where p1 is the fraction of the untreated course spent in phase 1. This script
fits it against the same anchors and compares.

Model output, not evidence. Anchors unverified against primary sources.
"""
import numpy as np
from dataclasses import replace
from scipy.optimize import minimize_scalar

from prion_model import PrionParams, growth_rate, x_crit, ANCHORS

FINITE = [(x, r) for _l, x, r, _s in ANCHORS if np.isfinite(r)]
NPM_SSR = 0.3590      # frozen in T0-4-prereg.md
NPM_K = 5             # x_crit + a, b, n, t_tox_frac


def sandberg_ratio(x, p1):
    return p1 + (1.0 - p1) / x


def ssr_sandberg(p1, data=FINITE):
    return sum((np.log(sandberg_ratio(x, p1)) - np.log(r)) ** 2 for x, r in data)


def aic(ssr, n, k):
    return n * np.log(ssr / n) + 2 * k


def aicc(ssr, n, k):
    a = aic(ssr, n, k)
    denom = n - k - 1
    return a + (2 * k * (k + 1) / denom if denom > 0 else np.inf)


def main():
    n = len(FINITE)
    res = minimize_scalar(ssr_sandberg, bounds=(0.0, 0.99), method="bounded")
    p1, s_sb = float(res.x), float(res.fun)

    print("=" * 76)
    print("Does the published dose-response support a PrP threshold at all?")
    print("=" * 76)
    print()
    print("1. TWO MODELS OF THE SAME DATA")
    print("-" * 76)
    print(f"   NPM (this repo)  survival = (1-f)*r(1)/r(x) + f     "
          f"k={NPM_K}  SSR={NPM_SSR:.4f}")
    print(f"     PrP-dependence in the REPLICATION phase; has a threshold x_crit")
    print(f"   Sandberg         survival = p1 + (1-p1)/x           "
          f"k=1  SSR={s_sb:.4f}")
    print(f"     PrP-dependence in the PLATEAU phase; NO threshold")
    print(f"     best-fit p1 = {p1:.3f}  (phase 1 = {100*p1:.0f}% of untreated course)")

    print()
    print("2. MODEL SELECTION")
    print("-" * 76)
    print(f"   {'model':<14}{'k':>3}{'SSR':>10}{'AIC':>10}{'AICc':>12}")
    for lab, ssr_, k in (("Sandberg 1/x", s_sb, 1), ("NPM", NPM_SSR, NPM_K)):
        print(f"   {lab:<14}{k:>3}{ssr_:>10.4f}{aic(ssr_, n, k):>10.2f}"
              f"{aicc(ssr_, n, k):>12.2f}")
    print()
    print(f"   dAIC  = {aic(NPM_SSR,n,NPM_K)-aic(s_sb,n,1):+.2f} favouring Sandberg")
    print(f"   dAICc = {aicc(NPM_SSR,n,NPM_K)-aicc(s_sb,n,1):+.2f} favouring Sandberg")
    print("   AICc is the correct criterion at n=7; the NPM's 5 parameters are")
    print("   barely estimable from 7 points and it is penalised accordingly.")

    print()
    print("3. WHERE THE MODELS DIVERGE — and where the data is")
    print("-" * 76)
    p = PrionParams(beta=6.665562105385422)      # repo's calibrated value
    xc, r1 = x_crit(p), growth_rate(p.x0, p)
    print(f"   calibrated NPM x_crit = {xc:.3f} residual PrP "
          f"({100*(1-xc):.1f}% knockdown)")
    print(f"   lowest anchor in the dataset = "
          f"{min(x for x, _ in FINITE if x > 0):.2f} residual PrP")
    print()
    print(f"   {'residual PrP':>13}{'NPM':>10}{'Sandberg':>11}{'NPM/Sandberg':>14}  data?")
    for x in (8.0, 1.0, 0.75, 0.50, 0.40, 0.35, 0.334, 0.25):
        rx = growth_rate(x * p.x0, p)
        npm = (r1 / rx) if rx > 0 else np.inf
        sb = sandberg_ratio(x, p1)
        has = any(abs(x - ax) < 0.02 for ax, _ in FINITE)
        npm_s = f"{npm:>10.2f}" if np.isfinite(npm) else f"{'inf':>10}"
        rat = f"{npm/sb:>14.1f}" if np.isfinite(npm) else f"{'inf':>14}"
        print(f"   {x:>13.3f}{npm_s}{sb:>11.2f}{rat}  {'YES' if has else '--'}")

    print()
    print("4. THE PROBLEM, STATED PLAINLY")
    print("-" * 76)
    print("   The NPM threshold sits at 33% residual PrP. The lowest anchor is")
    print("   49%. x_crit is therefore an EXTRAPOLATION beyond every data point,")
    print("   not an inference from them.")
    print()
    print("   Sandberg's own hedge is the same one: 'independent of PrP^C")
    print("   concentration OVER THE RANGE WE STUDIED' — and that range was")
    print("   Prnp+/o (50%), wild-type, and Tg20 (8x). Nobody has measured below")
    print("   ~50% residual PrP. Both models are fitted to, and validated on,")
    print("   exactly the region where they agree.")
    print()
    print("   Consequence for the headline: '65-90% knockdown required' is a")
    print("   property of an assumed model structure, not a measured quantity,")
    print("   and a simpler structure taken from the field's own kinetic data")
    print("   fits at least as well and implies NO threshold at all.")
    print()
    print("   RESEARCH_PLAN.md section 6 already named this as the falsification")
    print("   condition: 'If there is no inflection - if survival scales smoothly")
    print("   all the way down - the nucleated-polymerization framing is wrong.'")

    print()
    print("5. WHAT WOULD DISCRIMINATE")
    print("-" * 76)
    print("   The models are within ~1.6x of each other at every existing anchor")
    print("   and diverge 10-fold by 35% residual PrP. A depth titration must")
    print("   therefore go BELOW 40% residual to say anything:")
    for x in (0.40, 0.30, 0.25, 0.15):
        rx = growth_rate(x * p.x0, p)
        npm = (r1 / rx) if rx > 0 else np.inf
        npm_s = "divergent" if not np.isfinite(npm) else f"{npm:.1f}x"
        print(f"     at {100*x:>3.0f}% residual:  NPM predicts {npm_s:>10}   "
              f"Sandberg predicts {sandberg_ratio(x,p1):.1f}x")
    print()
    print("   Those are not subtle differences. One experiment settles it.")


if __name__ == "__main__":
    main()
