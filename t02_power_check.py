"""
T0-2 Phase 0 — power check for the PRNP penetrance refresh.

Answers one question before any gnomAD data is pulled: can a larger control
cohort actually tighten the 2016 confidence intervals, and by how much?

The intervals are recomputed here under the *repo's own* interval construction
(Kirov et al. 2014 opposite-corner, prnp_penetrance/src/generate_figures.r:120-129),
not under a delta-method approximation. Those two disagree: the corner method is
additive in the component log-widths, the delta method is root-sum-square.

This is NOT the reproduction gate. This is an independent Python recomputation
used to size the task. The reproduction gate (Phase A) runs the original R and
checks it against numbers printed in the paper.

Inputs are the 2016 supplementary tables, read from a local clone of
github.com/ericminikel/prnp_penetrance @ da681e0.
"""
import math

Z = 1.959963984540054  # qnorm(0.975), matching R binom.confint
BASELINE = 2e-4        # assumed_baseline_risk, generate_figures.r:114


def wilson(x, n, z=Z):
    """Wilson score interval, matching R binom.confint(method='wilson')."""
    if n == 0:
        return (0.0, 0.0, 0.0)
    center = (x + z * z / 2) / (n + z * z)
    half = (z / (n + z * z)) * math.sqrt(x * (n - x) / n + z * z / 4)
    return (max(0.0, center - half), x / n, min(1.0, center + half))


def clamp(v):
    return min(1.0, max(0.0, v))


def penetrance_confint(ac_case, n_case, ac_control, n_control):
    """Opposite-corner interval, generate_figures.r:120-129.

    lower = case_lower / control_upper ; upper = case_upper / control_lower.
    Both cohorts propagate. The allelic model doubles both denominators.
    """
    cl, cm, cu = wilson(ac_case, 2 * n_case)
    kl, km, ku = wilson(ac_control, 2 * n_control)
    lower = clamp(cl * BASELINE / ku) if ku > 0 else 0.0
    best = clamp(cm * BASELINE / km) if km > 0 else 1.0
    upper = clamp(cu * BASELINE / kl) if kl > 0 else 1.0
    return lower, best, upper


def logwidth(lo, hi):
    """Width on the log scale; inf if either end is clamp-censored at 0 or 1."""
    if lo <= 0 or hi <= 0:
        return math.inf
    return math.log(hi / lo)


def ceiling_tightening(ac_case, n_case, ac_control, n_control):
    """Max achievable tightening: control uncertainty -> 0, leaving a case-only
    interval. Independent of how large the new control cohort actually is."""
    lo, _, hi = penetrance_confint(ac_case, n_case, ac_control, n_control)
    w0 = logwidth(lo, hi)
    if ac_control == 0 or not math.isfinite(w0):
        return None
    cl, _, cu = wilson(ac_case, 2 * n_case)
    km = ac_control / (2 * n_control)
    return 1 - logwidth(clamp(cl * BASELINE / km), clamp(cu * BASELINE / km)) / w0


def min_multiplier(ac_case, n_case, ac_control, n_control, target):
    """Smallest control-cohort multiplier reaching `target` fractional tightening."""
    ceil = ceiling_tightening(ac_case, n_case, ac_control, n_control)
    if ceil is None or target > ceil:
        return None
    lo, _, hi = penetrance_confint(ac_case, n_case, ac_control, n_control)
    w0 = logwidth(lo, hi)
    m = 1.0
    while m < 1e7:
        l2, _, h2 = penetrance_confint(ac_case, n_case, ac_control * m, n_control * m)
        if 1 - logwidth(l2, h2) / w0 >= target:
            return m
        m *= 1.02
    return None


# The frozen variant list: every forest-plot row that uses ExAC as its control
# cohort (generate_figures.r:135-161). 23andMe rows are excluded because a gnomAD
# refresh does not touch them. This list is fixed by the 2016 paper's own
# selection and MUST NOT be extended with variants that look interesting in v4.
#
# Note the control cohorts: the 2016 analysis is already ancestry-matched, using
# single ExAC subpopulations (JPT n=663, TSI n=4795), not the full 60,706.
ROWS = [
    # label,              ac_case, n_case, ac_control, n_control, control cohort
    ("M232R",                  67,   1533,          5,       663, "ExAC JPT"),
    ("V180I",                 222,   1533,          2,       663, "ExAC JPT"),
    ("V210I",                 171,   1054,          2,      4795, "ExAC TSI"),
    ("Mendelians(pooled)",   1034,  10460,          0,     60706, "ExAC all"),
]

RULE = 78  # display width


def main():
    print("2016 ExAC forest rows, recomputed under the repo's corner method")
    print("-" * RULE)
    print(f"{'row':<20}{'ac_case':>8}{'ac_ctl':>7}{'n_ctl':>7}"
          f"{'penetrance (95% CI)':>32}")
    for label, ac, nc, kc, kn, _ in ROWS:
        lo, best, hi = penetrance_confint(ac, nc, kc, kn)
        ci = f"{best*100:6.2f}%  [{lo*100:5.2f}%, {hi*100:6.2f}%]"
        print(f"{label:<20}{ac:>8}{kc:>7}{kn:>7}{ci:>32}")

    print()
    print("Variance decomposition and tightening ceiling")
    print("-" * RULE)
    print(f"{'row':<20}{'1/x_case':>10}{'1/x_ctl':>9}{'dominated by':>14}{'ceiling':>10}")
    for label, ac, nc, kc, kn, _ in ROWS:
        inv_case = 1 / ac if ac else math.inf
        inv_ctl = 1 / kc if kc else math.inf
        ceil = ceiling_tightening(ac, nc, kc, kn)
        dom = "control" if inv_ctl > inv_case else "case"
        ceil_s = f"{ceil*100:>9.1f}%" if ceil is not None else f"{'n/a':>10}"
        ictl_s = f"{inv_ctl:>9.4f}" if math.isfinite(inv_ctl) else f"{'inf':>9}"
        print(f"{label:<20}{inv_case:>10.4f}{ictl_s}{dom:>14}{ceil_s}")

    print()
    print("Minimum control-cohort multiplier to reach a tightening target")
    print("-" * RULE)
    print(f"{'row':<20}{'25%':>10}{'50%':>10}{'75%':>10}")
    for label, ac, nc, kc, kn, _ in ROWS:
        cells = []
        for t in (0.25, 0.50, 0.75):
            m = min_multiplier(ac, nc, kc, kn, t)
            cells.append(f"{m:.1f}x" if m else "n/a")
        print(f"{label:<20}" + "".join(f"{c:>10}" for c in cells))

    print()
    print("Zero-count stratum: Mendelians(pooled), ac_control = 0")
    print("-" * RULE)
    print("Upper end is clamp-censored at 100%, so log-width is undefined and the")
    print("row cannot 'tighten'. More zero-allele controls push the LOWER bound UP:")
    for n_ctl in (60706, 150000, 303530, 730947):
        lo, best, hi = penetrance_confint(1034, 10460, 0, n_ctl)
        print(f"   n_control={n_ctl:>7}, ac=0  ->  [{lo*100:6.2f}%, {hi*100:6.1f}%]")
    print()
    print("So this row moves only when the zero resolves, i.e. when gnomAD shows")
    print("ac_control > 0. That is a different question from interval tightening.")

    print()
    print("Paper cross-check (generate_figures.r:288-289)")
    print("-" * RULE)
    lo, _, hi = penetrance_confint(1034, 10460, 0, 60706)
    ok = lo <= 0.60 and hi >= 0.90
    print(f"Paper states the Mendelians_ExAC CI contains published E200K survival")
    print(f"estimates of ~60-90%.  Computed: [{lo*100:.1f}%, {hi*100:.1f}%]  contains 60-90%: {ok}")


if __name__ == "__main__":
    main()
