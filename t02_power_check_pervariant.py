"""
T0-2 — per-variant power check for the unpooled Mendelian variants.

Run BEFORE pulling gnomAD. Unpooling changes the variance structure: the pooled
row had 1,034 case alleles, but the four variants individually span 33 to 571,
a 17-fold spread in 1/x_case. A variant that is case-dominated post-unpooling
has an interval that barely moves no matter what gnomAD shows, and that needs to
be known in advance rather than discovered afterwards.

Also computes the ExAC ancestry weights used for direct standardization, from
table_s08 -- pre-pull, so the weights cannot be chosen after seeing v4.
"""
import math
from t02_power_check import wilson, clamp, penetrance_confint, logwidth, BASELINE

N_CASE = 10460          # table_s01 TOTAL, prnp_sequenced
N_V4_EXOMES = 730947    # gnomAD v4.1 exomes; joint/coverage decision deferred to the pull

VARIANTS = [            # per-variant case counts, table_s01 TOTAL row
    ("P102L", 221),
    ("A117V", 33),
    ("D178N", 209),
    ("E200K", 571),
]

# ExAC ancestry composition as inferred by the 2016 paper (table_s08), collapsed
# to the gnomAD continental groups. FIN is separate from NFE in gnomAD.
EXAC_ANCESTRY = {
    "afr": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "amr": ["CLM", "MXL", "PEL", "PUR"],
    "eas": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "fin": ["FIN"],
    "nfe": ["CEU", "GBR", "IBS", "TSI"],
    "sas": ["BEB", "GIH", "ITU", "PJL", "STU"],
}
EXAC_POP_N = {
    "ACB": 2267, "ASW": 2151, "BEB": 483, "CDX": 19, "CEU": 14185, "CHB": 1553,
    "CHS": 1733, "CLM": 870, "ESN": 89, "FIN": 3977, "GBR": 10358, "GIH": 79,
    "GWD": 102, "IBS": 3534, "ITU": 1089, "JPT": 663, "KHV": 369, "LWK": 72,
    "MSL": 189, "MXL": 2658, "PEL": 1900, "PJL": 6300, "PUR": 579, "STU": 460,
    "TSI": 4795, "YRI": 232,
}


def exac_weights():
    """Direct-standardization weights: ExAC ancestry proportions over the groups
    shared with gnomAD v4. gnomAD's asj / mid / remaining have no ExAC counterpart
    in this labelling and are excluded; the excluded v4 fraction is reported at
    analysis time."""
    counts = {g: sum(EXAC_POP_N[p] for p in pops) for g, pops in EXAC_ANCESTRY.items()}
    total = sum(counts.values())
    return counts, total, {g: c / total for g, c in counts.items()}


def main():
    print("Per-variant case-side variance, post-unpooling")
    print("-" * 74)
    print(f"{'variant':<9}{'ac_case':>9}{'af_case':>10}{'1/x_case':>11}"
          f"{'x_ctrl for parity':>20}")
    for name, ac in VARIANTS:
        inv = 1 / ac
        print(f"{name:<9}{ac:>9}{ac/(2*N_CASE)*100:>9.3f}%{inv:>11.4f}{ac:>20}")
    print()
    print("'x_ctrl for parity' = control alleles needed before the control term")
    print("stops dominating. Below that count the interval is control-limited and")
    print("gnomAD helps; above it the frozen case count is the binding constraint.")

    print()
    print("Interval vs hypothetical v4 control count (n_control = %s)" % f"{N_V4_EXOMES:,}")
    print("-" * 74)
    hyp = [0, 1, 2, 5, 10, 20, 50, 100]
    print(f"{'variant':<9}" + "".join(f"{('ac='+str(h)):>8}" for h in hyp))
    for name, ac in VARIANTS:
        cells = []
        for h in hyp:
            lo, best, hi = penetrance_confint(ac, N_CASE, h, N_V4_EXOMES)
            cells.append(f"{best*100:>7.1f}%")
        print(f"{name:<9}" + "".join(cells))
    print("(point estimates; ac=0 clamps to 100% by construction)")

    print()
    print("If the zero holds at n=%s: lower bound on penetrance" % f"{N_V4_EXOMES:,}")
    print("-" * 74)
    for name, ac in VARIANTS:
        lo, best, hi = penetrance_confint(ac, N_CASE, 0, N_V4_EXOMES)
        lo16, _, _ = penetrance_confint(ac, N_CASE, 0, 60706)
        print(f"{name:<9} 2016-size lower {lo16*100:>6.1f}%   ->  v4-size lower {lo*100:>6.1f}%")
    print("A persisting zero is a tight lower bound on high penetrance, i.e. a")
    print("result. It is not a null and is not 'no data'.")

    print()
    print("Informativeness ceiling: interval width with control uncertainty -> 0")
    print("-" * 74)
    print(f"{'variant':<9}{'case-only log-width':>22}{'as fold-range':>16}")
    for name, ac in VARIANTS:
        cl, cm, cu = wilson(ac, 2 * N_CASE)
        w = math.log(cu / cl)
        print(f"{name:<9}{w:>22.3f}{math.exp(w):>15.1f}x")
    print("This is the narrowest interval obtainable for each variant at ANY")
    print("control cohort size. A117V's floor is set by 33 case alleles.")

    print()
    print("ExAC ancestry weights for direct standardization (pre-registered pre-pull)")
    print("-" * 74)
    counts, total, w = exac_weights()
    print(f"{'group':<8}{'ExAC n':>10}{'weight':>10}")
    for g in sorted(w, key=lambda k: -w[k]):
        print(f"{g:<8}{counts[g]:>10,}{w[g]*100:>9.2f}%")
    print(f"{'TOTAL':<8}{total:>10,}{100.0:>9.2f}%")
    assert total == 60706, total
    print("Sums to the published ExAC total of 60,706. gnomAD asj / mid /")
    print("remaining have no counterpart in this labelling and are excluded;")
    print("the excluded v4 fraction is reported when the standardization is run.")


if __name__ == "__main__":
    main()
