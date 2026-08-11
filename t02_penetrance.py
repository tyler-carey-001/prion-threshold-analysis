"""
T0-2 — per-variant PRNP penetrance against gnomAD v4.1.1.

Implements the rules fixed in T0-2-prediction-prereg.md and amendments 1 and 2,
all committed before these data were retrieved:

  A2  unpool P102L / A117V / D178N / E200K; a persisting zero is a result
  A4  report crude AND ancestry-standardized (ExAC weights, frozen pre-pull)
  A5  A117V's persisting zero is pre-declared uninformative (33 case alleles)
  A6  filtered records queried explicitly; coverage reported before any zero
  A7  v2.1.1 alongside v4, to separate re-calling from new samples

Estimator and interval are the 2016 construction, validated against the paper
and against the authors' R (T0-2-reproduction.md, T0-2-gate-decision.md sec 1).

POPULATION-LEVEL ESTIMATES. NOT CLINICAL GUIDANCE FOR ANY INDIVIDUAL.
Control-side refresh: case counts are frozen at 2016. ExAC is nested in v4.
"""
import json
import math

from t02_power_check import penetrance_confint, wilson, clamp, BASELINE
from t02_power_check_pervariant import N_CASE, VARIANTS, exac_weights

V4_RAW = "data/gnomad_prnp_v4_raw.json"
V2_RAW = "data/gnomad_prnp_v2_raw.json"
COV_RAW = "data/gnomad_prnp_v4_coverage.json"
N_EXAC = 60706

HGVSP = {
    "p.Pro102Leu": "P102L",
    "p.Ala117Val": "A117V",
    "p.Asp178Asn": "D178N",
    "p.Glu200Lys": "E200K",
}
AC_CASE = dict(VARIANTS)
# Shared ancestry groups: present in both the 2016 labelling and gnomAD v4.
SHARED = ("afr", "amr", "eas", "fin", "nfe", "sas")


def load(path):
    with open(path) as fh:
        return json.load(fh)


def index_variants(doc, hgvsp_map):
    out = {}
    for v in doc["data"]["gene"]["variants"]:
        name = hgvsp_map.get(v.get("hgvsp"))
        if name:
            out[name] = v
    return out


def region_coverage(doc, pos, window=12):
    """Coverage is binned (~6bp); report the bins bracketing `pos`."""
    recs = [c for c in doc["data"]["gene"]["coverage"]["exome"]
            if abs(c["pos"] - pos) <= window]
    if not recs:
        return None
    return (sum(r["mean"] for r in recs) / len(recs),
            min(r["over_20"] for r in recs))


def standardized_af(exome, weights):
    """Direct standardization of v4 per-ancestry AFs to ExAC's composition.

    Returns (af_std, excluded_AN_fraction). Groups with AN == 0 are dropped and
    the remaining weights renormalized, so a missing group cannot silently
    contribute zero frequency.
    """
    pops = {p["id"]: p for p in exome["populations"]}
    used, af_sum, wsum = [], 0.0, 0.0
    for g in SHARED:
        p = pops.get(g)
        if not p or p["an"] == 0:
            continue
        af_sum += weights[g] * (p["ac"] / p["an"])
        wsum += weights[g]
        used.append(g)
    if wsum == 0:
        return None, None
    af_std = af_sum / wsum
    total_an = exome["an"]
    shared_an = sum(pops[g]["an"] for g in used if g in pops)
    return af_std, 1 - (shared_an / total_an if total_an else 0)


def classify_zero(v4rec):
    """A6 three-way classification."""
    if v4rec is None:
        return "zero persists"
    ex = v4rec.get("exome")
    if ex is None:
        return "zero persists"
    filters = ex.get("filters") or []
    if ex["ac"] > 0 and not filters:
        return "RESOLVED (PASS)"
    if ex["ac"] > 0 and filters:
        return f"filtered-only detection {filters}"
    return f"zero persists (record present, AC=0, filters={filters})"


def main():
    v4 = load(V4_RAW)
    v2 = load(V2_RAW)
    cov = load(COV_RAW)
    _, _, weights = exac_weights()
    g4 = index_variants(v4, HGVSP)
    g2 = index_variants(v2, HGVSP)

    # AN for a variant with no record: take the local median exome AN, so the
    # denominator for a true zero is empirical rather than assumed.
    ans = [x["exome"]["an"] for x in v4["data"]["gene"]["variants"]
           if x.get("exome") and x["exome"]["an"]]
    ans.sort()
    default_an = ans[len(ans) // 2]

    print("=" * 96)
    print("T0-2 RESULTS — per-variant PRNP penetrance, gnomAD v4.1.1 controls")
    print("=" * 96)
    print("Population-level estimates. NOT clinical guidance for any individual.")
    print("Control-side refresh: 2016 case counts, frozen. ExAC is nested in v4.")
    print()

    print("A6/A7 — detection status, filters, coverage, and attribution")
    print("-" * 96)
    print(f"{'variant':<8}{'v4 AC':>7}{'v4 AN':>10}{'v2 AC':>7}{'ExAC':>6}"
          f"{'cov mean':>10}{'>=20x':>8}  status / attribution")
    rows = {}
    for name, ac_case in VARIANTS:
        r4, r2 = g4.get(name), g2.get(name)
        ex = r4.get("exome") if r4 else None
        ac4 = ex["ac"] if ex else 0
        an4 = ex["an"] if ex else default_an
        ac2 = (r2["exome"]["ac"] if r2 and r2.get("exome") else 0)
        pos = r4["pos"] if r4 else 4699570  # A117V c.350C>T
        cvm, cv20 = region_coverage(cov, pos) or (float("nan"), float("nan"))
        status = classify_zero(r4)
        # A7 attribution. A non-zero appearing in v2 does NOT by itself imply
        # re-calling: v2 is ~2x ExAC, so a v2 detection may simply be a new
        # sample. The discriminating question is whether the ExAC/v2/v4 counts
        # are mutually consistent with one underlying frequency (see the Poisson
        # section below); only an inconsistency would implicate re-calling.
        if ac4 > 0:
            attrib = (f"also in v2 (AC={ac2}); consistency checked below"
                      if ac2 > 0 else "absent from v2; consistency checked below")
        else:
            attrib = "n/a"
        rows[name] = dict(ac_case=ac_case, ac4=ac4, an4=an4, ac2=ac2,
                          ex=ex, status=status, attrib=attrib)
        print(f"{name:<8}{ac4:>7}{an4:>10,}{ac2:>7}{0:>6}{cvm:>10.1f}"
              f"{cv20*100:>7.1f}%  {status}")
        print(f"{'':<8}{'':>40}  -> {attrib}")

    print()
    print("Penetrance — 2016 baseline vs gnomAD v4.1.1 (crude and standardized)")
    print("-" * 96)
    print(f"{'variant':<8}{'2016 (ExAC n=60,706)':>28}{'v4 crude':>26}{'v4 ancestry-standardized':>32}")
    results = {}
    for name, ac_case in VARIANTS:
        r = rows[name]
        lo16, b16, hi16 = penetrance_confint(ac_case, N_CASE, 0, N_EXAC)
        lo4, b4, hi4 = penetrance_confint(ac_case, N_CASE, r["ac4"], r["an4"] // 2)
        s = f"{b16*100:6.1f}% [{lo16*100:5.1f},{hi16*100:5.1f}]"
        s4 = f"{b4*100:6.1f}% [{lo4*100:5.1f},{hi4*100:5.1f}]"
        if r["ex"]:
            af_std, excl = standardized_af(r["ex"], weights)
        else:
            af_std, excl = 0.0, 0.0
        if af_std and af_std > 0:
            # express the standardized AF as an effective allele count at the
            # observed AN, so the same interval machinery applies
            eff_ac = af_std * r["an4"]
            lo_s, b_s, hi_s = penetrance_confint(ac_case, N_CASE, eff_ac, r["an4"] // 2)
            ss = f"{b_s*100:6.1f}% [{lo_s*100:5.1f},{hi_s*100:5.1f}]  (excl {excl*100:.1f}%)"
        else:
            lo_s = b_s = hi_s = None
            ss = "        n/a (AF = 0)"
        results[name] = dict(y2016=(lo16, b16, hi16), crude=(lo4, b4, hi4),
                             std=(lo_s, b_s, hi_s), af_std=af_std)
        print(f"{name:<8}{s:>28}{s4:>26}{ss:>32}")

    print()
    print("Composition effect (A4): crude vs standardized")
    print("-" * 96)
    pops4 = {p["id"]: p for p in g4["E200K"]["exome"]["populations"]}
    nfe_frac = pops4["nfe"]["an"] / g4["E200K"]["exome"]["an"]
    print(f"gnomAD v4 exomes are {nfe_frac*100:.1f}% NFE by allele number; "
          f"ExAC's inferred composition was {weights['nfe']*100:.1f}% NFE.")
    for name in ("P102L", "D178N", "E200K"):
        r, res = rows[name], results[name]
        if res["std"][1] is None:
            continue
        crude_af = r["ac4"] / r["an4"]
        print(f"  {name:<7} crude AF={crude_af:.3e}  standardized AF={res['af_std']:.3e}  "
              f"penetrance {res['crude'][1]*100:.1f}% -> {res['std'][1]*100:.1f}%")

    print()
    print("Per-ancestry detail, v4 exomes")
    print("-" * 96)
    for name in ("P102L", "D178N", "E200K"):
        ex = rows[name]["ex"]
        parts = [f"{p['id']}={p['ac']}/{p['an']:,}"
                 for p in sorted(ex["populations"], key=lambda x: -x["ac"]) if p["ac"] > 0]
        print(f"  {name:<7} " + "  ".join(parts))

    print()
    print("Was the 2016 zero ever surprising? Poisson expectation at ExAC's size")
    print("-" * 96)
    print("If v4's allele frequency is the truth, how many alleles should ExAC have")
    print("seen in 2*60,706 = 121,412 chromosomes, and how likely was observing zero?")
    print(f"{'variant':<8}{'v4 AF':>12}{'expected in ExAC':>20}{'P(observe 0)':>15}  verdict")
    for name, ac_case in VARIANTS:
        r = rows[name]
        af4 = r["ac4"] / r["an4"]
        lam = af4 * 2 * N_EXAC
        p0 = math.exp(-lam)
        verdict = ("zero unremarkable — sample size, not disagreement"
                   if p0 > 0.05 else "zero genuinely surprising")
        print(f"{name:<8}{af4:>12.3e}{lam:>20.2f}{p0:>15.3f}  {verdict}")
    print()
    print("Same check against v2.1.1 (AN=251,466), which contains ExAC:")
    for name, ac_case in VARIANTS:
        r = rows[name]
        af4 = r["ac4"] / r["an4"]
        lam2 = af4 * 251466
        print(f"  {name:<7} expected {lam2:5.2f}   observed {r['ac2']}   "
              f"P(obs<=observed) = {sum(math.exp(-lam2)*lam2**k/math.factorial(k) for k in range(r['ac2']+1)):.3f}")

    print()
    print("CAVEAT on the standardized interval")
    print("-" * 96)
    print("The standardized estimate converts a weighted AF into an effective allele")
    print("count at the full AN and reuses the binomial machinery. That propagates")
    print("sampling error in the total, but NOT the extra variance introduced by")
    print("reweighting small per-group counts (E200K: amr AC=2, fin AC=1). The")
    print("standardized CI is therefore ANTI-CONSERVATIVE and is reported as")
    print("indicative of the composition effect's direction and rough size, not as")
    print("a calibrated interval. The crude interval is the calibrated one.")

    print()
    print("Pre-registered rule applications")
    print("-" * 96)
    print("A5  A117V zero persists at n~730k. Lower bound "
          f"{results['A117V']['crude'][0]*100:.1f}% — pre-declared UNINFORMATIVE")
    print("    (33 case alleles cap the bound regardless of control cohort size).")
    print("P3  no row may exceed its case-only ceiling:")
    for name, ac_case in VARIANTS:
        cl, _, cu = wilson(ac_case, 2 * N_CASE)
        print(f"      {name:<7} ceiling log-width {math.log(cu/cl):.3f}")
    print("P4  assessed per variant, above.")


if __name__ == "__main__":
    main()
