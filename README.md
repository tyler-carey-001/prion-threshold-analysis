# Prion disease: a kinetic threshold model, and a PRNP penetrance refresh

Independent, dry-lab computational work on two questions: **how deeply must PrP
be lowered to stop prion replication**, and **how penetrant are the PRNP variants
that cause genetic prion disease**.

---

## If you carry a PRNP variant, please read this part first

**Nothing in this repository is medical advice, and no number here describes your
risk.**

If you have landed here after a genetic test result, the honest summary is short:

1. **These are population averages, not personal risk.** They are computed from
   allele frequencies in large databases, not from anyone's family, age, or
   genotype at other sites.
2. **The E200K figure of 61.4% that appears in these files is not a reliable
   estimate of anything.** It rests on **13 alleles**. It disappears entirely if
   6 of those 13 were miscalled. Three separate sources of systematic error are
   identified in the write-up, they push in **opposite** directions, at least one
   cannot be quantified, and together they are larger than the statistical
   uncertainty. The most careful person working on this — Eric Minikel, who wrote
   the original 2016 analysis — considers **~90%** the more credible figure for
   E200K. This work does not overturn that and does not claim to.
3. **"Lifetime risk" here means risk from birth. It is not your remaining risk.**
   If you are 55 and unaffected, you have already passed through part of the risk
   window, and your remaining risk is *lower* than any from-birth number. This
   analysis does not compute conditional risk, because doing it properly needs
   onset-age data it does not have.
4. **The case counts are frozen at 2016.** Only the control side was updated.
   This is not a full re-analysis.
5. **Codon 129 matters and is not modelled here.** It modifies risk and
   determines whether D178N presents as fatal familial insomnia or as CJD.

**Where to go instead of here:**

- **[CJD Foundation](https://cjdfoundation.org)** — patient and family support,
  and connections to clinicians who actually see this disease.
- **[Prion Alliance](https://www.prionalliance.org)** — founded by two
  researchers, one of whom carries a PRNP variant.
- **[cureffi.org](https://www.cureffi.org)** — Minikel has written openly about
  this field, including penetrance and what it means for carriers, for over a
  decade. It is a better starting point than this repository.
- A **genetic counsellor**, for anything about your own result.

---

## Status and provenance

**This is unreviewed, unpublished, independent work.** It has not been through peer review. Corrections and disagreement are welcome — open an issue.

**Prior work you should read first.** Minikel published the gnomAD v4 PRNP allele
counts in April 2024:
[cureffi.org/2024/04/03/learn-prnp-gnomad-v4](https://www.cureffi.org/2024/04/03/learn-prnp-gnomad-v4/).
His counts and the ones here agree exactly (P102L 2, D178N 1, E200K 13). **The
allele counts in this repository are not a new finding.** That post was found
*after* this analysis was complete, which is a process failure and is recorded as
one. His treatment of age structure and of contaminating cohorts is better than
what is here.

**Written with AI assistance.** The analysis, code and documents in this
repository were produced in collaboration with Claude (Anthropic), working from
my direction and review. The commit trailers record this. I mention it plainly
because the right response to it is scepticism about verification, and the
answer to that scepticism is the audit trail rather than a disclaimer: every
external number is checked against a primary source, every prediction was
committed before the data was pulled, and every reversal is preserved in the
history rather than tidied away.

---

## The two threads

### 1. Kinetic threshold model (`prion_model.py`, `run_analysis.py`)

Nucleated polymerization (Masel–Jansen–Nowak) with an analytic self-sustaining
threshold `x_crit`, plus a toxicity layer.

- Required PrP knockdown to push replication below self-sustaining: **~65–90%**,
  depending on which published anchor is used. 50% is on the wrong side.
- Survival is strongly non-linear in knockdown depth, with a vertical asymptote
  at `x_crit` — which is why every 50%-lowering study looks similar and all of
  them fail.
- **T0-1 killed three of its own proposed experiments.** Flux-driven and
  neuronal-load-driven toxicity turn out to be *nested* — flux is the
  fast-clearance limit — so no reversal-timing or knockdown-depth experiment can
  separate them. See `T0-1-findings.md`. The retracted proposals are left in the
  history deliberately.
- Next step is `tasks/T0-4-dose-response-joint-fit.md`: fit the published
  dose–response jointly instead of to one anchor, which is why the 65–90% range
  is that wide.

### 2. PRNP penetrance refresh (`t02_*.py`, `T0-2-*.md`)

Reproduces Minikel et al. 2016 (*Sci Transl Med* 8:322ra9), then refreshes the
control side against gnomAD v4.1.1.

**What this contributes, stated narrowly:**

- Per-variant penetrance estimates for P102L, A117V, D178N and E200K, which the
  2016 paper pooled into a single zero-count row and which the 2024 blog post
  does not compute.
- The observation that **the 2016 zeros were never surprising**: at v4's
  frequencies, ExAC expected 0.17 / 0.08 / 1.08 alleles, with P(zero) = 0.85 /
  0.92 / 0.34. The original analysis was underpowered, not wrong.
- **No ancestry-matched control cohort exists for M232R or V180I in any current
  gnomAD release.** gnomAD v4 dropped subcontinental resolution; v2's `jpn` group
  contains 76 individuals against the 663 the 2016 analysis assigned to JPT.
  Since ExAC is nested inside v2, those ~587 people are still there, classified
  differently. The published estimates for those two variants rest on an ancestry
  assignment gnomAD's own inference no longer endorses.

**What it does not contribute:** the allele counts (see above), any
reclassification, or any defensible point estimate of E200K penetrance.

---

## Reading order

| file | what it is |
|---|---|
| `T0-2-results.md` | the penetrance results, with all caveats |
| `T0-2-reproduction.md` | reproducing 2016 against paper-printed numbers |
| `T0-2-power-check.md` | can more controls help? (asked before pulling data) |
| `T0-2-prediction-prereg.md` + `-amendment-1/-2` | predictions, committed pre-pull |
| `T0-2-gate-decision.md` | where the ancestry gate killed the primary estimand |
| `T0-1-findings.md` | the kinetic model, and three retracted discriminators |
| `RESEARCH_PLAN.md` | field overview and entry points |

## Method notes worth knowing

- The 2016 confidence interval is the **Kirov 2014 opposite-corner** construction,
  not a delta-method interval — log-width is additive in the component widths.
  Verified from source and against the paper's Methods.
- The reproduction is checked against **numbers printed in the paper**, not
  against the original code's own output. The estimator was committed *before*
  those numbers were retrieved, so it could not be tuned to match; `git log` is
  the evidence.
- The original R was also run (r-base 4.5.3 + CRAN `binom`) and agrees to ~1e-11.

## Run

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib pandas requests
python run_analysis.py            # kinetic model
python t02_power_check.py         # pre-pull power check
python t02_penetrance.py          # penetrance, from cached gnomAD responses
python t02_figure.py
```

Cached gnomAD responses are in `data/` so the analysis reruns without querying.
gnomAD data is released without restriction. The 2016 supplementary tables are
**not** redistributed here — that repository carries no licence and its
manuscript is posted "for personal use, not for redistribution." Clone it from
[github.com/ericminikel/prnp_penetrance](https://github.com/ericminikel/prnp_penetrance)
if you want to rerun the reproduction step.

## Key sources

- Minikel et al., *Sci Transl Med* 8:322ra9 (2016) — PRNP penetrance.
- Minikel, cureffi.org (3 April 2024) — PRNP in gnomAD v4.
- Mallucci et al., *Science* 302:871 (2003); *Neuron* 53:325 (2007) — reversal on
  neuronal PrP depletion.
- Masel, Jansen & Nowak, *Biophys Chem* 77:139 (1999) — nucleated polymerization.
- Büeler et al. (1993, 1994) — `Prnp-/-` resistance; `Prnp+/-` extended incubation.
- An, Davis et al., *Nat Med* (2025) — base editing, 50% PrP reduction, +52% lifespan.
- Gentile et al., *Nucleic Acids Res* (2026) — divalent siRNA, 2.7× survival at
  49% residual PrP.

Citations carried over from earlier stages are listed in `CITATIONS.md` with
their verification status. Anything not personally verified is marked.

## Caveat

Lumped-parameter toy model on the kinetics side; a control-frequency estimator
with uncharacterized systematic error on the genetics side. Both are hypothesis
generators. Neither is evidence in itself, and neither is clinical guidance.
