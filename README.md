# Prion disease: what the published data can and cannot tell us

Independent, dry-lab computational work on two questions:

- **How deeply must PrP be lowered to stop prion replication?**
- **How penetrant are the PRNP variants that cause genetic prion disease?**

The short version: **on both questions, the published data turns out not to
support the answer the field is working from** — and in both cases the reason is
specific enough to say what measurement would settle it. Those limits are the
content here, not a new estimate. See [`EXPERIMENT.md`](EXPERIMENT.md) for where
they point.

---

## If you carry a PRNP variant, please read this first

**Nothing here is medical advice, and no number here describes your risk.**

1. **These are population averages, not personal risk** — computed from allele
   frequencies in databases, not from anyone's family, age, or other genotypes.
2. **The E200K figure of 61.4% in these files is not a reliable estimate of
   anything.** It rests on **13 alleles** and disappears if 6 were miscalled.
   Three systematic errors are identified; they push in **opposite** directions,
   at least one cannot be quantified, and together they exceed the statistical
   uncertainty. Eric Minikel — who wrote the original 2016 analysis — considers
   **~90%** more credible. This work does not overturn that and does not claim to.
   The two numbers do not even measure the same thing: his comes from families
   that came to medical attention and is closer to *remaining* risk at a given
   age; ours is risk from birth in an unselected population, expected to be lower.
3. **"Lifetime risk" here means risk from birth — not your remaining risk.** At 55
   and unaffected, you have already passed through part of the risk window and
   your remaining risk is *lower*. This analysis does not compute that, because
   doing it properly needs onset-age data it does not have.
4. **Case counts are frozen at 2016.** Only the control side was updated.
5. **Codon 129 matters and is not modelled here.** It modifies risk and determines
   whether D178N presents as fatal familial insomnia or as CJD.

**Better places to go:** [CJD Foundation](https://cjdfoundation.org) ·
[Prion Alliance](https://www.prionalliance.org) ·
[cureffi.org](https://www.cureffi.org) · and a genetic counsellor for anything
about your own result.

---

## Status

**Unreviewed, independent work. Not peer reviewed.** No result here has been
checked by anyone who works on prion disease professionally, and several depend
on literature values not yet verified against their primary sources. Read it as a
set of arguments about what the published data can support, not as findings.

Corrections and disagreement are welcome — open an issue.

**Prior work you should read first.** Minikel published the gnomAD v4 PRNP allele
counts in [April 2024](https://www.cureffi.org/2024/04/03/learn-prnp-gnomad-v4/).
His counts and ours agree exactly (P102L 2, D178N 1, E200K 13). **The allele
counts here are not a new finding.** That post was found *after* the analysis was
complete — a process failure, recorded as one.

**Written with AI assistance.** The analysis, code and documents were produced in
collaboration with Claude (Anthropic), from my direction and review; commit
trailers record it. The right response to that is scepticism about verification,
and the answer is the audit trail rather than a disclaimer: every external number
is checked against a primary source, every prediction was committed before the
data was pulled, and **every reversal is preserved in the history rather than
tidied away.** There are four documented retractions. They are the most useful
thing in the repository.

---

## The data

| source | what | where |
|---|---|---|
| **gnomAD v4.1.1** | PRNP variants, exomes + genomes, ancestry and age breakdowns, filter status, coverage | `data/gnomad_prnp_v4_*.json` (cached) |
| **gnomAD v2.1.1** | same, for the nesting check | `data/gnomad_prnp_v2_raw.json` |
| **Minikel 2016 supplementary tables** | 2016 case and ExAC control counts | **not redistributed** — see LICENSE |
| **Published mouse dose–response** | 7 anchors, residual PrP vs survival ratio | `prion_model.ANCHORS` |

Cached responses mean every analysis reruns without network access. gnomAD is
released without restriction. The 2016 tables are not redistributed because that
repository carries no licence and its manuscript is posted "for personal use, not
for redistribution" — clone
[ericminikel/prnp_penetrance](https://github.com/ericminikel/prnp_penetrance)
directly to rerun the reproduction step.

**Anchor values in `ANCHORS` are literature figures not yet verified against
primary sources.** Conclusions here are about the *structure* of the inference,
which is why they stand ahead of that verification — but no fitted number should
be reported until it is done.

## The analysis

### Thread 1 — kinetics: how deep must PrP lowering go?

A nucleated-polymerisation replication layer (Masel–Jansen–Nowak) with an
analytic threshold `x_crit`, plus a toxicity layer. Headline output: **~65–90%
knockdown required**, with 50% — where every published survival study clusters —
on the wrong side.

Three findings, in the order they were established:

1. **The toxicity question is ill-posed.** ([`T0-1-findings.md`](docs/T0-1-findings.md))
   Flux-driven and neuronal-load-driven toxicity are *nested* — flux is the
   fast-clearance limit — so no reversal-timing or knockdown-depth experiment
   separates them. Three proposed discriminators were killed, including two of my
   own. Only a direct neuronal PrP-Sc clearance measurement works.
2. **The range is real but irreducible from published data.**
   ([`T0-4-identifiability.md`](docs/T0-4-identifiability.md)) A profile likelihood
   gives **62–92%**, reproducing the reported range — so it is the genuine
   likelihood width, not an artifact of anchor-picking. But studies at the *same*
   PrP level disagree by **35%**, and removing the Tga20 overexpression arm — the
   point that was supposed to tighten everything — moves the interval ~3 points.
   Reaching ±5 points needs that scatter cut to ~11%.
3. **The threshold's existence is structure-dependent.**
   ([`T0-4-structural-test.md`](docs/T0-4-structural-test.md)) — see below.

### Thread 2 — genetics: how penetrant are these variants?

Reproduces Minikel et al. 2016 (*Sci Transl Med* 8:322ra9) against
**paper-printed numbers**, then refreshes the control side against gnomAD v4.1.1.
([`T0-2-results.md`](docs/T0-2-results.md))

What it contributes, stated narrowly:

- **Per-variant estimates** for P102L, A117V, D178N, E200K — the 2016 paper pooled
  them into one zero-count row and the 2024 blog post does not compute them.
- **The 2016 zeros were never surprising.** At v4's frequencies ExAC expected
  0.17 / 0.08 / 1.08 alleles, with P(zero) = 0.85 / 0.92 / 0.34. The original
  analysis was **underpowered, not wrong**, and this supplies power rather than
  correcting an error.
- **No ancestry-matched control cohort exists for M232R or V180I** in any current
  gnomAD release. v4 dropped subcontinental resolution; v2's `jpn` group has 76
  individuals against the 663 the 2016 analysis assigned to JPT. Since ExAC is
  nested in v2 those ~587 people are still there, classified differently — so the
  published estimates rest on an ancestry assignment gnomAD no longer endorses.

What it does **not** contribute: the allele counts, any reclassification, or any
defensible point estimate of E200K penetrance.

## What's up with the model

**The headline number may not describe a real quantity, and that is the main
result.**

The model puts the PrP dependence in the *replication* phase and treats the toxic
phase as PrP-independent. **Sandberg et al. 2011** (*Nature* 470:540) — already
cited in this repo before any of this was noticed — measured the kinetics
directly and found the opposite assignment:

> "In phase 1, prions propagate exponentially, **not rate limited by** cellular
> prion protein (PrP^C) concentration ... to reach a maximal prion titre, which is
> also independent of PrP^C concentration **over the range we studied**."
>
> "This is followed by a plateau phase (phase 2), **which determines time to
> clinical onset**, the duration of which is **inversely proportional to** PrP^C
> concentration."

Written as a model that is **one parameter with no threshold at all** —
`survival(x) = p₁ + (1−p₁)/x` — which reproduces `Prnp-/-` resistance and fits
the published anchors adequately.

**But the load-bearing argument is not model selection.** It is this:

> **`x_crit` calibrates to 33% residual PrP. The lowest data point in existence is
> 49%.** The two structures agree within ~2× everywhere data exists and diverge
> **13-fold** below 40%. Both are validated precisely where they cannot be told
> apart.

A formal comparison nominally favours the threshold-free model, but **92% of that
margin is a small-sample arithmetic penalty** and the fit term actually favours
the original model. That comparison is weak support, not a result.

**And Sandberg's hedge cuts both ways.** A nucleated-polymerisation threshold
*predicts* that phase-1 PrP-independence must break down as PrP approaches
`x_crit`. Sandberg measured independence across 50–800% and found it held — which
is entirely compatible with a threshold below 50%. These are not competing
theories with a winner. They are **the same unmeasured region seen from two
sides.**

**Honest position: the required knockdown depth is unknown.** Not lower, not
higher. Unknown.

## What a future experiment would clarify

Full specification in **[`EXPERIMENT.md`](EXPERIMENT.md)**. In brief: infect mice
with RML, hold PrP at ~60 / 40 / 30 / 25% in separate arms, and measure **prion
titre over time** as well as survival — one lab, one strain, one route, one
endpoint definition, individual animal times.

Three things fall out of the analysis that sharpen the standard proposal:

- **The discriminating window is 40–25% residual PrP.** Above 40% the structures
  are within 2× and indistinguishable.
- **A 15% arm is uninterpretable.** The no-threshold model predicts ~787 days
  there, exceeding mouse lifespan — animals die of other causes either way, so a
  negative result means nothing. The widely-proposed 15% arm is confounded by
  design.
- **Titre-vs-time, not survival.** The structures disagree about *which phase*
  lengthens; survival confounds both and discards exactly the discriminating
  information.

The structural question is **cheap** — n ≈ 10–15 per arm, since it distinguishes
0/n from n/n. Locating `x_crit` precisely is expensive and should come second.

Every outcome is publishable and every outcome changes what the field does next.

## Reading order

| file | what it is |
|---|---|
| [`EXPERIMENT.md`](EXPERIMENT.md) | **the payoff** — what to run and why |
| [`T0-4-structural-test.md`](docs/T0-4-structural-test.md) | does the threshold exist at all? |
| [`T0-4-identifiability.md`](docs/T0-4-identifiability.md) | why the range is 62–92% and what would narrow it |
| [`T0-2-results.md`](docs/T0-2-results.md) | penetrance results, with all caveats |
| [`T0-2-reproduction.md`](docs/T0-2-reproduction.md) | reproducing 2016 against paper-printed numbers |
| [`T0-1-findings.md`](docs/T0-1-findings.md) | the kinetic model; three retracted discriminators |
| `docs/T0-*-prereg*.md`, `docs/T0-2-power-check.md`, `docs/T0-2-gate-decision.md` | predictions and gates, committed before the data |
| [`archive/README.md`](archive/README.md) | **the four retractions**, and why the superseded material is kept |
| [`RESEARCH_PLAN.md`](RESEARCH_PLAN.md) | field overview and entry points |

## Method notes

- The 2016 confidence interval is the **Kirov 2014 opposite-corner** construction,
  not delta-method — log-width is additive in the component widths. Verified from
  source and against the paper's Methods.
- The reproduction is checked against **numbers printed in the paper**, not the
  original code's own output. The estimator was committed *before* those numbers
  were retrieved, so it could not be tuned to match; `git log` is the evidence.
- The original R was also run (r-base 4.5.3 + CRAN `binom`) and agrees to ~1e-11.
- Commit trailers carry `Claude-Session:` URLs referencing private sessions; they
  do not resolve publicly. Left in place — the history is this repository's main
  asset and rewriting it to tidy a dead link would cost more than the link.

## Run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# run from the repository root — scripts read data/ and write figures/
python run_analysis.py              # kinetic model + figures 1-4
python t04_identifiability.py       # profile likelihood on x_crit
python t04_structural_test.py       # threshold vs no-threshold
python t04_bootstrap.py             # stability audit of that comparison (~5 min)
python t02_power_check.py           # pre-pull power check
python t02_penetrance.py            # penetrance, from cached gnomAD responses
python t02_figure.py                # penetrance figure
```

All seven run offline from cached data and reproduce every number quoted here.

## Layout

```
README.md          EXPERIMENT.md      <- start here
prion_model.py     t0*.py             analysis code, run from root
docs/                                 write-ups, pre-registrations, reproduction record
data/                                 cached gnomAD API responses
figures/                              current figures
archive/                              superseded and retracted material, with a README
RESEARCH_PLAN.md                      field overview (partly superseded — see its banner)
CITATIONS.md                          citation verification status
```

## Key sources

- Minikel et al., *Sci Transl Med* 8:322ra9 (2016) — PRNP penetrance.
- Minikel, cureffi.org (3 Apr 2024) — PRNP in gnomAD v4.
- **Sandberg et al., *Nature* 470:540 (2011); *Nat Commun* 5:4347 (2014)** —
  two-phase kinetics. The structural test turns on this.
- Mallucci et al., *Science* 302:871 (2003); *Neuron* 53:325 (2007) — reversal on
  neuronal PrP depletion.
- Masel, Jansen & Nowak, *Biophys Chem* 77:139 (1999) — nucleated polymerisation.
- Büeler et al. (1993, 1994) — `Prnp-/-` resistance; `Prnp+/-` incubation.
- An, Davis et al., *Nat Med* (2025) — base editing, 50% reduction, +52% lifespan.
- Gentile et al., *Nucleic Acids Res* (2026) — divalent siRNA, 2.7× at 49% residual.

`CITATIONS.md` lists carried-over citations with verification status. Anything not
personally verified is marked.

## Caveat

A lumped-parameter toy model on the kinetics side; a control-frequency estimator
with uncharacterised systematic error on the genetics side. Both are hypothesis
generators. Neither is evidence in itself, and neither is clinical guidance.
