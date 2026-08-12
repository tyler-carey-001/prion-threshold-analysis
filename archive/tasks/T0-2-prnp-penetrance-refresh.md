# T0-2 — Refresh PRNP penetrance estimates against current gnomAD

**Estimated session length:** one long session or two. Use plan mode first.
This is the highest-value task in Tier 0 and the one most likely to produce a
citable result.

## Objective

Reproduce the Bayesian penetrance analysis of Minikel et al., *Sci Transl Med*
8:322ra9 (2016), then re-run it against the current gnomAD release and report
how the penetrance estimates and their confidence intervals have moved.

## Why it matters

The 2016 paper used 60,706 ExAC exomes as population controls and showed that
PRNP missense variants previously reported as pathogenic are at least ~30× more
common in the population than genetic prion disease prevalence can explain —
meaning several "100% penetrant" variants actually confer low lifetime risk.
Some people have been told they carry a death sentence when their real risk may
be around 1% or lower.

gnomAD v4.1 spans 730,947 exomes and 76,215 genomes — roughly a 12-fold increase
in exomes over ExAC. For rare variants, confidence intervals scale steeply with
control cohort size. Several intermediate-penetrance variants that sat in a wide
uncertain band in 2016 should now resolve. Nobody appears to have published the
refresh.

**Verify the current gnomAD version yourself before starting** — v4.1 is the
latest I can confirm, but a newer release may exist.

## Method

1. **Reproduce before extending.** Clone `github.com/ericminikel/prnp_penetrance`
   and run the original analysis to completion first. Confirm you recover the
   published numbers. If you cannot, stop and report the discrepancy — do not
   proceed to new data on top of a broken reproduction.

2. **Understand the estimator before touching it.** The penetrance calculation
   is Bayes' theorem applied to allele frequencies in a case cohort versus a
   population cohort. Write out the derivation in your own words in
   `T0-2-method.md`, including every assumption:
   - assumed genetic prion disease prevalence in the general population;
   - assumption that population controls are below the age of onset;
   - assumption that the case cohort's ascertainment is unbiased.
   Each of these is a lever that can swing the answer. Say so.

3. **Pull current PRNP data.** Use the gnomAD GraphQL API
   (`gnomad.broadinstitute.org/api`) or the downloads page. You need, per variant
   across the PRNP coding region: allele count, allele number, and the
   exome/genome/joint breakdown. Cache the raw response to disk so the analysis
   is reproducible without re-querying.

4. **Handle the known subtleties.** Document how you deal with each:
   - **Exomes vs genomes vs joint frequencies.** v4 provides merged "joint"
     fields; using them naively can double-count. Check coverage over PRNP in
     both call sets before choosing.
   - **Cohort composition.** gnomAD excludes severe paediatric disease cohorts
     but includes adult disease cohorts. Check whether any contributing cohort
     is neurodegeneration-related, which would bias controls toward cases.
   - **Age distribution.** The estimator assumes controls have not yet reached
     onset. gnomAD publishes age distributions; report what fraction of carriers
     are above the mean onset age for each variant.
   - **Codon 129.** M129V genotype modifies risk and determines the D178N
     phenotype. Report codon 129 frequencies alongside variant counts. If phase
     is unavailable, say so rather than assuming.
   - **Truncating variants.** PRNP is single-exon, so premature stops truncate
     rather than trigger NMD, and the 2016 paper found position-dependent
     effects with genuine loss-of-function alleles in healthy older individuals.
     Keep truncating variants in a separate stratum.

5. **Re-estimate and compare.** For each variant analysed in 2016, produce a
   table: 2016 point estimate and CI, current point estimate and CI, and the
   change. Flag any variant whose classification would plausibly shift.

## Deliverables

- `T0-2-method.md` — the derivation and every assumption, written before results
- `data/gnomad_prnp_raw.json` — cached query response with retrieval date and
  gnomAD version recorded
- `t02_penetrance.py`
- `T0-2-results.md` — the comparison table plus a plain-language summary
- `fig_t02_penetrance_shift.png` — 2016 vs current estimates with CIs

## Acceptance criteria

- The original 2016 analysis is reproduced to within rounding before any new
  data is used. State the reproduced numbers explicitly.
- Every assumption in step 2 is listed with the direction it would bias results.
- Confidence intervals are reported everywhere. A point estimate without a CI is
  not an acceptable output for this task.
- The results file contains an explicit statement that these are population-level
  estimates and not clinical guidance for any individual.

## Traps

- **Do not reclassify variants.** Producing an updated estimate is the task.
  Declaring a variant benign or pathogenic is a clinical-genetics judgement that
  belongs to expert panels with access to segregation and functional data.
- Do not assume gnomAD v4 coordinates match the 2016 GRCh37 coordinates. v4 is
  native GRCh38 and lifted-over tracks are imperfect. Look up PRNP coordinates
  in the current build rather than reusing any from the old repo.
- If the original repo is in R and you would rather work in Python, port it —
  but reproduce it in its native language *first*, so a porting bug cannot be
  mistaken for a data-driven change.
