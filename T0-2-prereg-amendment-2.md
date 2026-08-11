# T0-2 — Pre-registration amendment 2

**Committed before any PRNP allele counts are retrieved from gnomAD.** Amends
`T0-2-prediction-prereg.md` and `T0-2-prereg-amendment-1.md`. Verify with
`git log` that this commit precedes any commit containing gnomAD PRNP data.

Covers: the composition problem in the surviving row, per-variant power, filter
handling and what counts as a resolved zero, and the v2-alongside-v4 pull.

---

## A4. The Mendelian row is unmatched, not matched — and that needs standardizing

`T0-2-gate-decision.md` §4 said the Mendelian row "never had an ancestry match to
degrade." That is true, and it was stated in a way that reads as a clean bill of
health. **The accurate version:** it is the one row where ancestry was never
controlled *in either epoch*. ExAC-all and gnomAD-v4-all are both unmatched
pools, and they are **differently composed** — v4 is UK Biobank-dominated, ExAC
was not. E200K carries founder enrichment (Libyan Jewish, Slovak), D178N also.
So the pooled control allele frequency is a composition-weighted average whose
weights changed between epochs.

The direction of that effect is not obvious, which is exactly why it is being
fixed in advance rather than interpreted afterwards. The single-variant rows got
this discipline in A1; applying it there but not here would let the asymmetry do
silent work.

### Rule: report two numbers per variant, both pre-specified

1. **Crude estimate** — v4 pooled control AF, directly comparable in construction
   to the 2016 ExAC-all figure. Carries the composition confound.
2. **Ancestry-standardized estimate** — v4 per-ancestry AFs directly standardized
   to **ExAC's** ancestry proportions:

   ```
   af_control_std = Σ_g  w_g^ExAC · af_g^v4        (g over shared ancestry groups)
   ```

   This isolates the sample-size effect from the composition effect: it answers
   "what would v4's control frequency be in a cohort composed like ExAC's."

**The weights are frozen here, computed pre-pull from `table_s08`** (see
`t02_power_check_pervariant.py`, which asserts they sum to the published 60,706):

| group | ExAC n | weight |
|---|---|---|
| nfe | 32,872 | 54.15% |
| sas | 8,411 | 13.86% |
| amr | 6,007 | 9.90% |
| afr | 5,102 | 8.40% |
| eas | 4,337 | 7.14% |
| fin | 3,977 | 6.55% |
| **total** | **60,706** | **100.00%** |

gnomAD's `asj`, `mid` and `remaining` have no counterpart in the 2016 labelling
and are excluded from the standardization; **the excluded fraction of v4 is
reported alongside every standardized estimate.** ExAC's `FIN` is kept separate
from `nfe`, matching gnomAD's own split.

Divergence between the crude and standardized estimates **is** the composition
effect, and is reported as such rather than absorbed into either number.

## A5. Per-variant power — checked before the pull, not after

Unpooling changes the variance structure. Per-variant case counts (`table_s01`
TOTAL) span 17-fold in `1/x_case`:

| variant | ac_case | 1/x_case | control alleles for parity | zero-persistence lower bound at n≈731k |
|---|---|---|---|---|
| P102L | 221 | 0.0045 | 221 | **70.5%** |
| A117V | 33 | **0.0303** | **33** | **8.6%** |
| D178N | 209 | 0.0048 | 209 | **66.4%** |
| E200K | 571 | 0.0018 | 571 | **100%** |

**The hypothesis that A117V may be case-dominated post-unpooling was checked and
does not hold in the plausible regime.** A117V becomes case-dominated only above
33 control alleles; at any realistic count for a genuine Mendelian prion variant
(0–10), `1/x_control ≥ 0.1` still exceeds `1/x_case = 0.030`, so all four
variants remain control-dominated and all four can be moved by gnomAD. Recording
this as a checked-and-rejected prediction rather than quietly dropping it.

**A real per-variant asymmetry does exist, and it is in the zeros.** If the zero
persists at v4 scale, it is a strong result for E200K (lower bound 100%), P102L
(70.5%) and D178N (66.4%), but a weak one for A117V (**8.6%**) — because 33 case
alleles cannot support a tight bound no matter how large the control cohort gets.

> **Pre-registered consequence:** a persisting zero for A117V is reported as
> **uninformative**, not as evidence of high penetrance. The other three are
> reported as informative lower bounds. This is fixed now so that A117V's weak
> bound cannot later be presented alongside the others as though comparable.

A117V's informativeness ceiling — the narrowest interval obtainable at *any*
control cohort size — is a 2.0-fold range, versus 1.2-fold for E200K.

## A6. Filter handling and what counts as a resolved zero

**This is the highest-risk item in the pull, and it lands exactly on the crux.** A
variant present in v4 but failing `AS_VQSR`, or flagged `AC0`, is indistinguishable
from a true absence if the query returns PASS-only records. The entire
reduced-scope deliverable turns on whether the zero resolves, so a false null
here would be the worst available outcome and the hardest to notice.

### Rules, fixed in advance

1. **Query filtered variants explicitly.** Retrieve every PRNP record regardless
   of filter status, and record the `filters` and `flags` arrays per variant per
   call set. A PASS-only query is not acceptable for this task.
2. **Report flags per variant**, verbatim, in the results table — including for
   variants that are absent, so "absent" is distinguishable from "present but
   filtered."
3. **Zero-resolution definition:**
   - **Resolved** — non-zero allele count in **PASS** records. Counts as new
     evidence.
   - **Filtered-only detection** — allele count zero in PASS but non-zero in
     filtered records (`AC0`, `AS_VQSR`, or other). Reported as its own category
     with the flags named. **Not** counted as resolved, and **not** counted as a
     persisting zero either.
   - **Zero persists** — no records in either, at stated coverage. Reported with
     the lower bound from A5, subject to A5's A117V carve-out.
4. **Coverage is reported before any zero is interpreted.** A zero at low
   callable coverage is not evidence of absence. Mean/median coverage over the
   PRNP coding region in the exome call set is reported per variant position.
5. **Exomes are primary** (ExAC was exomes, so this preserves the construction);
   genomes reported separately; **joint fields avoided** to prevent
   double-counting, per the original task spec.

## A7. Pull v2.1.1 counts alongside v4, to separate reprocessing from new evidence

ExAC ⊂ gnomAD v2 ⊂ (largely) gnomAD v4. If P102L or E200K shows non-zero in
**v2** where ExAC showed zero, the resolution is a **reprocessing / re-calling
change on largely the same samples**, not new evidence from new sequencing.

Those are different claims and only one of them is "12× more controls found
carriers." The distinction is cheap to check and is the first thing the Broad
group would ask about.

> **Rule:** every variant's row reports ExAC 2016 count, v2.1.1 count, and v4.1.1
> count, with the nesting stated. Any non-zero is attributed to one of:
> **(a) re-calling on overlapping samples** (appears in v2), or
> **(b) new samples** (appears in v4 but not v2). No resolution is described as
> new evidence without passing this check.

## A8. The 663 → 76 discrepancy is a reportable finding

ExAC is nested in gnomAD v2, so the ~587 individuals labelled Japanese by the
2016 analysis but not assigned `jpn` by gnomAD are still in the data — gnomAD's
ancestry inference simply classifies them otherwise.

**That means the published 2016 M232R and V180I penetrance estimates rest on a
population assignment that gnomAD itself no longer endorses.** This is a
limitation of the published analysis, independent of whether any refresh is
performed, and it was established without pulling any PRNP data.

> **Rule:** reported as a limitation paragraph in `T0-2-results.md`. Stated as a difference in ancestry
> inference between two methods, **not** as an error in the 2016 analysis —
> establishing which assignment is better is outside this task's scope and would
> require the individual-level data the project is not permitted to download.

## A9. Repurposed R validation

The R cross-check validated the four ExAC interval widths. Three of those belong
to rows whose tightening is now void under A1, so that part of its value has
lapsed. It remains the **port validation for the Mendelian computation** — the
`Mendelians_ExAC` width agreed to 2.9e-11, and the same Python code path produces
the per-variant estimates. Recorded so the R run is not later cited as validating
more than it does.
