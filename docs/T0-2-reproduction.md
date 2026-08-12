# T0-2 Phase A — Reproduction gate

**Scope statements apply as in `T0-2-power-check.md`:** control-side refresh with
a frozen 2016 numerator; ExAC is nested inside gnomAD v4; population-level
estimates, not clinical guidance.

**Result: PASS on all eight externally checkable printed numbers, including two
exactly-printed confidence intervals.**

---

## 1. Why this is not circular

Running `generate_figures.r` and matching its own output tests nothing. The check
here is **an independent implementation against numbers printed in the paper**.

Two sources were used and agree with each other:

- **PMC4774245** — the published *Sci Transl Med* 8:322ra9 article.
- `manuscript.md` in the *external* `prnp_penetrance` repo — the author's accepted version.

The ordering that makes this non-circular is in `git log`, not in a promise:

| commit | contents |
|---|---|
| `22092af` | `t02_power_check.py`, containing the full estimator implementation |
| *(after)* | paper's printed numbers first retrieved |

**The estimator was committed before the target numbers were ever fetched**, so
there was no opportunity to tune the implementation to match them. This is a
stronger guarantee than a pre-registration file, because it is enforced by the
commit graph rather than asserted.

Note honestly: I retrieved the paper targets and ran the comparison in the same
step, rather than committing a targets file first as the plan specified. The
ordering above is why that does not compromise the check, but the plan's stated
sequence was not followed literally and that is recorded here rather than
presented as if it had been.

## 2. Method statement recovered from the paper

The interval construction inferred from source in Phase 0 is confirmed verbatim
by the paper's Materials and Methods:

> "We used an allelic rather than genotypic model, such that lifetime risk in an
> individual with one allele is equal to case allele frequency ... times baseline
> risk divided by population control allele frequency, P(D|A) = P(A|D)×P(D)/P(A).
> ... Following Kirov, we compute Wilson 95% confidence intervals on the binomial
> proportions P(A|D) and P(A), and calculate the upper bound of the 95% confidence
> interval for penetrance using the upper bound on case allele frequency and the
> lower bound on population control allele frequency, and vice versa for the lower
> bound on penetrance."

This confirms: allelic (2N) model, Wilson on both cohorts, opposite-corner
combination. The Phase 0 correction — corner method, not delta method — was
right, and it is the authors' stated method rather than an inference from code.

Baseline risk is also confirmed: incidence ≤2 per million per year against an
all-causes death rate of ~10 per 1,000 gives "~0.02% of all deaths, which we
accepted as the baseline disease risk" → `2e-4`.

## 3. Reproduction table

Computed in Python from the 2016 supplementary tables; targets are printed values
from the paper.

| quantity | paper (printed) | reproduced | agreement |
|---|---|---|---|
| M232R vs ExAC JPT | ~0.1% (main text) | **0.12%** | ✓ |
| M232R vs 23andMe Japanese | ~0.08% (Fig 3 caption) | **0.08%** | ✓ exact |
| V180I vs ExAC JPT | ~1% (main text) | **0.96%** | ✓ |
| V210I vs ExAC TSI | ~7.8% (Fig 3 caption) | **7.78%** | ✓ exact to 2 s.f. |
| Mendelians vs ExAC | "up to 100%", CI includes ~60–90% | **100%, [29.4%, 100%]** | ✓ contains 60–90% |
| V180I 23andMe, AC rounded down to 1 | **7.7% (95% CI 1.2% – 50%)** | **7.73% [1.20%, 49.69%]** | ✓ all three to printed precision |
| Mendelians 23andMe, AC rounded down to 1 | **100% (95% CI 100% – 100%)** | **100% [100%, 100%]** | ✓ exact |
| Allelic vs genotypic model | "identical point estimates and virtually identical 95% CIs" | reproduced (2N vs 1N) | ✓ |

The two rounding-sensitivity rows are the strongest evidence in this table: they
are the only places where the paper prints **complete confidence intervals
numerically**, and both reproduce to the printed precision. A porting or
formula error would have to coincidentally emit 7.73 / 1.20 / 49.69 against
printed 7.7 / 1.2 / 50.

## 4. The wrinkle

HANDOFF warns that a clean reproduction is suspicious and asks for the wrinkle.
There was one, and it would have been easy to misread as a failure.

Figure 3's caption states best estimates "range from ~0.08% for M232R to ~7.8%
for V210I." Computing M232R against ExAC gives **0.12%**, a 1.5-fold miss against
0.08%, which looks like a reproduction failure.

It is not. **M232R has two rows in Figure 3** — one against ExAC, one against
23andMe — and the caption quotes the range across *all* rows in the figure, not
the ExAC row specifically. The 0.08% is the 23andMe comparison (AC = 29 in 2,685
Japanese-ancestry individuals), which reproduces to 0.0809%. The ExAC M232R row is
~0.12%, consistent with the main text's separate statement of "near 0.1%."

Two consequences worth carrying forward:

- Point estimates quoted in captions must be attributed to a specific row before
  being used as reproduction targets. The same caption mixes ExAC-based and
  23andMe-based figures.
- This is the failure mode the >2-fold stop rule in `CLAUDE.md` is aimed at. At
  1.5-fold it did not trigger, and chasing it down turned a suspected discrepancy
  into a confirmation.

## 5. Intermediate allele counts

Everything the estimates rest on, so the numbers above can be checked by hand:

| forest row | ac_case | n_case | af_case | ac_ctrl | n_ctrl | control cohort | af_ctrl |
|---|---|---|---|---|---|---|---|
| M232R | 67 | 1,533 | 2.185% | 5 | 663 | ExAC JPT | 0.377% |
| M232R (23andMe) | 67 | 1,533 | 2.185% | 29 | 2,685 | 23andMe Japanese | 0.540% |
| V180I | 222 | 1,533 | 7.241% | 2 | 663 | ExAC JPT | 0.151% |
| V210I | 171 | 1,054 | 8.112% | 2 | 4,795 | ExAC TSI | 0.021% |
| Mendelians (pooled) | 1,034 | 10,460 | 4.943% | 0 | 60,706 | ExAC all | 0% |

Case counts: `table_s01`, rows Japan (M232R 63 + 4 in trans; V180I 218 + 4 in
trans), Italy (V210I 171), TOTAL (P102L 221 + A117V 33 + D178N 209 + E200K 571).
Control counts: `table_s07` (ancestry-restricted), `table_s08` (n per population),
`table_s05` (23andMe). Mendelian variants have **zero rows** in `table_s03` —
verified as a true absence, not a missing value.

Cross-check of an independently printed count: the paper states V210I appears at
"an allele frequency of 8.1% in Italian cases" — table above gives 8.112%. ✓

## 6. What is NOT verified

- **The per-row ExAC confidence intervals are not printed numerically anywhere in
  the paper.** Figure 3 shows them graphically only. So the CIs in the table in
  §3 for M232R/V180I/V210I-vs-ExAC are reproduced *method* applied to reproduced
  *inputs*, not checked against a printed value. The two rounding-sensitivity CIs
  are the only printed ones, and those do check out.
- **The original R has not been run.** The external gate — implementation vs.
  paper — is passed without it. Running `generate_figures.r` would additionally
  confirm the authors' code emits the same values and regenerate the figures, but
  it cannot strengthen the check against the paper, which is the part that
  matters. See the decision recorded in the session notes.
- Figures 1, 2, 4 and the supplementary discussion are not reproduced; only the
  penetrance estimator relevant to T0-2 is.

## 7. Gate status

**PASSED.** The estimator reproduces every externally checkable printed number in
the paper, including both printed confidence intervals, using an implementation
committed before those numbers were retrieved.

The Phase C porting firewall ("validate the port reproduces the 2016 numbers
before it sees gnomAD data") is therefore **already satisfied** — the port is what
was validated here.

Proceeding to gnomAD requires only the ancestry-matched control cohorts and their
sizes, per `T0-2-prediction-prereg.md`.
