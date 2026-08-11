# T0-2 — Pre-registration amendment 1

**Committed before any gnomAD data is retrieved.** Amends
`T0-2-prediction-prereg.md`. Verify with `git log` that this commit precedes any
commit containing gnomAD data.

Three changes: a constraint on when tightening is interpretable at all, an
unpooling rule for the Mendelian row, and a decision rule for the
matched-vs-larger trade. All three are declared before the data exists so that
none of them can become a data-driven analysis choice.

---

## A1. Tightening is only interpretable within a fixed estimand

**The original prereg is wrong on this point and is superseded.** It claimed
interval tightening "does not require the allele frequency to be stable between
cohorts, only that more controls were sequenced," and used that to promote
tightening to primary estimand while demoting point movement to
secondary-and-confounded.

That justification fails in exactly the operation about to be performed. When the
control cohort **broadens** rather than merely growing, two things happen at once:

1. `n_control` rises steeply → `W_control` collapses → the interval tightens.
2. A founder-enriched variant measured across a broader ancestry group gets
   **diluted** → `af_control` falls → penetrance rises.

Both are artifacts of the denominator changing identity, not of learning anything
about the 2016 estimand. And the tightening artifact is the more dangerous of the
two, because movement was flagged as confounded while tightening was declared
clean.

### Worked demonstration — M232R

| control cohort | af_ctrl | penetrance | 95% CI | width | vs 2016 |
|---|---|---|---|---|---|
| 2016 ExAC JPT (matched) | 0.377% | 0.116% | [0.039, 0.343] | 2.169 | — |
| gnomAD EAS-wide, AF diluted to ~0.08% | 0.080% | **0.546%** | [0.305, 0.976] | 1.162 | **46.5% tighter** |
| gnomAD jpn-like, AF preserved | 0.378% | 0.116% | [0.078, 0.172] | 0.790 | 63.6% tighter |

The middle row is simultaneously **46.5% tighter and 4.7× higher**. It clears
P1's ≥25% bar and would likely clear P2's ±15 pp band. **Both pre-registered
predictions would pass, for the wrong reason**, and the output would be a
confidently narrow interval around a quantity that is not the 2016 estimand —
"lifetime risk in a Japanese heterozygote" quietly replaced by "lifetime risk in
an East Asian heterozygote," which is a different number about a different
population.

### Constraint

> **Tightening may be reported as a result only for rows where the control
> cohort's ancestry match to 2016 is preserved.** Where the match degrades,
> tightening is reported as **uninterpretable**, not as a hit, and the row is
> excluded from the P1/P2 hit rate rather than counted.

Match preservation is judged and stated per row before estimates are computed,
using the criteria in A3. A row cannot be rescued by noting that its interval got
narrower.

### Consequence for P1 and P2

P1 and P2 are hereby **conditional on match preservation**. A row whose match
degrades does not pass, fail, or count. If no row preserves its match, P1 and P2
are undefined and the primary estimand does not exist — see A3.

## A2. The pooled Mendelian row is unpooled

Pooling P102L, A117V, D178N and E200K was correct in 2016 because all four were
zero in ExAC and the row expressed a single upper-bound statement: *collectively
rare enough to be consistent with complete penetrance*.

Once any one of them resolves off zero, the pooled quantity becomes a
case-count-weighted blend across four variants with genuinely different
penetrances — weighted 571 : 221 : 209 : 33, so dominated by E200K — and
corresponds to no biological quantity at all.

### Rule

> **Report per-variant estimates for P102L, A117V, D178N and E200K separately**,
> using the per-variant case counts already present in `table_s01`:

| variant | ac_case | n_case | af_case |
|---|---|---|---|
| P102L | 221 | 10,460 | 1.056% |
| A117V | 33 | 10,460 | 0.158% |
| D178N | 209 | 10,460 | 0.999% |
| E200K | 571 | 10,460 | 2.730% |

> The pooled row is retained **only** as a reproduction artifact for continuity
> with 2016, never as a refreshed result, and is labelled as such.

This is declared now, before the data is seen. Pulling first, observing E200K's
control count, and then deciding to unpool would be a data-driven analysis choice.
Selection remains 2016-based — these are the same four variants the paper
identified on prior-pathogenicity grounds — so the frozen-list principle is
intact. This is an amendment to *presentation granularity*, not to *inclusion*.

If unpooling proves impossible for a variant (e.g. no interpretable matched
control cohort exists), that variant is reported as unresolvable rather than
folded back into a pool.

### A per-variant zero is a result, not a null

At n ≈ 730k, a variant that remains at `ac_control = 0` yields a **tight lower
bound on high penetrance**, which is a substantive finding and is reported as
one. It must not be described as "no data" or omitted for lack of a point
estimate. Prediction P4 is amended accordingly: it is assessed per variant, and
"still zero" is recorded as an informative outcome for that variant rather than
as P4 failing.

## A3. Ancestry-group availability is a gating question, answered before the pull

The 2016 controls are single subpopulations: **JPT** (Japanese in Tokyo) for
M232R and V180I, **TSI** (Toscani in Italia) for V210I. Preserving the estimand
requires control cohorts at comparable resolution.

gnomAD v2 carried subcontinental breakdowns — EAS into `jpn`/`kor`/`oea`, NFE
into `seu`/`nwe`/`est`/`swe`/`bgr`. gnomAD v3 moved to broader genetic-ancestry
groups. **Whether v4 exposes subcontinental resolution is not assumed in either
direction and is checked before the full pull** — one query per row, and it
determines whether the primary estimand exists.

### Decision rule, fixed in advance

> **Prefer ancestry match over cohort size whenever the matched option clears the
> kill gate.** The larger cohort is the wrong answer if it costs the estimand.

Concretely, per row, in priority order:

1. **v4 subcontinental group matching the 2016 subpopulation**, if exposed → use it.
2. **v2 subcontinental group** (`jpn`, `seu`) at ~2× the ExAC cohort → use it.
   The power table gives 25% tightening at 2.0–2.3×, which clears the kill gate,
   so the matched-but-smaller path is viable on its own terms and is **preferred
   over a ~30× broadened cohort**.
3. **Broadened continental group only** → the row's tightening is reported as
   uninterpretable per A1. A movement estimate may still be reported, explicitly
   labelled as a different estimand (e.g. "East Asian–wide" rather than
   "Japanese"), never compared numerically to the 2016 value as though it were
   the same quantity.

Note that option 2 reintroduces the nesting caveat more strongly, since v2 is
wholly contained in v4 and largely overlaps ExAC. Whatever is used, the overlap
fraction is stated, and no row is framed as an independent replication.

### Amended kill gate

> Full scope if ≥2 of the frozen rows show projected ≥25% tightening at a
> multiplier from a **match-preserving** cohort (tier 1 or 2 above).
> Rows resolvable only via tier 3 do not count toward the gate.
