# Archive — superseded material, kept deliberately

Nothing here is current. It is retained because **the record of what was tried and
retracted is more informative than a clean result**, and because deleting it would
make the git history harder to follow.

If you are reading the repository for its conclusions, you want the root
`README.md` and `EXPERIMENT.md`. Nothing in this directory should be cited.

| item | what it was | why it is here |
|---|---|---|
| `fig_t01_depth_discriminator.png`, `t01_depth_discriminator.py` | A proposed experiment to separate flux-driven from neuronal-load-driven toxicity by knockdown depth. | **Retracted.** The two hypotheses are *nested* — flux is the fast-clearance limit of neuronal load — so no depth titration can separate them. See `docs/T0-1-findings.md`. The figure depicts a discriminator that does not exist. |
| `STEP5_FINDINGS.md`, `step5_power_check.py` | A power check on whether recovery *latency* could discriminate the toxicity modes. | Superseded by the nesting result, which showed the question was ill-posed rather than underpowered. |
| `HANDOFF.md` | Instructions for running the Tier-0 tasks as a sequence of agent sessions. | Operational scaffolding. Those sessions are done; the instructions describe a workflow, not a finding. |
| `tasks/` | The original task specifications (T0-1 model falsification, T0-2 penetrance refresh, T0-3 fibril structures, T0-4 dose–response fit). | T0-1 and T0-2 are complete and written up in `docs/`. T0-3 was never run and was displaced as not decision-relevant. T0-4's central premise was falsified by the check that gated it — its spec carries a banner saying so. |

## The retractions, in order

Four claims were made and then withdrawn during this work. They are listed here
because a reader assessing reliability should be able to find them quickly.

1. **Recovery latency discriminates the toxicity modes** — withdrawn. Latency is
   largely non-identifiable from the training data (`STEP5_FINDINGS.md`).
2. **Knockdown depth discriminates the toxicity modes** — withdrawn. The
   hypotheses are nested (`docs/T0-1-findings.md`).
3. **The Tga20 overexpression arm will tighten `x_crit`** — withdrawn. Removing it
   entirely moves the interval ~3 percentage points
   (`docs/T0-4-identifiability.md`).
4. **ΔAICc = 64 shows the threshold-free structure fits better** — withdrawn. 92%
   of that margin is a small-sample arithmetic penalty and the fit term favours
   the original model (`docs/T0-4-structural-test.md` §3).

Several smaller corrections are recorded inline in the documents that made them,
including a mis-attributed penetrance figure, a wrong bias direction, an
overstated robustness margin, and an unfair model comparison.
