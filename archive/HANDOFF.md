# HANDOFF — running Tier 0 in Claude Code

## Setup

```bash
mkdir prion-tier0 && cd prion-tier0
# drop the downloaded files in here, then:
mkdir -p tasks
mv T0-*.md tasks/
git init && git add -A && git commit -m "baseline model and task specs"
python -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib pandas requests
python run_analysis.py   # confirm the baseline reproduces before changing anything
claude
```

`CLAUDE.md` sits at the repo root and Claude Code loads it automatically every
session, so the domain facts and integrity rules persist without you re-pasting.

## The three sessions

**Run one task per session and clear context between them.** These tasks share
almost no working state, and a long context full of T0-1's parameter sweeps will
make T0-2 worse, not better.

Order matters: T0-1 first (self-contained, no external data, fast feedback),
then T0-2 (highest value, needs the most care), then T0-3 (most exploratory).

### Session 1

> Read CLAUDE.md and tasks/T0-1-model-falsification.md. Before writing any code,
> tell me your plan and specifically what numeric criteria you propose for the
> pre-registration — I want to approve those before you fit anything.

### Session 2

Start in **plan mode** (`shift+tab` twice, or launch with `--permission-mode plan`).
This one touches external data and has the most ways to go quietly wrong.

> Read CLAUDE.md and tasks/T0-2-prnp-penetrance-refresh.md. Work in plan mode.
> Start by checking what the current gnomAD release is and confirming the
> ericminikel/prnp_penetrance repo still runs. Don't pull any new data until
> you've reproduced the 2016 numbers and shown me that they match.

### Session 3

> Read CLAUDE.md and tasks/T0-3-fibril-structures.md. Build structures.csv first
> and show it to me before doing any analysis — I want to check the strain and
> resolution metadata myself.

## Keeping the agent honest

The failure mode for this kind of work isn't bad code, it's a plausible-looking
result that quietly encodes a wrong assumption. Concretely:

- **Check `git log` ordering.** The pre-registration commit must predate the
  results commit. If it doesn't, the pre-registration is worthless.
- **Ask "what did you fit and what did you hold out?"** after any model claim.
  If the answer is vague, the comparison is probably circular.
- **Spot-check one citation per session.** Fetch the paper and confirm the
  number. Fabricated effect sizes are the most damaging possible error here and
  the hardest to notice downstream.
- **Be suspicious of clean results.** If T0-2 reproduces 2016 exactly on the
  first try, ask to see the intermediate allele counts. Real reproductions
  usually have at least one wrinkle worth understanding.
- **Reject "I couldn't find a counterexample."** T0-1 asks for a sampled
  parameter sweep with a stated pass fraction, not a manual search.

## After Tier 0

If T0-2 produces a meaningful shift in any penetrance estimate, that is worth
writing up — and worth emailing to the Broad group before you do, since they
built the original analysis and will tell you quickly if you've made an error
they already know about. Minikel has published openly for over a decade,
including negative results and their own IND filing; this is an unusually
approachable field.

If T0-1 kills the flux-driven model, that's the more interesting outcome and
changes what the reversibility claims in the README can say. Update them.
