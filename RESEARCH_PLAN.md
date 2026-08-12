# Getting into prion reversal research: what's settled, what's open, what you can do

> **Status: partly superseded.** This was the entry document, written before the
> analyses in `docs/`. Two of its claims have since been overtaken:
> - §2's "roughly 65–90% knockdown required" is a model output whose defining
>   feature (`x_crit` = 33% residual PrP) sits below every published data point
>   (lowest = 49%). See [`docs/T0-4-structural-test.md`](docs/T0-4-structural-test.md).
> - §6's experimental design is superseded by [`EXPERIMENT.md`](EXPERIMENT.md),
>   which corrects a lifespan confound in the proposed 15% arm.
>
> The field overview in §1–§5 is still a fair on-ramp, which is why it is kept.

## 1. The honest state of the question

**Reversal is already proven — in mice, early, by removing the substrate.**
Mallucci and colleagues (Science 2003; Neuron 2007) infected mice with scrapie,
let spongiform pathology, synaptic dysfunction and behavioural deficits develop,
then switched off neuronal PrP-C with Cre recombinase. Spongiosis reversed.
Cognitive and electrophysiological deficits recovered. The animals lived a normal
lifespan — *even though PrP-Sc kept piling up outside neurons*. That last detail is
the conceptual key: the plaque is not the poison. Ongoing conversion inside neurons is.

**What is not proven** is reversal after neurons are actually gone, and reversal
in humans with an achievable drug. Late-stage sporadic CJD patients have lost
tissue; nothing on the horizon regrows it.

So the field's real question is not "reversible or not" but **how deep and how early**.

## 2. The one number the whole programme turns on

PrP-C is the obligate substrate. `Prnp-/-` mice are healthy and totally resistant.
So lowering PrP is a target validated by nature — unusually strong for neurodegeneration.

But ~50% lowering is not enough. Every modality converges on the same disappointing
plateau:

| Intervention | Residual PrP | Survival effect |
|---|---|---|
| `Prnp+/-` heterozygote | ~50% | ~2x incubation |
| PrP-lowering ASO, chronic early | ~50% | up to ~3x |
| AAV base editing (R37X), humanized mice | ~50% | +52% lifespan |
| Divalent siRNA, chronic presymptomatic | 49% | 2.7x |
| Divalent siRNA, single dose at symptom onset | 49% | +64% |
| Divalent siRNA, single 348 µg dose | **17%** | (potency study, not yet a survival study) |

Every one of these animals still died. As Gentile et al. put it in the 2026 NAR paper:
halting or indefinitely delaying disease will require deeper than 50% lowering.

The accompanying model (`prion_model.py`) makes that intuition quantitative. Under
nucleated polymerization there is an exact critical substrate concentration
`x_crit` below which polymer creation cannot outrun clearance, and prion load
*declines*:

```
x_crit = [ a(a + b(2n-1)) + b² n(n-1) ] / (b · β)
```

Fitting β to the published dose-response puts the required knockdown at
**roughly 65–90%**, depending on which mouse anchor you trust. 50% sits on the
wrong side of that line. 83% (what di-siRNA already achieves) sits on the right side.

That is also, almost certainly, why Ionis reopened the PrProfile trial in March 2026
to test a third, higher dosing regimen — the public read is that ION717 looked safe
enough but did not lower PrP as deeply as hoped.

## 3. The open questions worth your time

Ranked by (importance x tractability for someone starting out):

1. **Where exactly is `x_crit`?** No one has run the depth-titration survival
   experiment. Every published survival study clusters at ~50%. A study at
   40 / 60 / 75 / 85% residual PrP would locate the threshold directly. This is the
   single most decision-relevant unrun experiment in the field.
2. **Does deep knockdown *clear* prions, or just stall them?** Serial-passage or
   RT-QuIC seeding-activity measurement in deeply-lowered animals over time
   distinguishes "load falling" from "load frozen."
3. **How much of the late-stage deficit is neuron loss vs. reversible dysfunction?**
   Determines whether symptomatic patients can recover function, not just stop declining.
4. **Biomarker-triggered prevention.** For carriers of high-penetrance PRNP variants,
   what CSF/plasma signal (RT-QuIC seeding, NfL, PrP levels) fires early enough to
   treat inside the reversible window? Regulatory path for prevention trials
   depends entirely on this.
5. **Strain-independence.** Small molecules failed largely through strain adaptation.
   Substrate lowering should be strain-agnostic — but that needs testing across
   sCJD MM1, E200K, vCJD, and CWD-derived isolates.

## 4. What you can actually do, by tier

### Tier 0 — this week, laptop only, no permission needed

- **Reproduce and break the model in this repo.** Change the toxicity assumption
  from flux-driven to load-driven and see whether it can still reproduce
  Mallucci's reversal. If it can't, that's evidence the flux hypothesis is
  load-bearing. Publishable as a short methods note if done rigorously.
- **Reproduce Minikel 2016 penetrance analysis.** The entire thing is open source
  at `github.com/ericminikel/prnp_penetrance`. Rerun it against current gnomAD
  (v4 has ~800k exomes vs the 60k ExAC they used). Penetrance estimates for
  intermediate PRNP variants will have tightened substantially. This is a real,
  citable contribution and requires nothing but R/Python and patience.
- **Mine cryo-EM fibril structures.** Ex vivo prion fibril structures (RML, ME7,
  263K, human isolates) are in the PDB. Ask a structural question nobody has:
  does the PIRIBS fold explain why certain PRNP polymorphisms (129M/V, 219E)
  are protective? Molecular dynamics on the fold with and without the polymorphism.
- **Read `cureffi.org` end to end.** Minikel has been blogging the field openly
  for over a decade, including negative results, failed compounds, and regulatory
  strategy. It is the best on-ramp that exists for any disease area.

### Tier 1 — months, still cheap

- **ASO/siRNA target-site design against PRNP.** Sequence walk, off-target
  screening against the transcriptome, secondary-structure accessibility.
  Purely computational; feeds directly into a real pipeline.
- **Reanalyse public transcriptomics.** GEO has prion-infected mouse brain
  time-courses. Question worth asking: which transcriptional signatures reverse
  after PrP knockdown and which do not? That is a molecular readout of the
  reversibility window.
- **Build the depth-titration power calculation.** Given mouse-to-mouse variance
  in incubation time, how many animals at how many dose levels do you need to
  localise `x_crit` to ±10%? Hand that to a lab and you've made yourself useful.

### Tier 2 — needs an institution

Non-infectious work first: recombinant PrP expression, thermodynamic stability
assays of mutants, RT-QuIC with recombinant substrate. Cell models (prion-infected
N2a) are next. Live rodent infection studies are the far end.

**Biosafety is not a formality here.** Émilie Jaumain, an INRAE technician, pricked
her thumb through two gloves in 2010 handling prion-infected humanized mouse brain
and died of vCJD in 2019, aged 33 — the first pathologically confirmed occupational
prion transmission. A second French case led five institutions to impose a
moratorium in 2021, and a Spanish prion researcher has also died of CJD. The Broad
team could not even complete a follow-up survival study because their institution
banned human prion challenges. Prions resist autoclaving, formalin, and most
disinfectants. Do not improvise this. Ever.

## 5. Where to plug in

- **Broad Institute Prion Alliance / Vallabh–Minikel lab** — unusually open group.
  They published their IND filing publicly, which almost no one does.
- **UCL MRC Prion Unit** (Collinge), **UCSF IND** (Prusiner lab lineage),
  **NIH Rocky Mountain Labs**, **Aguzzi lab, Zurich**, **Khvorova lab, UMass Chan**
  (oligonucleotide chemistry).
- **CJD Foundation** and **Prion Alliance** — funding and patient-community access.
- **Prion20XX** annual conference — small field, genuinely approachable.
- The divalent siRNA trial began enrolling in April 2026 via **NeuroNEXT** (NINDS).

## 6. The experiment this model says to run

> **Superseded by [`EXPERIMENT.md`](EXPERIMENT.md).** Two corrections since this
> was written:
> - **The 15% arm below is uninterpretable.** If the no-threshold structure is
>   right it predicts ~787 days there, exceeding mouse lifespan — animals die of
>   other causes either way, so a negative result means nothing. The
>   discriminating window is **40–25% residual PrP**.
> - **"The model predicts a sharp inflection" is not a prediction the data
>   support.** Whether an inflection exists is exactly what is in question; a
>   threshold-free structure taken from Sandberg 2011 fits the published anchors
>   adequately and has no inflection at all.
>
> Kept unedited below because the reasoning is sound and the record of what it
> got wrong is the useful part.

If you can persuade one lab to do one thing:

> Infect wild-type mice with RML. Titrate divalent siRNA dose to achieve
> **stable residual PrP of ~60%, ~40%, ~25%, and ~15%**, chronic dosing from
> before inoculation. Measure survival and, at matched timepoints, brain prion
> seeding activity by RT-QuIC.
>
> The model predicts a sharp inflection, not a smooth curve: somewhere between
> 25% and 40% residual PrP, survival should stop scaling and start diverging, and
> seeding activity should begin *falling* rather than plateauing.
>
> If that inflection exists, the field's target product profile changes from
> "lower PrP as much as tolerable" to "clear a specific, knowable bar." If there
> is no inflection — if survival scales smoothly all the way down — the
> nucleated-polymerization framing is wrong and something else limits propagation.

Either result is worth knowing. That is what a good experiment looks like.

### A second readout also discriminates the toxicity mechanism (added after T0-1)

T0-1 (`T0-1-findings.md`) falsified the standing-*total*-load toxicity model but
found that two hypotheses survive and are **nested**: neuronal conversion
**flux** is the fast-clearance limit of neuronal standing **load**. No
reversal-*timing* or reversal-*depth* experiment can separate them — they differ
only in how fast neuronal PrP-Sc clears, and in the fast limit the models are
identical. (A depth-titration with a behavioural readout, which an earlier draft
of this section proposed as the discriminator, does **not** work for this reason;
that proposal is withdrawn.)

The one measurement that separates them is a **direct neuronal PrP-Sc clearance
timecourse**, run as the neuronal analogue of Mallucci 2003. Take the mice in
which neuronal PrP has been depleted mid-infection (the Mallucci reversal
cohort), and in the window where behaviour recovers, measure **neuronal** PrP-Sc
— by RT-QuIC seeding activity or PrP-Sc immunostaining on neurons, distinct from
the extraneuronal compartment Mallucci already assayed:

> - Behaviour recovers **while neuronal PrP-Sc is still present** (has not yet
>   cleared) → **flux**: toxicity tracks ongoing conversion, not the standing
>   species. This is the intraneuronal version of the 2003 result that killed
>   total-load.
> - Behavioural recovery **tracks the fall** of neuronal PrP-Sc → **neuronal
>   load**.

The depth titration in §6 remains valuable — it locates `x_crit` and tests the
sharp-inflection prediction — but for the toxicity *mechanism* the readout that
matters is neuronal PrP-Sc kinetics versus behaviour, not knockdown depth.

---

*Sources for every number above are cited in the accompanying README. The model is
a lumped-parameter toy: use it to generate hypotheses and design experiments, never
as evidence in itself.*
