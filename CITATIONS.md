# Citation verification (V6)

Per CLAUDE.md integrity rule 4: claims are either fetched and quoted, or marked
`[UNVERIFIED]`. The design-critical facts are checked here before the criteria
that depend on them are written.

**Deliberately NOT recorded here yet:** the Mallucci 2007 electrophysiological
recovery *timecourse* (the held-out latency target). It is fetched and recorded
only in Step 6, *after* the power-check prediction intervals are committed, so
the pre-registered latency criterion cannot be tuned to the observed value. See
`archive/step5_power_check.py` / `docs/T0-1-findings.md`.

---

## Mallucci 2003 — the design premise [VERIFIED]

Mallucci G, et al. "Depleting neuronal PrP in prion infection prevents disease
and reverses spongiosis." *Science* 302(5646):871–874, 2003. PMID **14593181**.

Abstract, verbatim (fetched from PubMed 14593181, 2026-08-02):

> "We found that depleting endogenous neuronal PrPc in mice with established
> neuroinvasive prion infection reversed early spongiform change and prevented
> neuronal loss and progression to clinical disease. This occurred despite the
> accumulation of extraneuronal PrPSc to levels seen in terminally ill wild-type
> animals. Thus, the propagation of nonneuronal PrPSc is not pathogenic, but
> arresting the continued conversion of PrPc to PrPSc within neurons during
> scrapie infection prevents prion neurotoxicity."

This is the load-bearing premise for the whole of T0-1, and it is confirmed at
primary source:

- **Depletion was neuronal**, and **extraneuronal PrP-Sc kept accumulating to
  terminal levels** while the animal was rescued. → exactly the dissociation the
  single-compartment model cannot express, and the reason for the neuronal/glial
  split.
- "arresting the continued **conversion** ... within neurons ... prevents
  neurotoxicity" → direct support for the **flux** hypothesis over standing
  load, and the observation that **kills `load_total`** (total burden rose while
  toxicity reversed).
- It does **not** by itself separate `flux` from `load_neuronal`: both predict
  recovery, because neuronal load also falls once neuronal conversion stops.
  That separation is a *timing* question, handled by the sealed 2007 target.

Neuron-specificity of the Cre driver (NFH-Cre) is not in the abstract; the model
line is the standard NFH-Cre/tg37 conditional knockout, described as neuronal in
the Mallucci review (PMC2838383) and in CLAUDE.md's domain facts. Treated as
established.

## Prion mouse natural history — training anchors [PARTIALLY VERIFIED]

From the Mallucci reversibility review (PMC2838383, fetched 2026-08-02), for the
**tg37** line used in those studies (PrP-overexpressing, faster incubation):

> "Mice ... begin to show the first subtle clinical signs around 8 weeks
> post-infection" and "animals die at around 13 weeks post-inoculation."

→ earliest behavioural deficit at **~62% of incubation** (8/13 wk), consistent
with the general statement that early burrowing/nest-building deficits precede
frank illness by weeks.

**Caveat on transfer:** the kinetic model is calibrated to **wild-type RML**
(terminal 150 d, emergent onset 76 d = **51% of course**), not tg37 (~91 d
terminal). The *fraction* transfers better than absolute days. The model's 51%
onset falls just below the tg37 62% figure and inside the "early behavioural
deficit at ~50–60% of course" band. So interpreting the model's `D_sx` crossing
as the **early behavioural (burrowing) deficit** — not frank clinical onset — is
supported. See Step 4.

A precise wild-type-RML burrowing-deficit day, and a clean frank-onset→terminal
interval for wild-type RML specifically, were **not** located in an open source
and are marked `[UNVERIFIED]`; Step 4 states how each dependent choice is
loosened accordingly.

---

Sources consulted: PubMed 14593181; PMC2838383 (Mallucci review); cureffi.org
2013-04-08. Science/Neuron full texts are paywalled (HTTP 403).
