# T0-3 — Structural survey of prion fibrils and protective polymorphisms

**Estimated session length:** 2–3 hours. Most exploratory of the three; run last.

## Objective

Assemble every available ex vivo prion fibril structure, compare their folds,
and ask whether the known protective polymorphisms have a structural explanation.

## Why it matters

Cryo-EM structures of *ex vivo* prion fibrils (rather than recombinant amyloid)
appeared only recently, and they share a parallel in-register intermolecular
β-sheet architecture. Several PRNP polymorphisms are protective — codon 129
heterozygosity, E219K in Japanese populations, G127V which arose under selection
during the kuru epidemic. If protection maps onto specific structural positions
in the fibril core, that constrains which conversion step therapy should target.

This is a genuine open question with no guaranteed answer, which is why it goes
last. A well-documented negative result is still a useful contribution.

## Method

1. Search the PDB for ex vivo prion fibril structures — mouse-adapted strains
   (RML, ME7, 22L), hamster 263K, human isolates, and CWD where available.
   Build `structures.csv` with PDB ID, source organism, strain, resolution,
   ordered residue range, and citation. Verify each entry against the actual PDB
   entry rather than trusting a search summary.

2. Align the ordered cores across structures. Report which residues are ordered
   in all of them versus strain-specific.

3. Map onto the aligned cores: codon 129, E219K, G127V, and the pathogenic
   variants E200K, D178N, P102L. For each, record whether it sits in the ordered
   core, at an interface, or in a disordered region.

4. State the structural hypothesis, if any, and what would test it. If the
   positions show no coherent pattern, say that.

## Deliverables

- `structures.csv` with verified metadata
- `t03_structure_map.py`
- `fig_t03_variant_positions.png` — variant positions on the aligned fibril core
- `T0-3-findings.md`

## Acceptance criteria

- Every PDB ID is verified against the live entry; no IDs from memory.
- Residue numbering convention is stated explicitly (human PRNP numbering,
  and the offset used for non-human structures).
- Findings file distinguishes "observed in the structures" from "inferred".

## Traps

- Recombinant PrP amyloid fibrils adopt folds that differ from ex vivo prions.
  Do not mix them into the same comparison without flagging it prominently.
- Residue numbering differs between species. Off-by-a-few errors here will
  silently invalidate the whole variant mapping. Check twice.
