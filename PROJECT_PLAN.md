# Project plan

## Research question

How does the T790M gatekeeper mutation alter erlotinib binding to EGFR, and
what would that change imply for the dose required to maintain target
occupancy?

This is a self-directed learning project. The findings reproduce results
already established in the literature. The purpose is to work through the
argument from structure to dose using primary data, and to build the
computational skills the argument requires.

## Objectives

1. Characterise the erlotinib binding site in wild-type EGFR and locate the
   gatekeeper residue relative to the bound drug.
2. Compare the wild-type and T790M structures to describe what the mutation
   changes sterically.
3. Collect published binding affinities for erlotinib and osimertinib against
   wild-type and mutant receptor.
4. Build a receptor occupancy model that relates drug concentration, receptor
   turnover and downstream signal.
5. Use the model to ask what erlotinib dose would restore mutant occupancy to
   wild-type levels, and whether that dose is plausible.

## Phase 1 — Structural characterisation

- [x] Open wild-type EGFR + erlotinib structure (PDB 1M17)
- [x] Identify binding-site residues within 5 A of the ligand
- [x] Resolve the residue numbering offset (+24 vs clinical numbering)
- [x] Measure gatekeeper contacts with erlotinib (20 contacts, cutoff -2 A)
- [ ] Apply RCSB quality criteria to structure selection (resolution, ligand
      fit) and document the criteria used
- [ ] Consider replacing 1M17 with a higher-resolution structure if one meets
      the criteria better
- [ ] Obtain a T790M mutant structure; superimpose and report RMSD
- [ ] Obtain an osimertinib-bound structure; compare binding modes
- [ ] Produce three publication-quality figures

## Phase 2 — Parameter collection

- [ ] Published binding affinity (Kd, Ki or IC50) for erlotinib vs wild-type
- [ ] Same for erlotinib vs T790M
- [ ] Same for osimertinib vs both
- [ ] Receptor synthesis and degradation rates, if available in the literature
- [ ] Assemble a parameter table with full citations and assay conditions

Rule for this phase: every value is taken from a published measurement and
cited. No value is computed from structure. Assay conditions are recorded
because affinities measured by different methods are not directly comparable.

## Phase 3 — Model

- [ ] Working knowledge of R and the deSolve package
- [ ] Model 1: reversible binding only. Drug, free receptor, complex.
      Output: occupancy against concentration.
- [ ] Model 2: add receptor synthesis and degradation.
      Question: how does turnover change steady-state occupancy?
- [ ] Model 3: add a downstream signal driven by free receptor.
      Question: does 90% occupancy produce 90% loss of signal?
- [ ] Model 4: substitute mutant affinity values.
      Question: what dose restores wild-type occupancy, and is it tolerable?
- [ ] Sensitivity analysis on the least certain parameters

Each model is a separate script that runs on its own, so the progression is
visible and each step can be understood before the next is added.

## Phase 4 — Write-up

- [ ] README stating question, method, result and limitations
- [ ] Four to six page write-up
- [ ] questions.md — open questions raised during the work
- [ ] Optional: interactive version for teaching use

## Out of scope

Stated explicitly so the boundaries are clear:

- No molecular dynamics or free-energy calculation
- No prediction of binding affinity from structure
- No modelling of the covalent binding kinetics of osimertinib in the first
  version
- No pharmacokinetic model of absorption or distribution; drug concentration
  is treated as an input
- The ATP-affinity mechanism of T790M resistance is acknowledged but not
  modelled

## Success criteria

The project is complete when:

- The scripts run from a clean install and reproduce every figure
- Every parameter has a citation
- The limitations section names what the model cannot support
- The author can explain every line of code and every term in every equation
  without reference to notes

The last criterion is the binding one.

## Schedule

Approximately five hours per week alongside full-time teaching.

| Weeks | Focus |
|---|---|
| 1-2 | Phase 1 completion |
| 3 | Phase 2 |
| 4-5 | Phase 3, models 1-2 |
| 6-7 | Phase 3, models 3-4 |
| 8 | Phase 4 |

Realistic completion: eight to ten weeks.
