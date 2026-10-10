# Phase 1 specification — resistance shift predictor

This document fixes what the model predicts, in what units, against which
baselines, and what counts as success. Nothing in Phase 2 begins until the gate
at the end of this document is passed.

It extends the EGFR gatekeeper analysis in this repository from a single
mechanistic case study to a predictive model across kinases. The EGFR
erlotinib case remains the worked example the model is held to.

## The prediction

Given a protein kinase, a single point mutation in it, and a reversible
ATP-competitive inhibitor, predict how much potency that inhibitor loses
against the mutant compared with the reference form of the same kinase.

| | |
|---|---|
| Input 1 | Kinase identity, as a UniProt accession (EGFR is P00533) |
| Input 2 | One point mutation, in clinical numbering, e.g. `T790M` |
| Input 3 | Reference form the shift is measured against: wild type, or a named activating mutant such as `L858R` |
| Input 4 | Inhibitor, as a canonical SMILES string |
| Output 1 | Predicted log10 of the potency ratio, mutant over reference |
| Output 2 | A 90% prediction interval on that number |
| Output 3 | A resistant / not-resistant call at the 10-fold threshold, with a probability |

Worked example, from this repository's own data. Input: P00533, `T790M`,
reference `L858R`, erlotinib. The measured answer is Ki 41 nM against 4.7 nM,
a ratio of 8.7, so log10 = 0.94 (van Alderwerelt van Rosenburgh et al.,
Nat Commun 2022;13:6791). A model returning 0.94 with an interval of roughly
0.6 to 1.3 and a "not resistant at 10-fold, probability 0.4" call has answered
correctly.

The reference form is part of the input, not an assumption. T790M in a patient
almost always arises on top of an activating mutation, so the shift against
L858R is the clinically meaningful one; the shift against wild type is a
different number answering a different question.

## Scope

Every boundary below exists so that rows in the training table are physically
comparable. A model trained on incomparable rows learns the assay, not the
biology.

| In scope | Out of scope | Why the line is drawn here |
|---|---|---|
| Protein kinases | All other target classes | Kinases share a fold and an ATP site, so a mutation's position means the same thing across the family |
| Reversible ATP-competitive inhibitors | Covalent inhibitors (osimertinib, afatinib), allosteric inhibitors (EAI045) | A covalent drug's potency is governed by kinact/KI and depends on incubation time; it is not an equilibrium constant and cannot share a column with Ki |
| Single point mutations | Double and triple mutants, insertions, deletions | Effects of two mutations are not additive, and the data to learn the interaction does not exist |
| Purified-enzyme measurements | Cell-based IC50, patient outcomes | Cell assays add permeability, efflux and ATP concentration; those belong to a later phase |
| Ki, Kd, and IC50 with stated ATP conditions | IC50 with no stated ATP concentration | Without the ATP concentration an IC50 cannot be converted to Ki and is not comparable to anything |

Covalent drugs are excluded on the strength of a specific warning already
recorded in the project literature notes: Hoyt et al., J Med Chem 2024;67:2-16
(doi:10.1021/acs.jmedchem.3c01502) argues that covalent-inhibitor IC50 values
are assay-condition artefacts unless substrate conditions are stated, and that
mutant selectivity ratios must account for the variants' different ATP
affinities. Both points apply directly here.

Osimertinib therefore leaves the training table. It can return in a later phase
as a separate model with kinact/KI as the target, which is the honest way to
include it.

## The target variable

The thing predicted is a ratio, not an absolute potency. Ratios measured within
one laboratory cancel most assay bias; absolute values do not.

    y = log10( Ki_mutant / Ki_reference )

So y = 0 means no change, y = 1 means a 10-fold loss of potency, y = -0.3 means
the drug became about twice as potent against the mutant. The equivalent free
energy change, for comparison with the simulation literature, is
ddG = 1.36 * y kcal/mol at 25 C.

### Converting IC50 to Ki

Most published values are IC50, not Ki. For an ATP-competitive inhibitor the
Cheng-Prusoff relation applies:

    Ki = IC50 / (1 + [ATP] / Km_ATP)

Three rules govern this conversion, and they are the heart of the data phase.

1. Convert only when both the assay ATP concentration and that variant's own
   Km are known. A missing value is a dropped row, never an assumed one.
2. Use the Km of the variant being measured, not the wild-type Km. Using the
   wrong one is exactly the error that makes published selectivity ratios
   disagree.
3. Record in the table that the value was converted, with the ATP concentration
   and Km used. A converted value and a directly measured Ki are not the same
   kind of evidence.

A known hazard, already documented in this project: the ATP Km is disputed. One
source reports about 12 uM for wild-type EGFR, another roughly 100 uM, and the
conversion scales with that number. Where sources disagree, the row carries both
values and the model is trained and tested under each. If the conclusion changes
between them, that is a finding to report, not a parameter to pick.

Both halves of a ratio are taken from a single paper and a single assay wherever
possible. A mutant Ki from one laboratory divided by a reference Ki from another
is a last resort, and is flagged so it can be excluded in a sensitivity check.

## Baselines

A model is only interesting relative to what it beats. Three baselines are fixed
now, before any model exists, so the comparison cannot be chosen after the fact.

| Baseline | What it does | Why it is here |
|---|---|---|
| B0 — Constant | Predicts the training-set mean of y for every input | The floor. A model that does not beat this has learned nothing |
| B1 — Mechanistic | The existing competition model in `R/`: predicts the shift from the ATP competition term `1 + [ATP]/Ka` with published constants | The physics baseline. Beating it is the project's actual claim |
| B2 — Nearest neighbour | Returns the observed shift of the most chemically similar inhibitor against the same mutation | Catches the case where the model is only memorising close analogues |

B1 is already known to be wrong in at least one important case: this repository
found that simple ATP competition reproduces the measured erlotinib IC50 for
L858R but under-predicts the potency loss in L858R/T790M by at least tenfold.
The claim this project tests is whether a model trained across many kinases is
wrong less often, and whether it is wrong in interpretable ways.

If the learned model beats B0 but not B1, that is a real and reportable result:
it would say the mechanism carries the signal and the data adds nothing. This
specification treats that outcome as success of the experiment, not failure of
the project.

## Metrics and the evaluation split

### Metrics

| Metric | Definition | Why it is reported |
|---|---|---|
| RMSE | Root mean squared error in log10 units | The headline number. 0.7 means typically wrong by about 5-fold |
| Spearman rho | Rank correlation between predicted and observed shift | Does it rank variants correctly even if the magnitude is off? |
| Accuracy at 10-fold | Classification accuracy, resistant defined as y > 1, i.e. ddG > 1.36 kcal/mol | The threshold used in the Abl resistance literature, so the number is directly comparable |
| MCC | Matthews correlation coefficient on that same call | Accuracy is misleading when resistant cases are rare, and they are |
| Interval coverage | Fraction of truths falling inside the 90% interval | Checks that the uncertainty is honest. Should be close to 0.90 |
| Sign error rate | Fraction of cases where the direction is wrong | A drug predicted to improve that in fact fails is the costliest mistake |

### The split

Fixed here, and not to be changed after seeing results.

- Primary evaluation is leave-one-kinase-out. Train on all kinases except one,
  test on the held-out kinase, repeat for each. This measures the thing that
  matters: will it work on a kinase it has never seen?
- Secondary evaluation is leave-one-mutation-position-out within a kinase. This
  asks whether the model generalises to a new position or only interpolates
  between drugs at positions it already knows.
- A random row-wise split is reported too, and labelled as optimistic. It will
  look much better than the other two. Reporting all three makes the gap visible
  instead of hiding it.
- No tuning on the test fold. Hyperparameters are chosen by an inner split
  inside the training folds only.

One rule follows from the data's structure: all rows sharing a kinase go to the
same side of the split, always. Gatekeeper mutations in different kinases are
structurally analogous, so a model that saw Abl T315I can half-guess EGFR T790M.
That is a genuine kind of generalisation, but it must be measured across
kinases, not within one.

## Success criteria

| Level | Condition | What it means |
|---|---|---|
| Minimum | Beats B0 on leave-one-kinase-out RMSE, and interval coverage is between 0.85 and 0.95 | The model learned something and knows what it does not know |
| Target | RMSE <= 0.7 log10 units, accuracy at 10-fold >= 0.75, Spearman rho >= 0.5, and it beats B1 | Comparable to published single-kinase results, on more than one kinase |
| Stretch | Also beats B1 on the held-out EGFR L858R/T790M erlotinib case, within 3-fold | The project's own motivating case is predicted correctly by a model that never saw it |

### The condition that stops the project

If Phase 2 cannot assemble 200 rows across at least 4 kinases with verified
provenance, the learned-model approach is abandoned and the work becomes the
mechanistic analysis on the EGFR system alone. That decision is made on the row
count, not on a feeling, and it is made before any model is trained.

The roadmap's gate was 300 rows across 5 kinases. 200 across 4 is the hard floor
below which no honest leave-one-kinase-out evaluation is possible, since each
held-out fold needs enough rows to mean anything. Between 200 and 300 the
project continues with reduced ambition and says so in the write-up.

A negative result is written up with the same care as a positive one.

## Assumptions and threats to validity

| Assumption | If it is wrong | How it is checked |
|---|---|---|
| A potency ratio within one paper cancels assay bias | Ratios are not comparable across papers and the target is noise | Compare ratios for the same mutation-drug pair reported by two laboratories; the spread is the noise floor the model cannot beat |
| Cheng-Prusoff conversion is valid for these inhibitors | Converted Ki values are systematically wrong | Train once on directly measured Ki only, and compare |
| Resistance mutations are reported more often than neutral ones | The model learns the reporting bias, not the biology | Count how many rows have y near 0; if almost none, state it plainly as a limitation |
| Kinases share enough structure for cross-kinase transfer | Leave-one-kinase-out fails even though within-kinase works | This is exactly what the primary split measures |
| Published values are correct as printed | Garbage rows | Every row carries a verification flag, as `data/parameters.csv` already does |

The third row is the most serious and the least fixable. Nobody publishes a
paper reporting that a mutation changed nothing, so the training data will
over-represent large shifts and the model will over-predict resistance. This
cannot be corrected by cleverness, only disclosed, and the calibration check
will show it.

The fifth row is where the existing work gives this project an advantage. The
project literature notes document that in a summary of 16 citations, five
carried a wrong year, a wrong method, or a claim the paper does not make, and
two numbers were spliced from different papers. A project that assumes published
numbers are reliable would not survive contact with that reality. This one
assumes the opposite from the start.

## What Phase 2 starts with

The table schema is the contract between the data phase and everything after it.
It extends `data/parameters.csv` rather than replacing it.

| Column | Example | Note |
|---|---|---|
| `kinase_uniprot` | P00533 | |
| `kinase_name` | EGFR | |
| `mutation` | T790M | clinical numbering |
| `reference_form` | L858R | wild type, or the named parent mutant |
| `inhibitor_name` | erlotinib | |
| `inhibitor_smiles` | canonical SMILES | |
| `value_mutant`, `value_reference` | 41, 4.7 | |
| `value_units` | nM | |
| `value_type` | Ki | Ki, Kd or IC50 |
| `atp_conc_uM`, `km_atp_uM` | | needed for any conversion |
| `converted` | yes / no | whether Cheng-Prusoff was applied |
| `y_log10_ratio` | 0.94 | the target |
| `same_paper` | yes | both values from one study and one assay |
| `assay_conditions` | free text | as in the current CSV |
| `citation`, `doi` | | |
| `verified` | yes / no | checked against the source document |

The first rows already exist. `data/parameters.csv` holds the EGFR erlotinib
values against L858R and L858R/T790M, verified against the supplementary
information.

Phase 2 proceeds in this order: Platinum first, because it is already curated
for mutation effects; then Papyrus, which keeps variants separate by design;
then targeted manual curation of a few well-studied kinases — ABL1, EGFR, ALK,
KIT, BTK — where the resistance literature is dense.

## The Phase 1 gate

A reader can state what a correct prediction looks like, without asking the
author anything.

If any definition above is unclear, that is a defect in this document and it is
fixed before Phase 2 begins.
