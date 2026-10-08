# Methods

Sections marked **[done]** describe work completed. Sections marked
**[planned]** describe the intended approach and will be revised once the work
is carried out.

---

## 1. Software

| Tool | Version | Use |
|---|---|---|
| UCSF ChimeraX | fill in from Help > About | Structure visualisation, contact analysis, superposition |
| R | 4.5.2 | Numerical modelling |
| deSolve | 1.42 | Solving ordinary differential equations |
| RCSB PDB | accessed 2026 | Source of all structural data |

Structures are fetched from the PDB at run time rather than bundled with the
repository, so the analysis always uses the current deposited coordinates.

---

## 2. Structure selection **[done, criteria to be formalised]**

The initial structure was PDB 1M17: the EGFR tyrosine kinase domain
co-crystallised with erlotinib, chain A, covering EGFR_HUMAN residues 671-998.
Erlotinib is present as chemical component AQ4, stored as residue 999.

Selection criteria to be applied and documented in the next revision:

- Experimental method and resolution
- Ligand present, fully occupied, and well resolved in the electron density
- Binding-site residues complete, with no unmodelled gaps in the pocket
- Alternate conformations noted where present

1M17 contains 12 unmodelled residues and 14 atoms with alternate locations.
Neither falls within the binding site, but both are recorded here because they
constrain what the structure can support.

---

## 3. Residue numbering **[done]**

1M17 numbers residues from the mature protein. The clinical literature numbers
from the precursor, which includes a 24-residue signal peptide. The offset is
+24.

| Clinical | This structure | Role |
|---|---|---|
| T790 | THR 766 | gatekeeper |
| M793 | MET 769 | hinge |
| C797 | CYS 773 | covalent target of osimertinib |

The offset was established empirically rather than assumed. Atom names were
inspected with `info atoms`: residue 790 in this file carries CB, CG, CD1 and
CD2, which identifies a leucine; residue 766 carries OG1 and CG2, which
identifies a threonine. The binding-site residue list was then confirmed with
`info residues` on a 5 A zone around the ligand.

This check is recorded because applying clinical numbering directly to this
file would have produced a confident and entirely wrong analysis.

---

## 4. Binding site definition **[done]**

The binding site was defined as all residues with at least one atom within
5 A of any ligand atom, selected in ChimeraX with a zone selection on the
ligand.

Residues identified: LYS 721, GLU 738, MET 742, LEU 764, ILE 765, THR 766,
GLN 767, LEU 768, MET 769, PRO 770, PHE 771, GLY 772, CYS 773, LEU 820,
THR 830, ASP 831.

---

## 5. Contact analysis **[done]**

Contacts between the gatekeeper residue and erlotinib were identified with the
ChimeraX `contacts` command at an overlap cutoff of -2 A. This counts atom
pairs whose van der Waals surfaces come within 2 A of touching, so it captures
close approaches as well as direct contact.

Result: 20 contacts between THR766 and erlotinib.

Two constraints on interpretation:

- A contact count is descriptive geometry, not a thermodynamic quantity. It
  does not measure binding affinity or interaction energy.
- The count depends on the chosen cutoff. The value used is stated so the
  result can be reproduced or recomputed at a different threshold.

Before the measurement was accepted, the Models panel was checked to confirm a
single structure was loaded. An earlier run with the structure inadvertently
opened twice returned 80 contacts — four times the correct value, with no error
message. Duplicate models are a silent failure mode and are checked for
routinely.

---

## 6. Structural superposition **[planned]**

Wild-type and T790M structures will be superimposed with ChimeraX
`matchmaker`, which aligns by sequence before fitting coordinates. Root mean
square deviation will be reported over the aligned residues, together with the
number of residues used in the fit, since RMSD is not interpretable without it.

Residue numbering in the mutant structure will be verified independently. It
should not be assumed to share the offset used in 1M17.

---

## 7. Parameter collection **[in progress]**

All kinetic and affinity parameters will be taken from published experimental
measurements. None will be computed from structure.

For each parameter the following will be recorded: value, units, assay method,
experimental conditions, and full citation.

Assay conditions are recorded because affinity values measured by different
methods are not directly comparable. Where wild-type and mutant values are
compared, measurements from the same study and assay are preferred, and any
comparison across studies is flagged as such.

Values are held in `data/parameters.csv`, one row per parameter, with columns
for the construct measured, the assay method and conditions, the citation and
its DOI. A `verified` column records whether the value has been read off the
source paper itself rather than taken from a secondary note; model 4 prints a
warning listing every unverified row it uses. Rows still empty are open tasks,
and the scripts that need them fail with a message naming the missing symbol
instead of substituting a default.

---

## 8. Model formulation **[done, models 1-4]**

A compartmental model of receptor occupancy, expressed as ordinary
differential equations and solved numerically with `deSolve::ode` (lsoda).

State variables:

| Symbol | Meaning |
|---|---|
| D | free inhibitor concentration, a fixed input |
| A | free ATP concentration, a fixed input |
| R | free receptor |
| DR | inhibitor-receptor complex |
| AR | ATP-receptor complex |
| S | downstream signal |

Full equations, as implemented in model 3:

```
dR/dt  = ksyn - kdeg*R - konD*D*R + koffD*DR - konA*A*R + koffA*AR
dDR/dt = konD*D*R - koffD*DR - kint*DR
dAR/dt = konA*A*R - koffA*AR - kdeg*AR
dS/dt  = ktrans*R - kout*S
```

Fractional occupancy is `DR / (R + DR + AR)`.

Binding is mass action in both directions. Receptor synthesis is zero order
and degradation first order, applied to the ATP complex as well, since ATP
binding does not protect the receptor. The downstream signal is driven by free
receptor, because the inhibitor-bound receptor is inactive.

### Departure from the original plan: ATP competition

The plan put ATP competition out of scope. That was wrong for this question.
Erlotinib is ATP-competitive, and Yun et al. (2008) attribute T790M resistance
principally to increased ATP affinity, so a binding model without ATP cannot
represent the mechanism it is meant to explain, and any dose derived from it
would be wrong. ATP is therefore present from model 1 rather than added later.

The ATP term enters as `A/Ka` in the denominator of the occupancy expression,
and the dose required for occupancy f scales with `1 + A/Ka`.

### The four models

Each is a separate script that runs on its own.

| Script | Adds | Result |
|---|---|---|
| `R/01_binding_atp_competition.R` | reversible binding, ATP competition | `occupancy = (D/Kd)/(1 + D/Kd + A/Ka)`; dose for f is `f(1 + A/Ka)/(1 - f)` |
| `R/02_receptor_turnover.R` | synthesis, degradation, complex loss | same form with `Kd_app = (koffD + kint)/konD`; turnover shifts the curve along the dose axis without changing its shape |
| `R/03_downstream_signal.R` | lumped downstream signal | at occupancy f the retained signal is `1/(1 + (kint/kdeg)·f/(1 - f))`, independent of A |
| `R/04_mutant_dose.R` | measured values from `data/parameters.csv` | dose for 90% occupancy in L858R and L858R/T790M, and a comparison of predicted against measured IC50 |

Models 1 to 3 are dimensionless, so each runs before any parameter value is
available. Only model 4 reads measured values, and it reads them from the
parameter table rather than carrying them in code.

### Verification

Each script computes its steady state twice, once by integrating the ODEs and
once from the closed-form solution derived by hand, and stops with an error if
the two differ by more than a stated tolerance. Model 3 additionally checks the
full model against the one-line expression for retained signal. A wrong edit to
the equations therefore fails loudly rather than producing a plausible figure.

Reported differences on the author's machine (R 4.5.2, deSolve 1.42): 1.11e-16
for model 1, 1.11e-16 for model 2, and 3.13e-13 and 1.11e-16 for model 3's two
checks. The same equations were implemented independently in Python with SciPy
LSODA over randomised parameter sets, agreeing to 4.2e-15.

### Model 4 and what it does not yet settle

Model 4 compares L858R with L858R/T790M rather than wild-type with T790M,
because that is what the source measurements cover and because T790M arises
clinically on an already-mutant receptor.

Two approximations are carried and should be read as limitations. The ATP
Michaelis constant is used in place of an ATP dissociation constant, which
holds only when catalysis is slow relative to dissociation. And the inhibition
constant and the IC50 values come from different assay formats in the same
paper.

The script also tests the model against measurement, which is the only such
test in the project. Predicted IC50 is `Ki(1 + A/Km)`. For L858R this lands
close to the measured value; for L858R/T790M it is roughly tenfold too low. So
simple competitive binding with these two constants accounts for the sensitive
mutant and not for the resistant one. That gap is reported as a result.

Whether a calculated dose is achievable is not addressed, because the free
plasma concentration of erlotinib is not yet in the parameter table. The script
prints what is missing instead of estimating it.

---

## 9. Sensitivity analysis **[planned]**

Parameters will be varied individually across a plausible range and the effect
on predicted occupancy and dose recorded. Parameters with the weakest
literature support will be varied most widely.

---

## 10. Reproducibility

All analysis is scripted. ChimeraX commands are stored in `egfr_analysis.cxc`
and R code in numbered scripts. Software versions are recorded above.
Structures are fetched from the PDB rather than stored locally.

The intended standard is that a reader with the same software can run the
scripts and obtain the same figures and numbers without further instruction.

---

## 11. Limitations of the methods

- A crystal structure is a static, averaged snapshot under crystallisation
  conditions. It carries no information about conformational dynamics or
  binding kinetics.
- Contact counts are geometric descriptors and cannot be converted to binding
  energies.
- The mutant is inferred sterically from residue size rather than simulated.
- T790M also increases the receptor's affinity for ATP, and this contributes to
  resistance independently of any steric effect. The structural analysis does
  not address it, and the first version of the model does not include ATP
  competition.
- No experimental validation is performed. Every parameter is taken from the
  literature, and the model's predictions are not tested against data generated
  for this purpose. The one internal test available is the comparison of
  predicted against measured IC50 in model 4, and it fails for the double
  mutant by about tenfold.
- Free inhibitor and free ATP are treated as fixed inputs. Binding does not
  deplete either, and no absorption or distribution is modelled.
- The ATP Michaelis constant is used where the model calls for an ATP
  dissociation constant.
- Occupancy is not pathway inhibition. Model 3 shows the two differ by a factor
  that depends on `kint/kdeg`, a ratio this project has not obtained a value
  for, so no statement about signal loss in a cell is supported yet.
