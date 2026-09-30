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

## 7. Parameter collection **[planned]**

All kinetic and affinity parameters will be taken from published experimental
measurements. None will be computed from structure.

For each parameter the following will be recorded: value, units, assay method,
experimental conditions, and full citation.

Assay conditions are recorded because affinity values measured by different
methods are not directly comparable. Where wild-type and mutant values are
compared, measurements from the same study and assay will be preferred, and
any comparison across studies will be flagged as such.

---

## 8. Model formulation **[planned]**

A compartmental model of receptor occupancy, expressed as ordinary
differential equations and solved numerically.

State variables:

| Symbol | Meaning |
|---|---|
| D | free drug concentration |
| R | free receptor |
| DR | drug-receptor complex |
| S | downstream signal |

Equations:

```
dR/dt  = ksyn - kdeg*R - kon*D*R + koff*DR
dDR/dt = kon*D*R - koff*DR - kint*DR
dS/dt  = ktrans*R - kout*S
```

Fractional occupancy is calculated as `DR / (R + DR)`.

Binding is treated as mass action. Receptor synthesis is zero order,
degradation first order. The downstream signal is driven by free receptor,
since the drug-bound receptor is inactive.

The system will be solved with `deSolve::ode` using the lsoda method.

The model will be built in four stages, each a separate script that runs
independently:

1. Reversible binding only
2. Binding with receptor turnover
3. Turnover with downstream signal
4. Mutant affinity values substituted

Building in stages makes each addition's effect visible and keeps each step
explicable before the next is added.

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
  for this purpose.
