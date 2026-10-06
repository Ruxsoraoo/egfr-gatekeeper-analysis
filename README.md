# EGFR gatekeeper residue and erlotinib binding

A structural analysis of why the T790M mutation causes resistance to erlotinib
in EGFR-mutant non-small cell lung cancer.

This is a self-directed learning project. The findings reproduce established
results in the literature; the purpose is to work through the structural
argument from the primary data rather than to report anything new.

## Question

Erlotinib is a reversible EGFR tyrosine kinase inhibitor. Patients respond,
then relapse, and the T790M gatekeeper mutation accounts for a substantial
share of acquired resistance. Where does that residue sit relative to the drug,
and how close is the contact?

## Data

| | |
|---|---|
| Structure | PDB 1M17 |
| Description | EGFR tyrosine kinase domain with erlotinib |
| Chain | A, EGFR_HUMAN residues 671-998 |
| Ligand | AQ4 (erlotinib), residue 999 |
| Software | UCSF ChimeraX |

## A note on residue numbering

1M17 numbers residues from the mature protein. Clinical numbering includes the
24-residue signal peptide, giving an offset of +24:

| Clinical | This structure | Role |
|---|---|---|
| T790 | THR 766 | gatekeeper |
| M793 | MET 769 | hinge |
| C797 | CYS 773 | covalent target of osimertinib |

This was verified by inspecting atom names rather than assumed. Residue 790 in
this file carries CB/CG/CD1/CD2 (leucine); residue 766 carries OG1/CG2
(threonine). Using the clinical numbers directly would have produced a
confident but wrong analysis.

## Result

**20 contacts** between THR766 and erlotinib at an overlap cutoff of -2 A.

The gatekeeper threonine is in direct contact with the bound inhibitor. This proximity motivated the original steric hypothesis for T790M resistance. That hypothesis has since been shown to be inadequate — Yun et al. (2008) demonstrated that the mutation increases the receptor's affinity for ATP by more than an order of magnitude, outcompeting ATP-competitive inhibitors. The structure alone cannot distinguish between the two mechanisms.

## Reproducing it

Install [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/), then:

```
open egfr_analysis.cxc
```

The script fetches the structure from the PDB, performs the analysis, and
writes the figure. Every step is commented.

## Files

| File | Contents |
|---|---|
| `egfr_analysis.cxc` | Annotated ChimeraX script |
| `egfr_gatekeeper.png` | Binding site figure |
| `egfr_session.cxs` | ChimeraX session (not tracked; written by the script when it runs) |
| `notes.txt` | Working notes |

## Limitations

- A single crystal structure is a static snapshot. It says nothing about
  conformational dynamics or binding kinetics.
- Contact counts depend on the chosen cutoff and are descriptive, not
  thermodynamic. They do not measure binding affinity.
- The T790M mutant itself is not modelled here. The inference about steric
  encroachment is based on the wild-type structure and the known size
  difference between threonine and methionine.
- The ATP-affinity mechanism is the principal explanation for T790M resistance and cannot be addressed by structural analysis alone. Testing it would require modelling ATP competition explicitly, which is the next stage of this project.

## Next steps

- Compare against a structure of the T790M mutant
- Add an osimertinib-bound structure and superimpose
- Build a receptor occupancy model in R using published binding affinities, to
  ask what dose would be required to restore target engagement in the mutant

## References

Yun CH, Mengwasser KE, Toms AV, et al. The T790M mutation in EGFR kinase causes drug resistance by increasing the affinity for ATP. Proc Natl Acad Sci USA. 2008;105(6):2070-2075.

## Author

Ruxsorabonu Obidova — Assistant Lecturer in Pharmacology,
Central Asian Medical University, Fergana, Uzbekistan
