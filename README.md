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
- The ATP-affinity mechanism cannot be addressed by structural analysis alone. Testing it requires modelling ATP competition explicitly, which is what the models below do.

## The models

Four R scripts in `R/`, each running on its own, take the argument from
binding to dose. ATP competition is present from the first one, because
erlotinib is ATP-competitive and the structure alone cannot settle the
mechanism.

| Script | Question it answers |
|---|---|
| `01_binding_atp_competition.R` | How does ATP competition change the dose needed for a given occupancy? |
| `02_receptor_turnover.R` | How does receptor turnover change steady-state occupancy? |
| `03_downstream_signal.R` | Does 90% occupancy produce a 90% loss of signal? |
| `04_mutant_dose.R` | What erlotinib concentration restores occupancy in the mutant, and does the model match measurement? |

Two results so far. Occupancy is not signal loss: at occupancy f the retained
signal is `1/(1 + (kint/kdeg)·f/(1 - f))`, which equals `1 - f` only when the
drug-bound complex is cleared at the same rate as free receptor is degraded.
And the competitive model reproduces the measured erlotinib IC50 for L858R but
under-predicts the loss of potency in L858R/T790M by at least tenfold, so
simple ATP competition with these constants does not account for the
resistance.

Every parameter lives in `data/parameters.csv` with its assay conditions, its
citation and a flag recording whether it has been checked against the source.

## Next steps

- Compare against a structure of the T790M mutant
- Add an osimertinib-bound structure and superimpose
- Add the ATP affinities from Yun et al. 2008 and repeat the dose calculation
  under both labs' values, as a sensitivity analysis
- Add the free plasma concentration of erlotinib, without which no dose can be
  called achievable

## References

Yun CH, Mengwasser KE, Toms AV, et al. The T790M mutation in EGFR kinase causes drug resistance by increasing the affinity for ATP. Proc Natl Acad Sci USA. 2008;105(6):2070-2075.

van Alderwerelt van Rosenburgh IK, Lu DM, Grant MJ, Stayrook SE, Phadke M, Walther Z, Goldberg SB, Politi K, Lemmon MA, Ashtekar KD, Tsutsui Y. Biochemical and structural basis for differential inhibitor sensitivity of EGFR with distinct exon 19 mutations. Nat Commun. 2022;13:6791. doi:10.1038/s41467-022-34398-z

The source of the inhibition constants and ATP affinities in
`data/parameters.csv`. Its measurements are made in one laboratory with one
assay, so the mutant-to-mutant ratios are internally consistent.

These two papers disagree about the size of the effect. Yun et al. attribute
T790M resistance principally to a large increase in the receptor's affinity for
ATP. In the 2022 measurements the ATP affinity of L858R/T790M is only about
1.6-fold tighter than L858R, while the erlotinib inhibition constant shifts
roughly ninefold, so in that dataset the affinity change carries the
resistance rather than ATP competition. The project uses the 2022 values
because they come from a single internally consistent set, and the planned
sensitivity analysis repeats the dose calculation under Yun's ATP values rather
than treating either as settled.

## Author

Ruxsorabonu Obidova — Assistant Lecturer in Pharmacology,
Central Asian Medical University, Fergana, Uzbekistan
