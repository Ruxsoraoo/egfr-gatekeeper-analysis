# Parameter table

`parameters.csv` holds every numerical input to the models. It is deliberately
shipped with the `value` column empty: the project's rule for this phase is
that no parameter is assumed, computed from structure, or carried over from
memory. Each row is filled only when a published measurement has been read and
cited.

## Columns

| Column | Why it is recorded |
|---|---|
| `parameter` | What the quantity is, in words |
| `symbol` | The name used in the R scripts |
| `construct` | What was actually measured: full receptor, kinase domain, cell line |
| `value`, `units` | The measurement as published, not rescaled |
| `assay_method` | Affinities from different methods are not directly comparable |
| `assay_conditions` | ATP concentration in particular, since it shifts apparent inhibitor potency |
| `citation`, `pmid` | Full reference, so every number can be traced |
| `notes` | Anything that limits how the value may be used |

## Rules

- Where a wild-type and a mutant value are compared, values from the same study
  and the same assay are preferred. A comparison across studies is flagged as
  such in `notes`.
- An IC50 is not a Kd. If an IC50 is the only value available, the ATP
  concentration of that assay is recorded, because the conversion depends on it.
- A row with an empty `value` is an open task, not a missing file. The scripts
  that need a value fail with a clear message rather than substituting a
  default.
