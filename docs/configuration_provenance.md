# Workflow configuration and provenance

## Scope

The workflow implements the sequence of analysis stages described for the study:
MAFFT alignment, IQ-TREE inference, TreeTime refinement, ancestral sequence and
trait reconstruction, Auspice visualisation and subsequent iTOL annotation.
An Auspice JSON and mutation outputs do not contain every original command-line
option. The configuration therefore distinguishes recorded settings from values
that require checking before a new build.

## Settings supported by the available records

| Item | Recorded information | Use in this repository |
| --- | --- | --- |
| Augur | `generated_by` in the HA/NA outputs records version 27.2.0 | Pinned in the workflow environment |
| MAFFT | NA alignment log records version 7.526 | Pinned; confirm the HA run separately |
| IQ-TREE | NA tree log records version 2.3.6 | Pinned; confirm the HA run separately |
| NA tree model | Original NA command contains `-m GTR` | NA configuration uses `GTR` |
| NA IQ-TREE seed | Original NA log records 666962 | Used only for NA tree inference, not asserted as a TreeTime seed |
| HA reference tip | Study HA tree and manuscript identify `EPI_ISL_983345` | Configured for HA reference/patristic lookup |
| NA annotation record | Supplied NA GenBank record is `B_Aus_1359_NA`, with a CDS named `NA` | Reference name must be verified against the actual input FASTA and metadata |

The NA IQ-TREE command includes `--ninit 2 -n 2 --epsilon 0.05 -T AUTO --redo`,
matching the corresponding Augur IQ-TREE defaults. Thread count is a run argument.
The original NA reference is not retained as a tip under this name in the supplied
final NA Auspice dataset; its retention and rooting history must therefore be
checked rather than inferred from the GenBank identifier alone.

## Settings to confirm

`parameters_confirmed` is initially `false` for both segments. Confirm the actual
FASTA, metadata and reference GenBank files and set the following as applicable:

- `substitution_model`: the exact IQ-TREE model string for each segment. The
  original HA command was not supplied. Do not assume a rate-heterogeneity suffix
  from general manuscript wording; the available NA log specifies plain `GTR`.
- `min_length`, `max_length`, `max_ambiguous_fraction`: actual QC criteria.
- `min_date`, `max_date`, `exclude_ambiguous_dates_by`: date-selection criteria.
- `clock_filter_iqd`: the actual clock-outlier threshold. The manuscript names the
  option without giving its numeric value.
- `root`: an explicit rooting choice, or `null` to use `reference_tip` in this
  implementation. Reference-based rooting must be appropriate for the build.
- `tree_seed` and `seed`: separate IQ-TREE and TreeTime/ancestral seeds.
- `branch_length_inference`, `date_inference`, `coalescent`, `keep_polytomies`,
  `date_confidence`, `keep_ambiguous` and `keep_overhangs`: inference choices.

For optional threshold/date/seed settings, `null` means the corresponding option
is not passed, not that the historical study used a null value. The supplied
inference choices are explicit implementation defaults, not recovered historical
settings. A new deliberate configuration should be identified as a reanalysis.

## References and coordinate systems

The reference must be present as a unique terminal identifier in the input FASTA
and metadata. Its sequence must exactly match the GenBank sequence. That GenBank
record must have the correct CDS annotation and `/gene="HA"` or `/gene="NA"`.
The original HA GenBank file was not supplied with the existing output archive;
use the original annotation rather than deriving coordinates from mutation labels.
Do not create a collection date for a reference from its name or publication year.

`augur align --reference-name` strips alignment insertions relative to that
reference. `augur ancestral --root-sequence` defines the baseline for mutation
calling; it does not by itself force the inferred ancestral sequence to equal the
vaccine reference. These are distinct from the tree-rooting option.

## Input selection

Upstream cohort and background selection must be recorded separately. This workflow
does not sample a fixed number of non-UAE viruses or recreate an undocumented
comparator. It rejects duplicate identifiers, but does not collapse distinct
viruses merely because their sequences are identical. Optional inclusion lists
bypass Augur filters; the reference is included explicitly after passing input QC.

The automatic output metadata omit age, clinical category and facility fields.
The analysis copy marks blank geographic traits as `?` for reconstruction; the
separate observed metadata remain unchanged. Do not interpret inferred country
values as observed travel histories. Clinical analyses require their own controlled
input tables and analysis code.

## Environment and execution record

The environment specification preserves known versions but is not an original
lockfile. Record the resolved environment, configuration and `run_manifest.json`
for each execution. The manifest captures input checksums, executable versions,
commands, stage completion and errors. A successful new run documents that run;
it does not retrospectively verify previously reported results.
