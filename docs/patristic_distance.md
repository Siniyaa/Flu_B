# Patristic-distance analysis

## Purpose

`patristic_distance.py` calculates the tree-path distance from terminal sequences
to a chosen reference in an Auspice v2 JSON phylogeny. It does not infer a tree,
calculate an alignment distance, construct a country comparator, or estimate enrichment.

## Study data

The study analysis used a phylogeny inferred from actual influenza B/Victoria HA
sequences from the UAE and international comparison data. Distances were calculated
from the terminal sequences in this study tree to the selected reference sequence.
Synthetic data were not used to generate the study results.

The input path in the commands below refers to a locally supplied study tree.
The separate synthetic example distributed with the code is used only for automated
software tests; it is not the study dataset or a source of manuscript results.

## Inputs

The positional argument is the path to an Auspice JSON file with a top-level
`tree` object. Nodes use `name`, optional `children`, and `node_attrs` fields.

| Field | Use |
| --- | --- |
| `name` | Exact identifier; used for reference lookup and CSV output |
| `children` | Parent-child structure; nodes without children are leaves |
| `node_attrs.div` | Cumulative root-to-node divergence |
| `node_attrs.country.value` | Optional output field and country filter |
| `node_attrs.date.value` | Optional output date; copied without conversion |

Supply a unique reference **tip** name. The lookup searches all nodes, so choosing
an internal-node name would select an ancestor rather than a reference sequence.
The script reports an error when the reference is absent or has no divergence value.

Divergence values used in the calculation must be finite numbers on a common
cumulative scale, in substitutions per site for the supplied category thresholds.
They should not decrease from parent to child. The reference and the relevant
lowest common ancestors must have valid values. The script does not comprehensively
validate these conditions; check them in the upstream tree export.

The date field is read only from `node_attrs.date.value`. It is left blank when
absent; `num_date` is not converted to a calendar date.

## Calculation

For a terminal sequence `tip` and reference `ref`:

```text
d(tip, ref) = div(tip) + div(ref) - 2 * div(LCA(tip, ref))
```

`LCA` is their lowest common ancestor. Subtracting its cumulative divergence twice
removes the shared root-to-ancestor path. Internal nodes are used to locate this
ancestor but are not included as observations in the CSV.

For example, a reference divergence of 0.006, tip divergence of 0.007, and LCA
divergence of 0.004 give a patristic distance of 0.005, not 0.013.

The existing implementation clips negative computed distances to zero. This does
not validate non-monotonic or otherwise malformed divergence values, which must be
checked before interpreting results.

## Commands

Run from the repository root. Create the output directory first; the script does
not create parent directories.

```bash
mkdir -p results/patristic
python scripts/patristic_distance.py raw_data/ha_auspice.json \
  --reference EPI_ISL_983345 \
  --out results/patristic/all_tips.csv
```

Replace the input path and reference name with the actual tree and intended tip.
To keep only selected country labels:

```bash
python scripts/patristic_distance.py raw_data/ha_auspice.json \
  --reference EPI_ISL_983345 \
  --group-countries "UAE" "United Arab Emirates" \
  --out results/patristic/uae_tips.csv
```

The country filter controls output leaves, not the tree used to find ancestors.
It is exact and case-sensitive; no country labels are standardised automatically.
Omitting the filter scores all leaves with divergence values. A missing country
is retained in an unfiltered run but will not match a specified country label.

## Distance categories

| Category | Lower bound, inclusive | Upper bound, exclusive |
| --- | --- | --- |
| Very close | 0 | 0.005 |
| Close | 0.005 | 0.007 |
| Moderate | 0.007 | 0.010 |
| Distant | 0.010 | 0.020 |
| Very distant | 0.020 | No upper bound |

These thresholds are the supplied descriptive categories. They do not establish
vaccine mismatch, immune escape, antiviral resistance, or a transmission link.

## Outputs and denominators

The output CSV has these columns:

```text
name,country,date,div,patristic_distance,category
```

Rows are sorted by patristic distance. Both distance columns are rounded to six
decimal places for output. Category assignment uses the computed distance before
rounding; a rounded boundary value can therefore appear to fall in a neighbouring
category. The console summary uses the rounded distances written to the CSV.

The reported denominator is the number of selected leaves with a divergence
value. Leaves missing this value are omitted. The reference contributes a
zero-distance row if it is a leaf and passes the country filter. No further
reference exclusion, sampling, deduplication, or non-UAE selection is performed.

Document any upstream tree filtering and output filtering separately when reporting
study denominators. Country filtering alone does not produce a matched or
representative background sample.

## Software checks

```bash
python -m unittest discover -s tests -p 'test_patristic_distance.py' -v
```

These tests are separate from the study analysis and use a small synthetic tree
with known expected distances. They check shared-ancestor subtraction, reference
self-distance, terminal-only output, country filters, missing tip divergence,
category boundaries, expected CSV output, and missing-reference errors.
Test outputs are not study results. Passing these tests checks software behaviour;
it does not reproduce or validate the manuscript's numerical findings.
