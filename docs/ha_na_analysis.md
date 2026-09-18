# HA-NA substitution co-occurrence analysis

## Scope

The scripts in this component reproduce the aggregate clade and HA-burden
summaries and Figure 6 from the supplied count tables. A separate script checks
HA-NA matching at terminal-tip and all-node levels, and an optional R script
fits a Firth-penalised model using individual-level data.

The aggregate scripts do not reconstruct mutation states from sequence data,
recalculate the supplied top-pairwise-association table, or estimate global
burden correlations or geographic enrichment.

## Aggregate inputs

| File | Contents |
| --- | --- |
| `data/figure6_clade_counts.csv` | Clade counts in the matched UAE cohort |
| `data/figure6_high_burden_by_clade.csv` | High-burden and total counts by clade |
| `data/figure6_na_by_burden.csv` | Focal NA-substitution carriage by HA-burden class |
| `data/top_pairwise_associations_table.csv` | Supplied table of selected pairwise associations |

High HA-epitope burden is defined in these inputs as more than five substitutions.
The focal NA substitutions are P336T, G378E, I333R, D342N, and R295S. The scripts
use the aggregate categories as supplied; they do not derive the categories
from raw sequences.

## Reproduce the figure and aggregate summaries

Run both scripts:

```bash
bash scripts/run_all.sh
```

Or run them individually:

```bash
python scripts/create_figure6_clade_structure.py \
  --data-dir data \
  --outdir figures \
  --results-dir results \
  --basename Figure6_clade_structure \
  --formats png pdf svg \
  --dpi 600

python scripts/summarize_ha_na_clade_separation.py \
  --data-dir data \
  --outdir results
```

The summary script writes the proportion with high HA burden by clade and the
unadjusted high-versus-low HA-burden comparison for focal NA-substitution
carriage. It reports a two-sided Fisher exact P value and a Woolf confidence
interval for the odds ratio. It does not fit a clade-adjusted model.

## Terminal-pair checks

For local HA and NA JSON files:

```bash
python scripts/check_augur_terminal_pairs.py \
  --ha-aa-muts raw_data/ha_aa_muts.json \
  --ha-auspice raw_data/ha_auspice.json \
  --na-aa-muts raw_data/na_aa_muts.json \
  --na-auspice raw_data/na_auspice.json \
  --outdir results/unit_of_analysis
```

The script removes terminal `_HA` and `_NA` name suffixes before matching.
Terminal status is determined from the Auspice trees. The all-node comparison
is a record-name audit, not a matched-isolate dataset; shared internal-node
labels must not be treated as paired viruses. The script writes matched names
and count summaries, not mutation burdens or correlation estimates.

NA inputs can alternatively be supplied in a ZIP using `--na-zip`,
`--na-aa-muts-inner`, and `--na-auspice-inner`.

## Optional Firth-regression analysis

The R script requires a local CSV with the following columns:

```text
focal_NA_any,HA_epitope_burden,clade
```

`focal_NA_any` is a binary outcome. The script derives high-burden status using
the greater-than-five threshold and fits:

```text
focal_NA_any ~ high_HA_epitope_burden + clade
```

Install `logistf` in R, then run:

```r
install.packages("logistf")
```

```bash
mkdir -p results/firth
Rscript scripts/firth_logistic_sensitivity.R \
  raw_data/matched_ha_na.csv results/firth/model_summary.txt
```

This step is separate from `run_all.sh` and requires individual-level inputs.
A supplied analysis script is not a record that the model was run on the study
data; retain the actual inputs, output, software versions, and execution logs
when reporting an analysis.

## Interpretation

The unadjusted aggregate association is not a clade-adjusted estimate. A
sequence-defined dual-burden profile denotes high HA-epitope burden plus at
least one focal NA substitution; it describes co-carriage rather than a measured
functional phenotype. Any sequence-level extension must document its reference,
site definitions, sample matching, missing-data handling, and comparison groups.
