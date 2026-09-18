# Nextstrain phylogenetic workflow

## Starting point

This workflow starts with segment-specific consensus nucleotide FASTA files and
sample metadata, after the study's read processing, lineage validation and sample
selection. It runs HA and NA independently. ConSurf is a separately documented
analysis and is not called by this workflow.

## Automated stages

| Stage | Command | Main output |
| --- | --- | --- |
| Input validation | `prepare_phylogeny_inputs.py` | Prepared FASTA, analysis/observed metadata and QC summary |
| Sequence indexing | `augur index` | `sequence_index.tsv` |
| Filtering | `augur filter` | `filtered.fasta`, `metadata.tsv`, retained identifiers and filter log |
| Alignment | `augur align` with MAFFT | `aligned.fasta` |
| Tree inference | `augur tree` with IQ-TREE | `tree_raw.nwk` |
| Temporal refinement | `augur refine` with TreeTime | `tree.nwk`, `branch_lengths.json` |
| Nucleotide reconstruction | `augur ancestral` | `nt_muts.json` |
| Amino-acid translation | `augur translate` | `aa_muts.json`, `aa_HA.fasta` or `aa_NA.fasta` |
| Geographic reconstruction | `augur traits` | `traits.json` |
| Auspice export | `augur export v2` | `auspice/flu_b_HA.json` or `auspice/flu_b_NA.json` |
| Dataset validation | `augur validate export-v2` | Validation log |
| iTOL preparation | `export_itol.py` | Divergence/time Newick and annotation files |
| Optional HA distance analysis | `patristic_distance.py` | `patristic_distances.csv` |

The precise commands are generated from `config/analysis.json` and recorded in
`work/<segment>/commands.sh`. Each stage must complete successfully before the next
starts. Outputs are not reused silently; choose a new `output_root` or archive
existing outputs with `--force` before rebuilding.

## Running

```bash
conda activate flu-b-phylogeny
python scripts/run_nextstrain.py --segment HA --dry-run
bash scripts/run_phylogenetic_pipeline.sh config/analysis.json 4
```

A dry run prints the plan even when local study files are absent. It is not a
successful analysis. The actual run requires confirmed parameters, input files
and the configured versions of Augur, MAFFT and IQ-TREE.

## Output interpretation

`tree_raw.nwk` is the IQ-TREE maximum-likelihood output. The later `tree.nwk` and
Auspice dataset reflect TreeTime refinement and the configured branch-length mode.
Do not describe the refined and unrefined branch lengths as interchangeable.
Divergence and temporal outputs have different units: the driver requests
per-site divergence, while the temporal export records branch durations in years.

The translated FASTA and mutation JSONs can contain internal ancestral records.
`aa_muts.json` lists changes along branches; a terminal branch list is not the
complete reference-relative genotype of that virus. For HA-NA matching, restrict
to terminal records and use verified isolate/segment identifiers. For reference-
relative mutation burden, compare aligned terminal sequences to the reference and
handle missing calls, successive substitutions and reversions explicitly. Do not
union root-to-tip mutation strings and count every historical event as a distinct
present-day substitution.

Clades are retained from supplied metadata. The workflow does not assign influenza
clades from an unspecified definition, estimate vaccine effectiveness, infer direct
transmission pairs or run individual clinical association models automatically.
The included aggregate scripts reproduce their named summaries only.

## Geographic annotations

The default trait columns are `country` and `region`. iTOL country/clade strips use
observed tip metadata when supplied, rather than substituting inferred ancestral
states. Study-specific definitions of introductions, exports and clusters must be
reported with their own filtering rules; trait reconstruction alone is not proof
of direct transmission or individual travel.

## Visualisation

Explore the exported datasets in Auspice. Then upload the corresponding Newick
and annotation files to iTOL. This follows the study order: phylogeny inference
first, iTOL visualisation afterward. No iTOL credentials or upload API keys are
required by the scripts.

## Documentation

- [Augur 27.2.0](https://docs.nextstrain.org/projects/augur/en/27.2.0/)
- [Alignment](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/align.html)
- [Tree inference](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/tree.html)
- [Refinement](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/refine.html)
- [Ancestral sequences](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/ancestral.html)
- [Translation](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/translate.html)
- [Traits](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/traits.html)
- [Export](https://docs.nextstrain.org/projects/augur/en/27.2.0/usage/cli/export.html)
