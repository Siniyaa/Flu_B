# Influenza B/Victoria Genomic and Structural Analysis

Code and documentation for influenza B/Victoria analyses using UAE and international comparison sequences. Sequence processing and phylogenetic inference use **Nextstrain/Augur**. Trees are explored in **Auspice** and subsequently visualised and annotated in **iTOL**.

The workflow starts with consensus sequences and metadata and produces phylogenies, mutation annotations, Auspice datasets and iTOL-ready tree files. Separate components calculate HA patristic distances and summarise aggregate HA-NA co-occurrence data.

## Analysis components

| Component | Documentation |
| --- | --- |
| Nextstrain phylogenetic workflow | [Nextstrain workflow](docs/nextstrain_workflow.md) |
| Input settings and execution records | [Configuration and provenance](docs/configuration_provenance.md) |
| iTOL tree visualisation | [iTOL visualisation](docs/itol_visualisation.md) |
| HA patristic distances | [Patristic distances](docs/patristic_distance.md) |
| HA-NA substitution co-occurrence | [HA-NA analysis](docs/ha_na_analysis.md) |
| Conservation scores and structural visualisation | [Conservation-score README](README_Conservation_Score.md) |

The existing `README_Conservation_Score.md`, `pymol_code.py` and root-level `patristic_distance.py` are retained unchanged. The pipeline uses an identical copy of the patristic script in `scripts/`.

## Workflow

```text
UAE and global consensus sequences + metadata
    -> Input validation and quality filtering
    -> MAFFT alignment                         (augur align)
    -> IQ-TREE maximum-likelihood phylogeny    (augur tree)
    -> TreeTime temporal refinement           (augur refine)
    -> Ancestral sequences and translation    (augur ancestral / translate)
    -> Ancestral geographic traits            (augur traits)
    -> Auspice dataset                        (augur export v2)
    -> iTOL tree upload, annotation and visualisation
```

HA and NA are processed separately. HA-NA matching uses terminal virus identifiers, not shared internal-node labels in independently inferred trees.

## Repository contents

```text
config/                    Input paths, segment parameters and Auspice settings
scripts/                   Workflow drivers, export and analysis scripts
data/                      Aggregate clade, burden and association tables
results/                   Supplied aggregate statistical summaries
docs/                      Analysis instructions and interpretation
tests/                     Automated software checks
examples/                  Small artificial fixtures for software tests
.github/workflows/         Python continuous-integration checks
environment.yml            Phylogenetic workflow environment
requirements.txt           Python helper dependencies
README_Conservation_Score.md
pymol_code.py
patristic_distance.py
```

Rendered figures are generated locally with the supplied figure script rather than included in this source update.

## Installation

```bash
conda env create -f environment.yml
conda activate flu-b-phylogeny
```

The environment specifies Augur, MAFFT and IQ-TREE versions documented in the available analysis records. It is not an original environment lockfile. After installation, record the resolved environment:

```bash
conda env export --no-builds > environment.resolved.yml
```

For Python-only helper tests and aggregate summaries, use `pip install -r requirements.txt`. This does not install the external phylogenetic tools. Auspice viewing requires a configured Nextstrain CLI/runtime or a separate Auspice installation. iTOL is used through its web interface.

## Study inputs and configuration

Supply the actual study inputs locally:

```text
raw_data/
  HA/
    sequences.fasta
    metadata.tsv
    reference.gb
  NA/
    sequences.fasta
    metadata.tsv
    reference.gb
```

FASTA headers must match metadata `strain` identifiers. The reference must be present in the FASTA and metadata, and its sequence must match the reference GenBank record. Metadata require `strain` and `date`; the default trait analysis also requires `country` and `region`. Existing `clade` annotations are optional.

Review `config/analysis.json` before running. Confirm input paths, references, QC thresholds, the substitution model and clock settings against the analysis records. Set `parameters_confirmed` to `true` only after those settings have been checked. Unavailable historical settings are not silently treated as known. See [configuration and provenance](docs/configuration_provenance.md).

## Run the Nextstrain workflow

Print a command plan without running tools or reading study inputs:

```bash
python scripts/run_nextstrain.py --segment HA --dry-run
python scripts/run_nextstrain.py --segment NA --dry-run
```

Run both segment builds after configuration:

```bash
bash scripts/run_phylogenetic_pipeline.sh config/analysis.json 4
```

Run one segment:

```bash
python scripts/run_nextstrain.py --config config/analysis.json --segment HA --threads 4
```

The final argument to the launcher is the number of threads. Outputs go to `work/HA/` and `work/NA/`. Each build records input checksums, tool versions, commands, logs and a run manifest. Existing outputs are not overwritten; `--force` archives a segment output directory before rebuilding.

## View in Auspice and iTOL

With a configured Nextstrain runtime:

```bash
nextstrain view work/HA/auspice
nextstrain view work/NA/auspice
```

Upload `work/HA/itol/HA_divergence.nwk` or the corresponding NA file to iTOL, then add the generated country/clade annotation files. Temporal Newick files are written separately only when source dates pass validation.

Export an existing study Auspice tree without rerunning inference:

```bash
python scripts/export_itol.py \
  --auspice raw_data/ha_auspice.json \
  --metadata raw_data/HA/metadata.tsv \
  --outdir work/HA/itol \
  --prefix HA \
  --columns country clade
```

The exporter creates local files. Uploading to iTOL, arranging annotations and exporting the final figure are separate visualisation steps.

## HA patristic distances

The study analysis used actual influenza B/Victoria sequence data. Supply the study HA tree and its exact reference identifier:

```bash
mkdir -p results/patristic
python scripts/patristic_distance.py raw_data/ha_auspice.json \
  --reference EPI_ISL_983345 \
  --out results/patristic/all_tips.csv
```

The script calculates tree-path distances on an existing phylogeny; it does not infer a tree. Its categories are descriptive genetic-distance bins, not vaccine-effectiveness estimates. The artificial example in `examples/` is exclusively a software-test fixture and was not used to generate study results.

## HA-NA aggregate analyses

```bash
bash scripts/run_all.sh
```

This generates Figure 6 and aggregate burden summaries from `data/`. It does not reconstruct mutation profiles, calculate global burden correlations or perform the optional Firth regression. Terminal-pair checks and Firth regression require separate individual-record inputs. See [HA-NA analysis](docs/ha_na_analysis.md).

## Software checks

```bash
python -m compileall -q scripts tests
python -m unittest discover -s tests -v
bash -n scripts/run_phylogenetic_pipeline.sh scripts/run_all.sh
```

The tests check helpers and command construction; they do not run the external inference chain or establish reproduction of all study results. See [software checks](docs/software_checks.md).

## Data and reproducibility

The automated workflow begins with consensus FASTA files. Laboratory sequencing, read-level QC, assembly, lineage validation and study-specific background selection are upstream steps. A full study rebuild requires the original inputs and confirmed configuration.

Raw study sequences, clinical metadata, manuscripts and private correspondence are not included. Review identifiers, metadata and applicable data-use permissions before publishing generated trees or uploading them to an external service. Conservation-score analysis remains documented separately in [its existing README](README_Conservation_Score.md).

## Software documentation

- [Augur 27.2.0](https://docs.nextstrain.org/projects/augur/en/27.2.0/)
- [Nextstrain viewing](https://docs.nextstrain.org/projects/cli/en/stable/commands/view/)
- [iTOL help](https://itol.embl.de/help.cgi)
