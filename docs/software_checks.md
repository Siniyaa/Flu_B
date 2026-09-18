# Software checks

## Automated checks

```bash
python -m compileall -q scripts tests
python -m unittest discover -s tests -v
bash -n scripts/run_phylogenetic_pipeline.sh scripts/run_all.sh
```

The 50 Python helper tests passed in the local validation environment. The GitHub
workflow reruns them under Python 3.11. They cover:

- Patristic calculation, reference self-distance, terminal-only output, country
  filters and category boundaries.
- FASTA/metadata matching, duplicate identifiers, reference sequence/CDS checks,
  ambiguity filtering and exclusion of clinical columns from browser exports.
- Newick path lengths, terminal annotation rows, deterministic colour strips,
  authoritative observed metadata and invalid branch-scale handling.
- Command construction, reference selection, parameter-confirmation guards,
  IQ-TREE versus TreeTime seeds, geographic trait options and separate HA/NA builds.

The artificial trees and sequences are software fixtures, not inputs used to
produce the study's patristic-distance or HA-NA results. The original patristic
script is retained without analytical changes. The tests do not execute the
external Augur, MAFFT or IQ-TREE inference chain.

## Checks using existing study outputs

The iTOL divergence exporter was checked against the supplied HA and NA Auspice
files during package preparation. It retained the terminal identifiers and checked
reference-to-tip distances against the source tree. Study sequence records,
metadata and generated trees are not included in this repository update.

Some supplied point node dates decreased along parent-child edges. For those
inputs, the exporter omits temporal Newick output rather than changing dates or
writing negative-duration branches. Divergence export is independent of temporal
export. These checks concern conversion of point node dates to branch durations,
not a reanalysis of clock-model fit. Review the temporal inference and source export
before interpreting calendar-time branch lengths in iTOL.

## Execution scope

The complete external inference chain has not been run as part of this source
upload. The environment specification is not an original lockfile, and its
resolution must be checked on the execution platform. Helper tests do not establish
reproduction of fitted study trees, adjusted models, statistics or figures.

For a full execution, supply the actual inputs, complete the configuration,
create the environment and run the phylogenetic driver. Retain the resolved
environment, `commands.sh`, stage logs and `run_manifest.json`. iTOL upload and
figure styling remain separately recorded visualisation steps. The optional Firth
script and aggregate HA-NA analyses are separate from the phylogenetic workflow.
