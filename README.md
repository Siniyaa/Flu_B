# Patristic_Distance
Patristic_Distance analysis
# patristic_distance

Compute **patristic (tree-path) distance** from every leaf in an
[Auspice v2](https://docs.nextstrain.org/projects/auspice/) phylogeny
(Nextstrain / Augur `*.json`) to a chosen reference tip, and bin the results
into divergence categories.

The tool was written for influenza B/Victoria genomic surveillance — measuring
how far each circulating HA sequence sits from a vaccine reference strain on the
tree — but it works for any Auspice v2 JSON and any reference tip.

## What it computes

Patristic distance between two tips *A* and *B* is the sum of branch lengths
along the unique path connecting them on the tree. Using the cumulative
root-to-tip divergence that Auspice stores on every node
(`node_attrs.div`, in substitutions per site), this reduces to:

```
d(A, B) = div(A) + div(B) − 2 × div(LCA(A, B))
```

where **LCA** is the lowest common ancestor of *A* and *B*. The `− 2 × div(LCA)`
term removes the shared `root → LCA` path that would otherwise be
double-counted.

This is the standard tree-path distance used across phylodynamics pipelines
(Nextstrain, TreeTime, BEAST, ete3, dendropy).

### What it is *not*

- **Not** raw Hamming / alignment distance between sequences.
- **Not** independent of the input tree — it reads the branch lengths already
  estimated by your Augur/TreeTime run, so the distances are only as good as
  that tree.

## Requirements

- Python 3.8+
- No third-party dependencies (standard library only)

## Installation

```bash
git clone https://github.com/<your-username>/patristic_distance.git
cd patristic_distance
```

## Usage

```bash
python patristic_distance.py <auspice.json> --reference <TIP_ID> [options]
```

### Example

```bash
python patristic_distance.py flu.json \
    --reference EPI_ISL_983345 \
    --group-countries "UAE" "United Arab Emirates" \
    --out uae_patristic_distances.csv
```

Drop `--group-countries` to score **every** leaf in the tree:

```bash
python patristic_distance.py flu.json --reference EPI_ISL_983345
```

### Arguments

| Argument | Required | Description |
|---|---|---|
| `json_file` | yes | Auspice v2 phylogeny JSON (e.g. `flu.json`). |
| `--reference` | yes | Reference tip name/id to measure distance to (e.g. `EPI_ISL_983345`). |
| `--group-countries` | no | One or more country values; only leaves matching these are reported. Omit to score all leaves. |
| `--out` | no | Output CSV path. Default: `patristic_distances.csv`. |

> **Note on country matching:** Auspice exports sometimes carry the same
> country under two different strings (e.g. `UAE` and `United Arab Emirates`).
> Pass every variant you need to `--group-countries`, or harmonize the metadata
> at the source before export.

## Output

A CSV sorted by ascending distance, one row per leaf:

| Column | Description |
|---|---|
| `name` | Leaf name / accession. |
| `country` | Country value from `node_attrs`. |
| `date` | Collection date from `node_attrs`. |
| `div` | Root-to-tip divergence of the leaf (subs/site). |
| `patristic_distance` | Tree-path distance to the reference (subs/site). |
| `category` | Divergence band (see below). |

The script also prints a summary (n, min, median, mean, max) and the
category distribution to stdout.

## Divergence categories

Distances are binned into five bands (substitutions/site):

| Category | Range |
|---|---|
| Very close | `< 0.005` |
| Close | `0.005 – 0.007` |
| Moderate | `0.007 – 0.010` |
| Distant | `0.010 – 0.020` |
| Very distant | `≥ 0.020` |

**These thresholds are dataset-specific, data-driven cutoffs — not a published
or clinical standard.** They were chosen to separate the observed modes in a
particular B/Victoria dataset. If you reuse this tool on other data, review the
distance distribution (e.g. as a histogram) and adjust the `CATEGORIES` list in
`patristic_distance.py` to fit your own data, and document your choice.

## Reproducibility

Because the script reads the `div` branch lengths from a specific Auspice JSON,
reproducing the numbers requires the **same** JSON. To make a repo
self-contained, either commit the exact `*.json` used or point to the precise
Augur build (config + input sequences) that produced it.

## License



## Citation

If you use this tool in a publication, please cite the associated manuscript

