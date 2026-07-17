# conservation-scores-pymol
# Influenza Protein Conservation and Structural Visualization Pipeline

## Overview

This repository contains a reproducible workflow for analyzing influenza protein mutations using sequence analysis, multiple sequence alignment, conservation scoring, and structural visualization.

The workflow integrates:

- DNA to protein translation using Python
- Multiple sequence alignment using MAFFT
- Evolutionary conservation analysis using ConSurf
- Structural visualization using PyMOL
- Mutation and antigenic site mapping

This pipeline can be applied to influenza A proteins such as Hemagglutinin (HA) and Neuraminidase (NA).

---

# Workflow Summary

```text
DNA Sequence
    ↓
Protein Translation (Python)
    ↓
Protein FASTA File
    ↓
Multiple Sequence Alignment (MAFFT)
    ↓
Aligned FASTA File
    ↓
ConSurf Analysis
    ↓
Conservation Scores + PDB Files
    ↓
PyMOL Structural Visualization
```

---

# Step 1 — Convert DNA Sequence to Protein Sequence

## Input

- DNA sequence in FASTA format

Example:

```fasta
>sample_sequence
ATGGAGAAA...
```

---

## Python Script for DNA to Protein Translation

Save as:

```text
dna_to_protein.py
```

```python
from Bio.Seq import Seq
from Bio import SeqIO

input_fasta = "dna_sequence.fasta"
output_fasta = "protein_sequence.fasta"

with open(output_fasta, "w") as outfile:
    for record in SeqIO.parse(input_fasta, "fasta"):
        protein_seq = Seq(str(record.seq)).translate(to_stop=True)

        outfile.write(f">{record.id}\n")
        outfile.write(str(protein_seq) + "\n")

print("Protein translation completed.")
```

---

# Install Required Python Package

```bash
pip install biopython
```

---

# Run the Script

```bash
python dna_to_protein.py
```

---

# Output

```text
protein_sequence.fasta
```

This FASTA file contains translated protein sequences for downstream analysis.

---

# Step 2 — Multiple Sequence Alignment Using MAFFT

## Install MAFFT

Official website:

https://mafft.cbrc.jp/alignment/software/

---

# Input

```text
protein_sequence.fasta
```

---

# Run MAFFT Alignment

```bash
mafft protein_sequence.fasta > aligned_sequences.fasta
```

---

# Output

```text
aligned_sequences.fasta
```

This aligned FASTA file will be used for conservation analysis.

---

# Step 3 — ConSurf Conservation Analysis

## ConSurf Website

https://consurf.tau.ac.il/

---

# Upload Files

Upload:

- Protein structure (PDB)
- Aligned FASTA file

Example PDB structures:

| Protein | PDB ID |
|---|---|
| H3 Hemagglutinin | 8G5B |
| H1 Hemagglutinin | 1RU7 |
| N1 Neuraminidase | 3TI6 |

---

# ConSurf Processing

ConSurf performs:

- Homology search
- Multiple sequence analysis
- Phylogenetic analysis
- Conservation scoring

---

# Download ConSurf Results

After analysis, download:

- Conservation score table
- Colored PDB structure
- PyMOL session/script
- ConSurf grades file

Important files include:

```text
consurf_grades.txt
consurf_grades.pdb
```

---

# Step 4 — Structural Visualization Using PyMOL

## Load Structure

```python
load consurf_grades.pdb, HA
```

---

# Basic Visualization

```python
hide everything
show cartoon
color gray80
bg_color white
```

---

# Select Chain

```python
select HA_A, chain A
show cartoon, HA_A
```

---

# Visualize Mutation

Example mutation:

```python
select mut, HA_A and resi 302
show sticks, mut
color red, mut
zoom mut, 8
```

---

# Label Mutation

```python
label mut and name CA, "I302S"
set label_size, 18
set label_color, black
```

---

# Highlight Nearby Region

```python
select epi, byres (HA_A within 6 of mut)
show surface, epi
color yellow, epi
set transparency, 0.4
```

---

# Apply Conservation Coloring

```python
spectrum b, blue_white_red
```

Typical interpretation:

- Blue → variable residues
- White → intermediate conservation
- Red → highly conserved residues

---

# Distance Measurements

```python
distance mut_epi, mut, epi
color cyan, mut_epi
```

---

# Rendering Publication-Quality Figures

```python
set antialias, 2
set cartoon_transparency, 0.15

ray 2400,1800

png HA_I302S_consurf.png, dpi=300
```

---

# Example Complete PyMOL Workflow

```python
load consurf_grades.pdb, HA

hide everything
show cartoon
color gray80
bg_color white

select HA_A, chain A

select mut, HA_A and resi 302
show sticks, mut
color red, mut

label mut and name CA, "I302S"
set label_size, 18

select epi, byres (HA_A within 6 of mut)
show surface, epi
color yellow, epi
set transparency, 0.4

spectrum b, blue_white_red

orient

ray 2400,1800

png HA_I302S_consurf.png, dpi=300
```
---
---

# Software and Tools

| Tool | Purpose |
|---|---|
| Python | Sequence processing |
| Biopython | DNA to protein translation |
| MAFFT | Multiple sequence alignment |
| ConSurf | Conservation analysis |
| PyMOL | Structural visualization |
| RCSB PDB | Protein structure source |

---

# Applications

This workflow can be used for:

- Influenza mutation analysis
- Antigenic drift studies
- Conservation analysis
- Structural bioinformatics
- Computational virology
- Mutation interpretation
- Epitope mapping

---

# References

## PyMOL

https://pymol.org/

## MAFFT

https://mafft.cbrc.jp/alignment/software/

## ConSurf

https://consurf.tau.ac.il/

## RCSB Protein Data Bank

https://www.rcsb.org/

---

# License

This project is intended for academic and research purposes.
