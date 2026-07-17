"""
PyMOL Commands Used for Influenza Structural Visualization
==========================================================

This file contains PyMOL commands used throughout the workflow for:

- HA and NA visualization
- Mutation mapping
- Antigenic site visualization
- ConSurf integration
- Publication-quality rendering

Author: Generated for influenza structural mapping workflow
"""

# =========================================================
# BASIC STRUCTURE LOADING
# =========================================================

# Load HA structure
load 8g5b.pdb, HA

# Load NA structure
load 3ti6.pdb, NA

# Fetch structures directly from PDB
fetch 1RU7, HA
fetch 3TI6, NA

# =========================================================
# BASIC VISUALIZATION
# =========================================================

hide everything
show cartoon

color gray80
color gray70, HA
color gray50, NA

bg_color white

# =========================================================
# CHAIN SELECTION
# =========================================================

# Select chain A
select HA_A, chain A

# Show only selected chain
hide everything
show cartoon, HA_A
color gray80, HA_A

# =========================================================
# MUTATION VISUALIZATION
# =========================================================

# Example mutation: I302S
select mut, HA_A and resi 302

# Other example mutations
select HA_mut, HA and resi 145
select NA_mut, NA and resi 275

# Show mutation residues
show sticks, mut
show sticks, HA_mut
show sticks, NA_mut

# Color mutations
color red, mut
color red, HA_mut
color blue, NA_mut

# Zoom into mutation
zoom mut, 8

# =========================================================
# LABELING
# =========================================================

label mut and name CA, "I302S"
label HA_mut and name CA, "HA mutation"
label NA_mut and name CA, "NA mutation"

set label_size, 18
set label_size, 20

set label_color, black
set label_color, white

# =========================================================
# ANTIGENIC SITE VISUALIZATION
# =========================================================

# Automatic nearby antigenic region
select epi, byres (HA_A within 6 of mut)

# HA antigenic neighborhood
select HA_epi, byres (HA within 6 of HA_mut)

# NA antigenic neighborhood
select NA_epi, byres (NA within 6 of NA_mut)

# Example literature-based antigenic residues
select siteA, resi 121+122+124+133+135+137

# Manual antigenic site example
select HA_epi, HA and resi 140-150+155+158
select NA_epi, NA and resi 270-285

# =========================================================
# SURFACE VISUALIZATION
# =========================================================

show surface, epi
show surface, HA_epi
show surface, NA_epi
show surface, siteA

# Surface coloring
color yellow, epi
color yellow, HA_epi
color orange, NA_epi
color magenta, siteA

# Transparency
set transparency, 0.4
set transparency, 0.4, epi
set transparency, 0.4, HA_epi
set transparency, 0.4, NA_epi

# =========================================================
# DISTANCE MEASUREMENTS
# =========================================================

distance mut_epi, mut, epi
distance HA_dist, HA_mut, HA_epi
distance NA_dist, NA_mut, NA_epi

# Color distance lines
color cyan, mut_epi
color cyan, HA_dist
color magenta, NA_dist

# =========================================================
# ORIENTATION AND POSITIONING
# =========================================================

orient

# Separate HA and NA structures
translate [-60,0,0], HA
translate [60,0,0], NA

# WT vs MUT comparison
translate [-30,0,0], WT
translate [30,0,0], MUT

# =========================================================
# TRANSPARENCY SETTINGS
# =========================================================

set cartoon_transparency, 0.2
set cartoon_transparency, 0.15

# =========================================================
# IMAGE QUALITY SETTINGS
# =========================================================

set antialias, 2
set ray_trace_mode, 1

# =========================================================
# RENDERING
# =========================================================

ray 2000,1500
ray 2400,1800

# Save publication-quality figures
png HA_mutation.png, dpi=300
png HA_I302S.png, dpi=300
png HA_I302S_consurf.png, dpi=300

# =========================================================
# CONSURF INTEGRATION
# =========================================================

# Load ConSurf-colored structure
load consurf_grades.pdb

# Apply conservation coloring
spectrum b, blue_white_red

# Highlight mutation on conservation structure
select mut, resi 302
show sticks, mut
color yellow, mut

# =========================================================
# WILDTYPE VS MUTANT COMPARISON
# =========================================================

create WT, HA_A
create MUT, HA_A

# =========================================================
# MUTAGENESIS WIZARD WORKFLOW
# =========================================================

"""
GUI Workflow:

1. Wizard -> Mutagenesis
2. Select residue
3. Choose amino acid
4. Select rotamer
5. Apply
6. Done
"""

# =========================================================
# COMPLETE HA WORKFLOW
# =========================================================

"""
load 8g5b.pdb, HA
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

orient
ray 2400,1800
png HA_I302S.png, dpi=300
"""

# =========================================================
# COMPLETE HA + NA WORKFLOW
# =========================================================

"""
fetch 1RU7, HA
fetch 3TI6, NA

hide everything
show cartoon, HA
show cartoon, NA

color gray70, HA
color gray50, NA

select HA_mut, HA and resi 145
select NA_mut, NA and resi 275

show sticks, HA_mut
show sticks, NA_mut

color red, HA_mut
color blue, NA_mut

label HA_mut and name CA, "HA mutation"
label NA_mut and name CA, "NA mutation"

select HA_epi, byres (HA within 6 of HA_mut)
select NA_epi, byres (NA within 6 of NA_mut)

show surface, HA_epi
show surface, NA_epi

color yellow, HA_epi
color orange, NA_epi

translate [-60,0,0], HA
translate [60,0,0], NA

ray 2400,1800
"""
