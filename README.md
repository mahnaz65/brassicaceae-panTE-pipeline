# brassicaceae-panTE-pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-green.svg)](https://www.python.org/)
[![R 4.1+](https://img.shields.io/badge/R-4.1+-276DC3.svg)](https://www.r-project.org/)

Reproducible pan-TE annotation, clustering, and comparative phylogenomic scaling pipeline across **13 Brassicaceae genomes** (*Aubrieta deltoidea*, *Bornmuellera tymphaea*, *Draba lasiocarpa*, *Iberis sempervirens*, *Isatis tinctoria*, *Rorippa pyrenaica*, *Aethionema saxatile*, *Alyssum gmelinii*, *Arabidopsis thaliana*, *Cardamine glauca*, *Erysimum linariifolium*, *Noccaea praecox*, and *Odontarrhena muralis*) using EDTA, MMseqs2, RepeatMasker, Python, and R.

---

## Overview

Single-species *de novo* repeat libraries introduce discovery bias when comparing repeat landscapes across divergent plant genera. This pipeline implements a standardized **panEDTA + MMseqs2** workflow followed by Smith–Waterman-prioritized overlap resolution and comparative phylogenomic analyses:

1. **Step 1 (`scripts/01_run_edta.sh`):** Per-species structural *de novo* TE discovery using `EDTA v2.2.0` (with optional CDS filtering).
2. **Step 2 (`scripts/02_build_pan_library.sh`):** Cross-species consensus clustering (`80-80-80` rule) with `MMseqs2` to build a non-redundant 42,688-family pan-Brassicaceae TE library.
3. **Step 3 (`scripts/03_run_repeatmasker.sh`):** Uniform genome-wide annotation across all 13 assemblies with `RepeatMasker v4.1.2-p1`.
4. **Step 4 (`scripts/04_parse_nonoverlapping_te.py`):** Per-scaffold interval deduplication (resolving nested/overlapping hits by highest Smith–Waterman alignment score) and quantitative abundance filtering (`39,785` primary families vs. `42,688` unfiltered families).
5. **Step 5 (`scripts/05_phylogeny_and_scaling.R`):** Outgroup-rooted Neighbor-Joining (NJ) phylogenetic reconstruction (Hellinger distance, 500 bootstrap replicates), genome-size scaling regressions (all $N = 13$ species and $<500$ Mbp sensitivity subset $N = 10$), and generation of publication figures and summary tables.

---

## Repository Structure

```text
brassicaceae-panTE-pipeline/
├── README.md
├── LICENSE
├── CITATION.cff
├── environment.yml
├── data/
│   └── species_metadata.tsv              # Assembly accession/filename mapping & taxonomy
├── scripts/
│   ├── 01_run_edta.sh                    # Step 1: De novo TE discovery per assembly (EDTA)
│   ├── 02_build_pan_library.sh           # Step 2: Pan-library clustering (MMseqs2)
│   ├── 03_run_repeatmasker.sh            # Step 3: Uniform RepeatMasker annotation
│   ├── 04_parse_nonoverlapping_te.py     # Step 4: Overlap resolution & abundance matrices
│   └── 05_phylogeny_and_scaling.R        # Step 5: NJ tree, Pearson scaling & figures
└── results/
    ├── 04_abundance_matrices/            # Lightweight summary matrices (included for fast reproduction)
    └── 05_figures_and_trees/             # Generated trees, summary tables, and PDF/PNG figures\

## Installation & Dependencies
Create the Conda environment containing all required bioinformatics, Python, and R packages:\
```
conda env create -f environment.yml
conda activate brassicaceae-pante
```
Core Software Versions Used:

EDTA v2.2.0

MMseqs2 v13.45111

RepeatMasker v4.1.2-p1

Python v3.8+ (numpy, pandas, scipy)

R v4.1.3+ (ape, ggplot2, dplyr, tidyr, patchwork)

Data Availability
Because raw genome assemblies and RepeatMasker alignment files exceed GitHub's 100 MB file-size limit, large input and intermediate datasets are archived externally:

Genome assemblies & CDS annotations: NCBI / ENA BioProject [ACCESSION] (see data/species_metadata.tsv)

Curated 42,688-family pan-TE library & RepeatMasker .out alignments: Zenodo DOI [10.5281/zenodo.XXXXXXX]\
