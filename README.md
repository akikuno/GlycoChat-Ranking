# Immune Checkpoint Lectin Ranking in PDAC Single-Cell Analysis

## Overview

This repository contains a computational framework for identifying and ranking lectins related to immune checkpoint pathways in pancreatic ductal adenocarcinoma (PDAC) using single-cell multi-modal data analysis. The analysis integrates glycan binding patterns with RNA expression profiles to discover cell-type-specific lectin interactions that may play crucial roles in immune checkpoint regulation.

## Research Objective

The primary goal is to identify lectins that show specific binding patterns to immune cells in the PDAC tumor microenvironment, which could serve as:
- Novel immune checkpoint targets
- Biomarkers for immune cell subtypes
- Therapeutic targets for cancer immunotherapy

## Methodology

### Multi-Modal Data Integration
- **RNA Expression Data**: Single-cell RNA sequencing data capturing receptor gene expression
- **Glycan Binding Data**: Lectin microarray data measuring glycan-lectin interactions
- **Cell Types Analyzed**: 
  - Cancer cells (Classical, Basal-like, Intermediate)
  - Immune cells (T cells, TAMs, MDSCs, dendritic cells, B cells, etc.)

### Scoring Algorithm
The analysis employs a sophisticated scoring system that evaluates:
1. **RNA Specificity Score**: How specifically a lectin receptor is expressed in immune cells
2. **RNA Expression Level**: The magnitude of receptor expression in immune cells
3. **Glycan Specificity Score**: How specifically a lectin binds to cancer cells
4. **Glycan Binding Level**: The strength of lectin binding to cancer cells

### Key Features
- Emphasizes single cell-type specificity using coefficient of variation and ratio metrics
- Customizable weighting system for different scoring components
- Identifies top-expressing cell types for both RNA and glycan data

## Key Findings

The analysis identifies several lectins with strong cell-type specificity:
- **P-Selectin (SELP)**: Shows high specificity for certain immune cell populations
- **Langerin (CD207)**: Demonstrates specific binding patterns relevant to immune responses
- **LSECtin (CLEC4G)**: Exhibits cancer cell-specific interactions

## Installation

```bash
conda create -n scglycan python=3.12
conda install -y -n scglycan -c conda-forge r-essentials r-base r-seurat r-pheatmap r-patchwork r-ggplotify r-languageserver
conda activate scglycan
```

## Usage

1. Ensure the PDAC single-cell dataset (`pdac_ctype.RData`) is placed in the `data/` directory

>[!NOTE]
> The `pdac_ctype.RData` file should be requested from the corresponding author (Dr. Hiroaki Tateno: h-tateno[at]aist.go.jp)

2. Run the main analysis:
   ```r
   # In R or RStudio
   quarto::quarto_render("immune_checkpoint_lectin_ranking.qmd")
   ```
3. Results will be saved to:
   - `data/glycan_ranking.csv`: Ranked list of lectins with scores
   - `data/glycan_ranking_top10.png`: Visualization of top-ranked lectins

## Repository Structure

```
├── immune_checkpoint_lectin_ranking.qmd  # Main analysis notebook
├── immune_checkpoint_lectin_ranking.R    # Generated R script
├── scripts/
│   └── load.R                           # Data exploration script
└── data/
    ├── pdac_ctype.RData                # Input: Seurat object with multi-modal data
    └── glycan_ranking.csv              # Output: Lectin rankings
```

## Biological Significance

This analysis provides insights into:
- Lectin-mediated immune checkpoint mechanisms in PDAC
- Cell-type-specific glycosylation patterns in the tumor microenvironment
- Potential therapeutic targets for enhancing anti-tumor immunity

<!-- ## Citation

If you use this analysis in your research, please cite:

```
[Author information to be added]
``` -->

## License

This project is licensed under the MIT License - see the LICENSE file for details.
