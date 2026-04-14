# HELIUS: Baseline gut microbiota and incident cardiometabolic diagnoses

This repository contains the analysis code for a longitudinal study examining associations between baseline gut microbiota composition and incident cardiometabolic diagnoses (type 2 diabetes, hypertension, dyslipidemia) in a multi-ethnic urban cohort. Data are from the [HELIUS study](https://www.heliusstudy.nl/) (Healthy Life in an Urban Setting), Amsterdam, the Netherlands.

**Contact:** b.j.verhaar@amsterdamumc.nl

---

## Study overview

- **Cohort:** HELIUS — a multi-ethnic cohort including Dutch, South-Asian Surinamese, African Surinamese, Ghanaian, Turkish, and Moroccan participants
- **Design:** Longitudinal (baseline + follow-up), with 16S rRNA gut microbiota profiling at baseline and cardiometabolic outcomes assessed at follow-up
- **Outcomes:** New-onset type 2 diabetes, hypertension, and dyslipidemia
- **Microbiome data:** 16S rRNA amplicon sequencing, rarefied phyloseq objects (ASV-level)
- **Additional data:** Untargeted plasma metabolomics (Metabolon platform, paired baseline/follow-up)

---

## Repository structure

```
├── scripts/
├── data/       # Input data (not tracked in git)
└── results/    # Output figures and tables
```

---

## Scripts

| Script | Description |
|---|---|
| `1a_datacleaning.R` | Merges and cleans HELIUS clinical data (SPSS), constructs incident diagnosis flags (`DM_new`, `HT_new`, `Dyslip_new`), and prepares the rarefied phyloseq object matched to participants with follow-up data |
| `1b_metabolomics.R` | Cleans and log10-transforms Metabolon untargeted metabolomics data; aligns subject IDs across baseline and follow-up |
| `1c_datacleaning_allbaseline.R` | Prepares the complete baseline dataset including participants without follow-up, used for sensitivity analyses |
| `1d_rarefaction_curve.R` | Generates rarefaction curves to assess sequencing depth adequacy |
| `2_table1.R` | Descriptive statistics (Table 1 by timepoint/sex; ethnicity-stratified table; metabolomics subset table) using `tableone` |
| `3_diversity_ordination.R` | Alpha diversity (Shannon, richness, Faith's PD) and Bray-Curtis PCoA with PERMANOVA, stratified by incident diagnosis |
| `4a_microbes_diagnoses.R` | ASV-level logistic regression (CLR-transformed) for incident diabetes, hypertension, and dyslipidemia; two models (age-adjusted; fully adjusted); FDR correction; forest plots |
| `4b_microbes_diagnoses_ethnicity.R` | As 4a, with ethnicity interaction terms |
| `4c_microbes_diagnoses_sex.R` | As 4a, with sex interaction terms |
| `5_metabolites_microbes.R` | Linear regression of plasma metabolites on diagnosis-associated ASVs |

---

## Dependencies

The analyses use R (≥ 4.1). Package dependencies are managed with `renv`; restore the environment with:

```r
renv::restore()
```

Key packages:
- **Data handling:** `tidyverse`, `haven`, `rio`
- **Microbiome:** `phyloseq`, `vegan`, `picante`
- **Statistics:** `broom`, `tableone`, `compositions`
- **Visualization:** `ggplot2`, `ggpubr`, `patchwork`, `ggsci`, `ggthemes`, `ggrepel`

---

## Execution order

Scripts should be run in the order indicated by their filename numbering.
