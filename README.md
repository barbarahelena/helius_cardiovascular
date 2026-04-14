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
│   ├── 1a_datacleaning.R              # Clinical data cleaning + phyloseq preparation
│   ├── 1b_metabolomics.R              # Metabolomics data cleaning and preprocessing
│   ├── 1c_datacleaning_allbaseline.R  # Data cleaning for all baseline participants
│   ├── 1d_rarefaction_curve.R         # Rarefaction curve generation
│   ├── 2_table1.R                     # Descriptive statistics (Table 1 and supplementary tables)
│   ├── 3_diversity_ordination.R       # Alpha diversity and beta diversity (Bray-Curtis PCoA)
│   ├── 4a_microbes_diagnoses.R        # ASV-level associations with new diagnoses (overall)
│   ├── 4b_microbes_diagnoses_ethnicity.R  # Interaction analyses by ethnicity
│   ├── 4c_microbes_diagnoses_sex.R    # Interaction analyses by sex
│   └── 5_metabolites_microbes.R       # Associations between microbiota and metabolomics
├── data/                              # Input data (not tracked in git; see Data section below)
├── results/                           # Output figures and tables
└── helius_cardiovascular.Rproj        # R project file
```

---

## Scripts

### 1a — Data cleaning (`1a_datacleaning.R`)

Merges and cleans multiple HELIUS data extracts (SPSS `.sav` files). Key steps:
- Selects and renames variables across domains: demographics, physical examination, cardiometabolic diagnoses, medications, laboratory values, dietary intake, fecal sample characteristics
- Recodes Dutch-language factor levels to English
- Constructs derived variables: dyslipidemia (high TG or low HDL), delta variables (change from baseline to follow-up), incident diagnosis flags (`DM_new`, `HT_new`, `Dyslip_new`)
- Excludes participants who used antibiotics at baseline
- Pivots data to both wide (one row per participant, baseline + follow-up columns) and long (one row per timepoint) formats
- Loads rarefied 16S phyloseq object, filters to bacterial taxa and to participants with matched microbiome data, and saves cleaned phyloseq and taxonomy table
- Exports a baseline ASV table, taxonomy table, reference sequences (FASTA), and phylogenetic tree for analyses in the CBS (Statistics Netherlands) secure environment

**Output:** `data/clinicaldata_wide.RDS`, `data/clinicaldata_long.RDS`, `data/phyloseq_rarefied_cleaned.RDS`, `data/taxtable_rarefied_cleaned.RDS`

---

### 1b — Metabolomics cleaning (`1b_metabolomics.R`)

Processes untargeted metabolomics data from Metabolon:
- Maps Metabolon IDs to chemical names using the provided info file
- Aligns subject IDs across baseline and COVID follow-up samples
- Removes metabolites with zero variance
- Log10-transforms abundance values with pseudocount
- Exports paired baseline/follow-up metabolomics matrix

---

### 1c — All-baseline data cleaning (`1c_datacleaning_allbaseline.R`)

Prepares the complete baseline dataset (all participants with 16S data, including those without follow-up), used for sensitivity analyses.

---

### 1d — Rarefaction curve (`1d_rarefaction_curve.R`)

Generates rarefaction curves to assess sequencing depth adequacy across samples, using the unrarefied phyloseq object.

---

### 2 — Descriptive statistics (`2_table1.R`)

Produces:
- **Table 1:** Baseline and follow-up characteristics stratified by timepoint and sex
- **Supplementary Table (ethnicity-stratified):** Baseline characteristics by ethnic group, including incident diagnosis rates
- **Metabolomics subset table:** Characteristics of the subset with paired metabolomics data

Uses the `tableone` package. Non-normally distributed variables (SCORE CVD mortality, triglycerides) are reported as median (IQR).

---

### 3 — Diversity and ordination (`3_diversity_ordination.R`)

#### Alpha diversity
Calculates three metrics from the rarefied phyloseq object:
- **Shannon index** (diversity)
- **Species richness** (number of observed ASVs)
- **Faith's phylogenetic diversity** (PD)

Compares each metric between participants with and without incident diagnoses (Wilcoxon test), producing violin + boxplot figures.

#### Beta diversity (Bray-Curtis)
For each diagnosis separately:
- Computes Bray-Curtis dissimilarity matrix
- Performs principal coordinate analysis (PCoA) with Cailliez correction
- Runs PERMANOVA (`adonis2`) to test for community-level differences
- Plots PCoA with 95% confidence ellipses annotated with PERMANOVA R² and p-value

**Output:** `results/alphadiversity/`, `results/ordination/`, `results/diversityplots_diagnoses.pdf`

---

### 4a — Microbiota–diagnosis associations (`4a_microbes_diagnoses.R`)

ASV-level logistic regression for each of the three incident diagnoses:
- CLR (centered log-ratio) transformation of ASV counts (with pseudocount 0.5) prior to filtering
- ASVs retained if present at >5 reads in >20% of samples
- **Model 1:** ASV + age
- **Model 2:** ASV + age + sex + BMI + smoking + alcohol

Results filtered to ASVs nominally significant in Model 1, then FDR-corrected within each diagnosis. Forest plots display ORs with 95% CIs, colored by bacterial family, with shape indicating FDR significance.

**Output:** `results/diabetesassociations_smoking_alcohol_bmi.csv`, `results/dyslipidemiaassociations_smoking_alcohol_bmi.csv`, `results/hypertensionassociations_smoking_alcohol_bmi.csv`, `results/microbiome_diagnoses/diagnosispred_long_adjusted.pdf`

---

### 4b — Ethnicity interaction (`4b_microbes_diagnoses_ethnicity.R`)

Repeats the ASV-level logistic regression from script 4a with ethnicity interaction terms, testing whether microbiota–diagnosis associations differ across ethnic groups.

---

### 4c — Sex interaction (`4c_microbes_diagnoses_sex.R`)

Repeats the ASV-level logistic regression from script 4a with sex interaction terms.

---

### 5 — Metabolomics–microbiota associations (`5_metabolites_microbes.R`)

Examines whether ASVs associated with incident diagnoses (from the CBS secure environment analyses: MACE, MACEplus) also associate with plasma metabolites:
- Selects diagnosis-associated ASVs (FDR < 0.05, Model 1)
- Log10-transforms metabolomics data
- Runs linear regression of each metabolite on each selected ASV
- Produces volcano plots or summary figures of significant associations

---

## Data

Raw data files are not included in this repository due to privacy regulations (HELIUS participant data). Data access can be requested through the HELIUS study management at Amsterdam UMC.

Expected data files:
| File | Description |
|---|---|
| `data/240411_HELIUS data Barbara Verhaar.sav` | Main HELIUS clinical data extract |
| `data/210517_HELIUS data Ulrika Boulund_2.sav` | PPI, metformin, statin variables |
| `data/220712_HELIUS data Barbara Verhaar - Cov1 groep en datum.sav` | COVID-19 follow-up dates |
| `data/16s/phyloseq/rarefied/phyloseq_rarefied.RDS` | Rarefied 16S phyloseq object |
| `data/16s/phyloseq/rarefied/taxtable_rarefied.RDS` | Taxonomy table |
| `data/metabolomics/HELIUS_metabolomics_abundance.xlsx` | Metabolon abundance matrix |
| `data/metabolomics/HELIUS_metabolomics_metadata.xlsx` | Metabolomics sample metadata |
| `data/metabolomics/Info_HELIUS_metabolomics.xlsx` | Metabolon compound information |

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

Run scripts in numbered order:

```
1a → 1b → 1c → 1d → 2 → 3 → 4a → 4b → 4c → 5
```

Scripts 1c and 1d are independent of 1b and can be run in parallel with it. Scripts 4b and 4c are independent of each other and can be run in parallel with 4a.
