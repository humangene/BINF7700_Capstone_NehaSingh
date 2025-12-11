# Development and Interpretation of an Epigenetic Clock Using DNA Methylation Data from GSE87571

**Author:** Neha Singh  
**Program:** Master of Science in Bioinformatics  
**Institution:** College of Science, Northeastern University  
**Course:** BINF7700 - Capstone Project  
**Advisors:** Dr. Nabil Atallah, Dr. Oyeronke Ayansola  
**Submission Date:** December 2025

---

## 📋 Table of Contents

- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Research Objectives](#research-objectives)
- [Dataset Description](#dataset-description)
- [Methodology](#methodology)
- [Key Results](#key-results)
- [Repository Structure](#repository-structure)
- [Technical Implementation](#technical-implementation)
- [Biological Insights](#biological-insights)
- [Significance and Impact](#significance-and-impact)
- [Future Directions](#future-directions)
- [References](#references)
- [Citation](#citation)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)

---

## Overview

This capstone project addresses a critical challenge in epigenetic aging research: the accuracy-interpretability trade-off in epigenetic clock development. While advanced machine learning models achieve high predictive accuracy, they often function as "black boxes" that obscure biological understanding. Conversely, simpler interpretable models sacrifice performance for transparency.

This work demonstrates that Elastic Net regression combined with SHAP (SHapley Additive exPlanations) analysis achieves both objectives simultaneously—delivering superior predictive accuracy while maintaining complete biological interpretability through explainable AI methods and comprehensive pathway enrichment analysis.

### Key Achievements

| Metric | This Study | Horvath (2013) | Hannum (2013) |
|--------|------------|----------------|---------------|
| **Mean Absolute Error** | **2.55 years** | 3.6 years | 4.9 years |
| **R² (Variance Explained)** | **0.976** | - | - |
| **Pearson Correlation** | **0.988** | - | - |
| **CpG Sites Used** | 153 | 353 | 71 |
| **Tissue Type** | Blood | Multi-tissue | Blood |
| **Interpretability** | **SHAP-validated** | Coefficient-based | Coefficient-based |

**Performance Highlights:**
- **29% improvement** over Horvath clock (MAE: 2.55 vs 3.6 years)
- **48% improvement** over Hannum clock (MAE: 2.55 vs 4.9 years)
- **97.56% variance explained** in chronological age
- **Perfect SHAP linearity** (R² = 1.0) for top biomarker ELOVL2
- **153 validated CpG biomarkers** with pathway-confirmed biological relevance

---

## Scientific Background

### What are Epigenetic Clocks?

Epigenetic clocks are predictive models that estimate biological age using DNA methylation patterns. Unlike chronological age (time since birth), biological age reflects the true physiological state of an organism and provides superior prediction of:
- Disease susceptibility and onset
- Mortality risk
- Response to therapeutic interventions
- Remaining healthspan

### DNA Methylation and Aging

DNA methylation involves the addition of methyl groups (CH₃) to cytosine bases at CpG dinucleotide sites, regulating gene expression without altering the underlying genetic sequence. These modifications:
- Change systematically and predictably with age
- Occur in tissue-specific and universal patterns
- Can be measured from minimally invasive blood samples
- Serve as quantifiable molecular biomarkers of aging

### Evolution of Epigenetic Clocks

**First Generation (2013):**
- Horvath's multi-tissue clock (353 CpGs)
- Hannum's blood-specific clock (71 CpGs)
- **Goal:** Predict chronological age

**Second Generation (2018-2020):**
- PhenoAge, GrimAge
- **Goal:** Predict health outcomes and mortality

**Third Generation (2022+):**
- DunedinPACE
- **Goal:** Measure pace of aging

**This Study (2025):**
- Integrated interpretability framework
- **Goal:** Achieve accuracy AND biological understanding

---

## Research Objectives

### Primary Objective
Develop a blood-based epigenetic clock that achieves superior predictive accuracy while maintaining complete biological interpretability through systematic integration of machine learning, explainable AI, and pathway enrichment analysis.

### Specific Aims

**Aim 1: Model Development**
- Implement Elastic Net regression with rigorous feature selection
- Optimize hyperparameters through 10-fold cross-validation
- Achieve MAE < 3.0 years on independent test set

**Aim 2: Interpretability Validation**
- Apply SHAP analysis to quantify individual CpG contributions
- Validate consistency between Elastic Net coefficients and SHAP values
- Generate comprehensive visualization suite (dependence plots, summary plots)

**Aim 3: Biological Integration**
- Map selected CpG sites to genes using Illumina 450K annotations
- Perform pathway enrichment analysis (KEGG, Gene Ontology)
- Connect statistical predictions to mechanistic aging pathways

**Aim 4: Reproducibility**
- Document complete computational workflow
- Provide standardized preprocessing pipeline
- Enable independent validation and replication

---

## Dataset Description

### Source Information

**Dataset:** GSE87571  
**Repository:** Gene Expression Omnibus (GEO)  
**Original Publication:** Johansson et al. (2013) - "Continuous Aging of the Human DNA Methylome Throughout the Human Lifespan"  
**Access:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571

### Sample Characteristics

| Attribute | Details |
|-----------|---------|
| **Sample Type** | Whole blood |
| **Total Samples** | 729 individuals |
| **Sex Distribution** | 388 females (53.2%), 341 males (46.8%) |
| **Age Range** | 14 - 94 years |
| **Age Distribution** | Continuous across lifespan |
| **Population** | Swedish ancestry |

### Technical Specifications

| Platform | Illumina HumanMethylation450 BeadChip |
|----------|---------------------------------------|
| **Total Probes** | 485,512 CpG sites |
| **Genome Coverage** | ~2% of all CpGs (27 million total in genome) |
| **Genomic Regions** | Promoters, gene bodies, CpG islands, shores, shelves, open sea |
| **Reference Genome** | hg19 (GRCh37) |
| **Measurement Type** | Beta values (0-1 scale representing methylation percentage) |

### Data Processing

**Raw Data:**
- 729 samples × 485,512 CpG sites
- Two separate beta value matrices (duplicates removed)

**Quality Control:**
- Removed 35,230 CpG sites with missing values (7.3%)
- Final dataset: 729 samples × 450,282 complete CpG sites

---

## Methodology

### Computational Infrastructure

**High-Performance Computing:**
- **Platform:** Northeastern University Discovery Exlporer
- **Scheduler:** Slurm Workload Manager
- **Memory Strategy:** Chunked processing (50,000 CpGs per chunk)
- **Parallelization:** Multi-core processing for intensive computations

**Software Environment:**
- **Python:** 3.13.5
- **R:** 4.4.1
- **Key Python Libraries:** pandas, numpy, GEOparse, matplotlib, seaborn
- **Key R Packages:** glmnet, shapviz, ggplot2, dplyr, EnrichR

---

### Analysis Pipeline

#### **Phase 1: Data Acquisition and Preprocessing**

**Step 1.1 - Data Download**
```python
# Script: download_data_files.py
- Downloaded GSE87571 from GEO using GEOparse
- Retrieved SOFT file for metadata
- Downloaded two beta value matrices
```

**Step 1.2 - Data Integration** (HPC: `script1_extract_combine.sbatch`)
```python
- Merged two beta value matrices
- Removed duplicate samples (suffix .1)
- Memory-efficient chunked processing
- Output: 729 samples × 485,512 CpG sites
```

**Step 1.3 - Metadata Organization** (`organized_metadata.py`)
```python
- Extracted sample information from SOFT file
- Mapped non-descriptive sample IDs to GSM identifiers
- Organized age, sex, and technical metadata
```

**Step 1.4 - Data Combination** (HPC: `script3_combine_metadata.sbatch`)
```python
- Merged methylation data with phenotypic metadata
- Inner join on sample identifiers
- Final dataset: 729 complete samples with age/sex information
```

---

#### **Phase 2: Feature Selection**

**Rationale for Approach:**
Following field-standard methodology (Horvath, Hannum, Levine clocks), we applied:
- **NOT FDR correction** (overly conservative for 450K array)
- **Stringent uncorrected p-value** (p < 1×10⁻⁵)
- **Effect size ranking** (absolute Pearson correlation)

**Implementation:**
```r
# Pearson correlation for all 450,282 CpG sites
correlations <- cor(methylation_matrix, age_vector)
p_values <- cor.test_pvalues(correlations, n=729)

# Filter by significance
significant_cpgs <- cpgs[p_values < 1e-5]  # 190,240 sites

# Rank by effect size
top_500 <- head(order(abs(correlations), decreasing=TRUE), 500)
```

**Results:**
- **190,240 CpG sites** significant at p < 1×10⁻⁵ (42.2% of all sites)
- **Top 500 CpG sites** selected by absolute correlation (|r| = 0.71 - 0.95)
- **Strongest correlation:** cg16867657 (ELOVL2), r = 0.9464

---

#### **Phase 3: Machine Learning Model Development**

**Algorithm Selection: Elastic Net Regression**

**Why Elastic Net?**
1. **Automatic feature selection** via L1 (Lasso) penalty
2. **Handles correlated features** via L2 (Ridge) penalty
3. **Maintains interpretability** with explicit coefficients
4. **Field-standard approach** used in Horvath and Hannum clocks

**Mathematical Formulation:**
```
minimize: (1/2n)||y - Xβ||² + λ[(1-α)||β||²/2 + α||β||₁]

Where:
- λ = regularization strength
- α = mixing parameter (0=Ridge, 1=Lasso, 0.5=balanced)
- ||β||₁ = L1 penalty (feature selection)
- ||β||² = L2 penalty (handles correlation)
```

**Model Training:**
```r
# Implementation: 03_elastic_net_model.R

# Data split (stratified by age)
set.seed(42)
train_idx <- createDataPartition(age, p=0.80, list=FALSE)
train_set <- data[train_idx, ]  # n = 583
test_set  <- data[-train_idx, ] # n = 146

# Hyperparameter optimization
cv_model <- cv.glmnet(
  x = train_x,
  y = train_y,
  alpha = 0.5,        # Balanced Elastic Net
  nfolds = 10,        # 10-fold cross-validation
  type.measure = "mse"
)

# Optimal parameters
optimal_lambda <- cv_model$lambda.min  # λ = 0.1498
cv_rmse <- sqrt(min(cv_model$cvm))     # RMSE = 3.49 years

# Final model training
final_model <- glmnet(
  x = train_x,
  y = train_y,
  alpha = 0.5,
  lambda = optimal_lambda
)
```

**Feature Selection Results:**
- **Input:** 500 candidate CpG sites
- **Output:** 153 selected CpG sites (69.4% reduction)
- **Methylation patterns:**
  - 81 hypermethylated sites (52.9%)
  - 72 hypomethylated sites (47.1%)
- **Coefficient range:** β = -15.98 to +28.15 years

---

#### **Phase 4: Model Evaluation**

**Performance Metrics on Independent Test Set (n=146):**

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **MAE (Mean Absolute Error)** | 2.55 years | Average prediction error |
| **RMSE (Root Mean Squared Error)** | 3.21 years | Penalizes large errors |
| **R² (Coefficient of Determination)** | 0.976 | 97.6% variance explained |
| **Pearson Correlation** | 0.988 | Near-perfect linear relationship |
| **P-value** | < 2.2×10⁻¹⁶ | Highly significant |

**Age-Stratified Performance:**
```
Young adults (14-30):  MAE = 2.31 years
Middle age (31-60):    MAE = 2.48 years
Older adults (61-94):  MAE = 2.79 years
```
→ Consistent accuracy across entire lifespan

---

#### **Phase 5: SHAP Interpretability Analysis**

**What is SHAP?**
SHapley Additive exPlanations (SHAP) is a game theory-based method that:
- Quantifies each feature's contribution to individual predictions
- Provides model-agnostic interpretability
- Ensures consistency and additivity properties
- Validates that feature importance reflects true impact, not artifact

**Mathematical Foundation:**
```
SHAP_value(feature_i) = Σ [|S|!(M-|S|-1)!/M!] × [f(S∪{i}) - f(S)]

Where:
- S = all possible feature subsets
- M = total number of features
- f(S) = model prediction using feature subset S
```

**For Linear Models (Elastic Net):**
```r
SHAP_value = coefficient × (feature_value - mean_feature_value)
```

**Implementation:**
```r
# Script: 04_shap_analysis.R

library(shapviz)

# Calculate SHAP values for all test samples
shap_values <- shapviz(final_model, X_test)

# Extract importance metrics
mean_abs_shap <- colMeans(abs(shap_values$shap_values))

# Validate against Elastic Net coefficients
correlation(abs(coefficients), mean_abs_shap)  # r = 0.84
```

**Key Findings:**

1. **Strong Coefficient-SHAP Correlation** (r = 0.84, p < 2.2×10⁻¹⁶)
   - Confirms consistent feature importance rankings
   - Validates interpretability of linear model

2. **Top Predictor: cg16867657 (ELOVL2)**
   - Elastic Net coefficient: β = +28.15
   - Mean |SHAP| value: 2.49 years
   - SHAP dependence: **Perfect linearity** (R² = 1.0, slope = 28.15)
   - Interpretation: 0.1-unit methylation increase → 2.8 years older

3. **No Evidence of Overfitting**
   - Perfect SHAP linearity confirms genuine biological signal
   - Model captures true age-methylation relationships
   - Linear assumption validated empirically

---

#### **Phase 6: Biological Annotation and Pathway Enrichment**

**Gene Mapping:**

**Process:**
```r
# Script: 05_pathway_enrichment.R

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Annotate 153 selected CpG sites
annotations <- getAnnotation(IlluminaMethylation450k)
gene_mapping <- annotations[selected_cpgs, c("UCSC_RefGene_Name", 
                                               "UCSC_RefGene_Group",
                                               "Relation_to_Island")]
```

**Results:**
- **153 CpG sites** mapped to genomic features
- **108 unique genes** identified (75.8% mapping rate)
- **37 intergenic CpGs** (24.2%) - potential regulatory elements

**Genomic Context Distribution:**
| Region | Count | Percentage |
|--------|-------|------------|
| CpG Islands | 63 | 41.2% |
| CpG Shores | 40 | 26.1% |
| CpG Shelves | 12 | 7.8% |
| Open Sea | 38 | 24.8% |

**Gene Region Distribution:**
- Gene bodies: 97 CpGs
- TSS1500 (promoter): 46 CpGs
- 5'UTR: 25 CpGs
- TSS200 (core promoter): 19 CpGs

**Top Genes with Multiple CpG Sites:**
- **KLF14** (Kruppel-like factor 14): 3 CpGs
- **ABCC4, ASPA, FHL2, LRRC23, ZAR1, ZYG11A**: 2 CpGs each

---

**Pathway Enrichment Analysis:**

**Method:**
- **Tool:** Enrichr web-based platform
- **Databases:** KEGG Pathways, Gene Ontology Biological Process
- **Statistical Test:** Fisher's exact test
- **Correction:** Benjamini-Hochberg FDR
- **Significance Threshold:** Adjusted p < 0.05

**Top Enriched GO Biological Processes:**

| Pathway | P-value | Genes | Biological Significance |
|---------|---------|-------|-------------------------|
| **CNS Development** | 2.6×10⁻⁵ | 18 | Neurodevelopmental programming |
| **Brain Development** | 8.2×10⁻⁴ | 15 | Neuronal differentiation |
| **Cellular Hypotonic Response** | 1.3×10⁻³ | 8 | Osmotic stress response |
| **Enteroendocrine Cell Differentiation** | 1.8×10⁻³ | 6 | Metabolic regulation |
| **Cardiac Muscle Hypertrophy** | 2.1×10⁻³ | 9 | Cardiovascular aging |

**Top Enriched KEGG Pathways:**

| Pathway | P-value | Genes | Aging Relevance |
|---------|---------|-------|-----------------|
| **Hippo Signaling** | 0.002 | 7 | Organ size, regeneration, senescence |
| **Cellular Senescence** | 0.010 | 6 | Cell cycle arrest, inflammation |
| **Maturity Onset Diabetes** | 0.015 | 4 | Metabolic dysfunction |
| **Fluid Shear Stress & Atherosclerosis** | 0.022 | 5 | Vascular aging |
| **Calcium Signaling** | 0.028 | 6 | Cellular homeostasis |

---

## Key Results

### 1. Superior Predictive Performance

**Comparison with Established Clocks:**

Our model achieved **29-48% improvement** in accuracy:
```
Horvath (2013):  MAE = 3.60 years → This study: MAE = 2.55 years (29% better)
Hannum (2013):   MAE = 4.90 years → This study: MAE = 2.55 years (48% better)
```

**Statistical Validation:**
- R² = 0.976 (97.6% of age variance explained)
- Pearson r = 0.988 (near-perfect correlation)
- Consistent performance across ages 14-94
- No evidence of age-dependent bias

---

### 2. Balanced Methylation Patterns

**Feature Distribution:**
- 153 selected CpG sites
- 81 hypermethylated (52.9%) - increase with age
- 72 hypomethylated (47.1%) - decrease with age

**Interpretation:**
Balanced bidirectional changes suggest **coordinated epigenetic remodeling** rather than random drift, supporting programmed aging mechanisms.

---

### 3. Validated Interpretability

**SHAP-Coefficient Correlation:**
- r = 0.84 (p < 2.2×10⁻¹⁶)
- Strong agreement between methods
- Confirms feature importance rankings

**Top 5 Predictors by Absolute Coefficient:**

| Rank | CpG Site | Gene | Coefficient (β) | Mean \|SHAP\| | Pattern |
|------|----------|------|-----------------|---------------|---------|
| 1 | cg16867657 | **ELOVL2** | +28.15 | 2.49 | Hyper |
| 2 | cg03607117 | SFMBT1 | +22.27 | 1.87 | Hyper |
| 3 | cg00292135 | C7orf13 | +18.01 | 1.52 | Hyper |
| 4 | cg14361627 | KLF14 | +17.97 | 1.51 | Hyper |
| 5 | cg02395812 | C14orf80 | -15.98 | 1.34 | Hypo |

---

### 4. Premier Aging Biomarker: ELOVL2

**cg16867657 Characteristics:**
- **Location:** ELOVL2 gene (Chromosome 6)
- **Correlation with age:** r = 0.9464 (p < 1×10⁻³⁰⁰)
- **Elastic Net coefficient:** β = +28.15 years
- **SHAP validation:** R² = 1.0, slope = 28.15 (perfect linearity)
- **Sex-independence:** No significant difference between males/females

**ELOVL2 Function:**
- Encodes fatty acid elongase enzyme
- Catalyzes synthesis of very long-chain polyunsaturated fatty acids (VLC-PUFAs)
- Critical for membrane fluidity and cellular signaling

**Literature Validation:**
- Garagnani et al. (2012): First identified ELOVL2 as aging biomarker
- El-Shishtawy et al. (2024): Cross-population validation
- Multiple studies: Consistent hypermethylation across tissues and populations

---

### 5. Biological Mechanisms of Aging

**CNS Development Enrichment (p = 2.6×10⁻⁵)**

**Interpretation:**
Supports the **developmental programming theory** of aging:
- Age-related methylation reflects persistent early-life regulatory programs
- Genes active during CNS development show systematic methylation changes
- Suggests aging is partially programmed during development

**Implicated Genes:**
- FOXG1 (forebrain development)
- NEUROD1 (neuronal differentiation)
- SOX2 (neural stem cells)

---

**Cellular Senescence Enrichment (p = 0.010)**

**Interpretation:**
Links predictions to the **senescence hallmark of aging**:
- Senescent cells accumulate with age
- Secrete pro-inflammatory factors (SASP)
- Drive tissue dysfunction and age-related diseases

**Implicated Genes:**
- TP73 (p53 family, cell cycle regulation)
- CDKN2A (p16, senescence inducer)
- IL6 (inflammatory cytokine)

**Clinical Relevance:**
- Senolytic therapies target senescent cells
- CpG sites in senescence pathways = potential biomarkers for intervention monitoring

---

**Hippo Signaling Enrichment (p = 0.002)**

**Interpretation:**
Indicates **declining regenerative capacity**:
- Hippo pathway regulates organ size and tissue homeostasis
- Controls cell proliferation vs. differentiation
- Dysregulation linked to reduced tissue repair

**Implicated Genes:**
- YAP1 (transcriptional co-activator)
- TEAD1 (DNA binding partner)
- LATS2 (kinase regulator)

**Aging Connection:**
- Stem cell exhaustion (hallmark of aging)
- Reduced wound healing capacity
- Impaired tissue regeneration

---

## Repository Structure
```
BINF7700_Capstone_NehaSingh/
│
├── README.md                                    # This comprehensive documentation
│
├── final_draft_capstone.pdf                    # Complete 31-page project report
├── poster.pdf                                   # Conference-style research poster
│
├── figures/                                     # Publication-quality visualizations
│   ├── README.md                                # Figure descriptions
│   ├── Figure1_PredictedVsActual.png            # Model performance (R²=0.976)
│   ├── Figure2_Heatmap_Top20CpGs.png            # Age-associated methylation patterns
│   ├── Figure3_ELOVL2_Correlation.png           # Top biomarker validation
│   ├── Figure4_Top20Genes.png                   # Feature importance ranking
│   └── Figure5_PathwayEnrichment.pdf            # Biological pathway analysis
│
├── scripts/                                     # Python preprocessing pipeline (HPC)
│   ├── README.md                                # Script documentation
│   ├── download_data_files.py                   # GEO data acquisition
│   ├── organized_metadata.py                    # Metadata extraction & organization
│   ├── script1_extract_combine.sbatch           # Combine beta value matrices
│   ├── script2_rename_transpose.sbatch          # Sample ID mapping & transpose
│   └── script3_combine_metadata.sbatch          # Merge methylation + phenotype data
│
└── analysis/                                    # R analysis workflow
    ├── README.md                                # Analysis documentation
    └── analysis_code.Rmd                        # Complete R analysis notebook
        # Contains:
        # - Feature selection (correlation analysis)
        # - Elastic Net model training & cross-validation
        # - SHAP interpretability analysis
        # - Pathway enrichment (KEGG, GO)
        # - Statistical analysis & visualization
```

---

## Technical Implementation

### Python Preprocessing Pipeline (HPC)

**Computational Requirements:**
- Memory: 64-128 GB RAM (chunked processing)
- Storage: ~50 GB for intermediate files
- Time: ~2-3 hours total processing

**Key Scripts:**

**1. `download_data_files.py`**
```python
import GEOparse
import pandas as pd

# Download GSE87571 from GEO
gse = GEOparse.get_GEO(geo="GSE87571", destdir="./data")

# Extract methylation matrices
beta_matrix_1 = gse.phenotype_data['GSM*']
beta_matrix_2 = gse.phenotype_data['GSM*.1']  # Duplicate samples
```

**2. `script1_extract_combine.sbatch`**
```bash
#!/bin/bash
#SBATCH --job-name=combine_matrices
#SBATCH --time=02:00:00
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=8

# Chunked processing to manage memory
python3 << EOF
import pandas as pd
chunk_size = 50000

for chunk in pd.read_csv('beta1.csv', chunksize=chunk_size):
    # Process and merge with beta2
    combined = merge_remove_duplicates(chunk)
    combined.to_csv('combined.csv', mode='a')
EOF
```

---

### R Analysis Workflow

**Environment Setup:**
```r
# Required packages
required_packages <- c(
  "glmnet",        # Elastic Net regression
  "shapviz",       # SHAP analysis
  "ggplot2",       # Visualization
  "dplyr",         # Data manipulation
  "EnrichR",       # Pathway enrichment
  "pheatmap",      # Heatmap generation
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"  # Annotation
)

# Install if needed
install.packages(required_packages)
```

**Key Analysis Steps:**

**Feature Selection:**
```r
# Pearson correlation for all CpG sites
cors <- apply(methylation_data, 2, function(cpg) {
  cor.test(cpg, age)$estimate
})

p_vals <- apply(methylation_data, 2, function(cpg) {
  cor.test(cpg, age)$p.value
})

# Filter and rank
significant <- names(p_vals[p_vals < 1e-5])
top_500 <- names(sort(abs(cors[significant]), decreasing=TRUE)[1:500])
```

**Model Training:**
```r
# Cross-validation
cv_fit <- cv.glmnet(
  x = as.matrix(train_data),
  y = train_age,
  alpha = 0.5,
  nfolds = 10,
  type.measure = "mse",
  standardize = TRUE
)

# Final model
final_model <- glmnet(
  x = as.matrix(train_data),
  y = train_age,
  alpha = 0.5,
  lambda = cv_fit$lambda.min
)

# Extract coefficients
coefficients <- coef(final_model)
selected_cpgs <- rownames(coefficients)[coefficients != 0]
```

**SHAP Analysis:**
```r
library(shapviz)

# Calculate SHAP values
shap_obj <- shapviz(final_model, X_pred = test_data)

# Summary plot (top 20 features)
sv_importance(shap_obj, kind = "beeswarm", max_display = 20)

# Dependence plot for ELOVL2
sv_dependence(shap_obj, v = "cg16867657")
```

**Pathway Enrichment:**
```r
library(enrichR)

# Get gene list
gene_list <- unique(gene_annotations$UCSC_RefGene_Name)

# Enrichment analysis
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
enriched <- enrichr(gene_list, dbs)

# Extract significant pathways
kegg_results <- enriched$KEGG_2021_Human
kegg_sig <- kegg_results[kegg_results$Adjusted.P.value < 0.05, ]
```

---

## Biological Insights

### Coordinated Epigenetic Remodeling

**Observation:**
Near-balanced methylation changes (52.9% hyper vs. 47.1% hypo)

**Interpretation:**
- NOT random methylation drift
- Suggests **programmed epigenetic aging**
- Coordinated regulation by transcription factors
- Potential master regulators controlling multiple CpG sites

**Implications:**
- Aging may be partially reversible through targeted interventions
- Key regulatory nodes could be therapeutic targets

---

### Developmental Origins of Aging

**Finding:**
Strong enrichment of neurodevelopmental pathways (p = 2.6×10⁻⁵)

**Theory:**
**Antagonistic Pleiotropy**
- Genes beneficial early in life (development) may become detrimental later (aging)
- Early developmental programs leave lasting epigenetic marks
- These marks gradually drift or accumulate changes

**Evidence from This Study:**
- FOXG1: Critical for forebrain development, methylation changes with age
- SOX2: Neural stem cell marker, shows age-related hypermethylation
- Genes active during CNS development disproportionately affected

**Clinical Relevance:**
- Understanding developmental programming could inform early-life interventions
- Critical periods for establishing healthy epigenetic trajectories

---

### Senescence as Central Aging Mechanism

**Finding:**
Cellular senescence pathway enrichment (p = 0.010)

**Senescent Cell Characteristics:**
1. Permanent cell cycle arrest
2. Resistance to apoptosis
3. Senescence-Associated Secretory Phenotype (SASP)
   - Pro-inflammatory cytokines (IL-6, IL-8)
   - Matrix metalloproteinases
   - Growth factors

**Consequences:**
- Chronic low-grade inflammation ("inflammaging")
- Tissue dysfunction
- Stem cell exhaustion
- Increased disease risk

**Therapeutic Potential:**
- **Senolytics:** Drugs that selectively eliminate senescent cells
- CpG biomarkers in senescence genes could:
  - Monitor senescent cell burden
  - Track senolytic therapy effectiveness
  - Predict treatment response

**Genes Identified:**
- TP73 (p53 family tumor suppressor)
- CDKN2A (encodes p16^INK4a senescence inducer)
- Components of inflammatory signaling

---

### Declining Regenerative Capacity

**Finding:**
Hippo signaling pathway enrichment (p = 0.002)

**Hippo Pathway Function:**
- Controls organ size during development
- Regulates stem cell proliferation vs. differentiation
- Maintains tissue homeostasis

**Age-Related Dysregulation:**
- Reduced YAP/TAZ nuclear localization
- Impaired tissue repair after injury
- Diminished organ regeneration capacity

**Connection to Aging Hallmarks:**
- **Stem cell exhaustion:** Reduced proliferative capacity
- **Altered intercellular communication:** Dysregulated growth signals
- **Loss of proteostasis:** Impaired stress responses

**Potential Interventions:**
- YAP/TAZ activators to boost regeneration
- Tissue engineering informed by developmental signals
- Targeted therapies to restore homeostatic balance

---

## Significance and Impact

### Methodological Contributions

**1. Bridges Accuracy-Interpretability Gap**
- Demonstrates that sophisticated explainable AI (SHAP) can validate simpler models
- Eliminates false choice between performance and understanding
- Provides framework applicable to other biological prediction tasks

**2. Establishes Field-Standard Workflow**
- Combines best practices from landmark studies
- Reproducible, well-documented pipeline
- Template for future epigenetic clock development

**3. Validates Biological Relevance**
- Goes beyond statistical prediction
- Connects CpG sites to mechanistic pathways
- Provides targets for experimental validation

---

### Clinical Translation Potential

**1. Precision Medicine Applications**

**Biological Age Assessment:**
- More accurate than chronological age (MAE = 2.55 years)
- Individualized health risk stratification
- Early detection of accelerated aging

**Examples:**
```
Patient A: Chronological age 45, Biological age 38
→ "Younger" phenotype, lower disease risk

Patient B: Chronological age 45, Biological age 53
→ Accelerated aging, higher intervention priority
```

**2. Therapeutic Monitoring**

**Anti-Aging Interventions:**
- Caloric restriction
- Exercise programs
- Senolytic drugs
- NAD+ boosters
- Metformin, rapamycin

**Monitoring Strategy:**
- Baseline epigenetic age measurement
- Serial assessments during intervention
- Track rate of biological aging (Δ biological age / Δ chronological age)

**Success Criteria:**
- Decelerated aging: Δ bio age < Δ chrono age
- Rejuvenation: Actual decrease in biological age

**3. Disease Risk Prediction**

**Age Acceleration as Risk Factor:**
```
Epigenetic Age Acceleration = Predicted Age - Chronological Age

Positive acceleration (older biological age) predicts:
- Cardiovascular disease
- Type 2 diabetes
- Alzheimer's disease
- Cancer
- All-cause mortality
```

**Preventive Medicine:**
- Identify high-risk individuals before disease onset
- Implement targeted preventive strategies
- Monitor effectiveness of prevention

---

### Research Impact

**1. Mechanistic Hypothesis Generation**

**Testable Questions:**
- Does ELOVL2 knockdown accelerate cellular aging?
- Can modulating Hippo signaling rejuvenate aged tissues?
- Do senolytic drugs reverse age-related methylation changes?

**2. Biomarker Validation**

**153 CpG Sites:**
- Candidates for focused methylation panels
- Targets for CRISPR-based epigenetic editing
- Readouts for aging intervention trials

**3. Multi-Omics Integration**

**Future Directions:**
- Combine with transcriptomics (gene expression)
- Integrate with proteomics (protein abundance)
- Add metabolomics (metabolic state)
- Create systems-level aging models

---

## Future Directions

### 1. External Validation

**Cross-Population Studies:**
- **African ancestry populations**
- **Asian ancestry populations**
- **Hispanic/Latino populations**
- **Admixed populations**

**Goal:** Ensure clock generalizes across genetic backgrounds

**Expected Challenges:**
- Population-specific methylation patterns
- Environmental influences
- Gene-environment interactions

**Solution:** Develop population-specific corrections or universal clocks

---

### 2. Longitudinal Tracking

**Current Limitation:**
Cross-sectional design prevents causal inference

**Proposed Study:**
- **N = 1000 participants**
- **Baseline:** Age 40-60 years
- **Follow-up:** Every 2 years for 10 years
- **Measurements:** Methylation, health outcomes, lifestyle

**Objectives:**
- Measure pace of aging (rate of biological age change)
- Identify factors that accelerate/decelerate aging
- Test whether epigenetic age predicts incident disease

---

### 3. Multi-Omics Integration

**Proposed Integrative Model:**
```
DNA Methylation (Epigenome)
    ↓ regulates
Gene Expression (Transcriptome)
    ↓ translates to
Protein Abundance (Proteome)
    ↓ drives
Metabolic State (Metabolome)
    ↓ determines
Phenotype (Clinical Outcomes)
```

**Analytical Approach:**
- Multi-omic factor analysis
- Network-based integration
- Causal inference methods

**Expected Insights:**
- Identify master regulators coordinating multi-omic aging
- Prioritize intervention targets based on centrality
- Predict response to interventions

---

### 4. Causal Validation Through Epigenetic Editing

**CRISPR-Based Approaches:**

**dCas9-DNMT3A (Methylation Writer):**
```
Target specific CpG sites → Induce hypermethylation → Measure cellular phenotypes
```

**dCas9-TET1 (Methylation Eraser):**
```
Target specific CpG sites → Induce hypomethylation → Measure cellular phenotypes
```

**Experimental Design:**
1. **Target:** cg16867657 (ELOVL2)
2. **Cell Type:** Primary human fibroblasts
3. **Intervention:** Induce premature hypermethylation
4. **Readouts:**
   - Cellular senescence markers (SA-β-gal, p16)
   - Proliferation rate
   - SASP factor secretion
   - Gene expression changes

**Key Questions:**
- Does forced ELOVL2 methylation cause aging phenotypes?
- Can reversing age-related methylation rejuvenate cells?
- Which CpG sites are drivers vs. passengers?

---

### 5. Clinical Trial Implementation

**Proposed Trial: "EpiAge-Monitor"**

**Design:**
- **Phase:** II clinical trial
- **Intervention:** Senolytic drug (Dasatinib + Quercetin)
- **Duration:** 12 months
- **N:** 200 adults (age 60-75) with metabolic syndrome

**Arms:**
1. **Treatment:** D+Q (monthly pulse dosing)
2. **Placebo:** Matched control

**Primary Endpoint:**
- Change in epigenetic age (Δ biological age)

**Secondary Endpoints:**
- Senescent cell burden (p16/p21 expression in blood)
- Inflammatory markers (IL-6, CRP)
- Physical function (6-minute walk test)
- Cognitive function (MoCA score)

**Methylation Measurements:**
- Baseline
- Month 3
- Month 6
- Month 12

**Success Criteria:**
- Treatment group shows decelerated or reversed epigenetic aging
- Correlation between Δ biological age and clinical improvements

---

### 6. Tissue-Specific Clock Development

**Rationale:**
- Blood captures systemic aging
- Specific tissues may have unique aging signatures

**Proposed Tissue Clocks:**

**Brain:**
- Neurodegenerative disease risk
- Cognitive decline prediction
- Alzheimer's/Parkinson's biomarkers

**Liver:**
- Metabolic dysfunction
- Drug metabolism capacity
- Regenerative potential

**Skin:**
- Visible aging assessment
- UV damage accumulation
- Cosmetic intervention monitoring

**Challenges:**
- Tissue biopsy invasiveness
- Cross-tissue validation
- Integration of multi-tissue clocks

---

### 7. Artificial Intelligence Advancement

**Deep Learning Approaches:**

**Convolutional Neural Networks:**
- Learn higher-order CpG interaction patterns
- Capture non-linear age relationships
- Potential for improved accuracy

**Attention Mechanisms:**
- Identify key CpG sites dynamically
- Context-dependent feature importance
- Enhanced interpretability

**Transfer Learning:**
- Pre-train on large datasets (thousands of samples)
- Fine-tune for specific populations/diseases
- Improve performance with limited data

**Explainable AI Methods:**
- Extend beyond SHAP
- Integrated Gradients
- Layer-wise Relevance Propagation
- Counterfactual explanations

**Goal:**
- Achieve MAE < 2.0 years
- Maintain biological interpretability
- Enable personalized aging predictions

---

## References

### Core Methodology Papers

1. **Horvath, S.** (2013). DNA methylation age of human tissues and cell types. *Genome Biology*, 14(10), R115. https://doi.org/10.1186/gb-2013-14-10-r115

2. **Hannum, G., et al.** (2013). Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular Cell*, 49(2), 359-367. https://doi.org/10.1016/j.molcel.2012.10.016

3. **Johansson, Å., Enroth, S., & Gyllensten, U.** (2013). Continuous aging of the human DNA methylome throughout the human lifespan. *PLoS ONE*, 8(6), e67378. https://doi.org/10.1371/journal.pone.0067378

### Interpretability & Machine Learning

4. **Lundberg, S.M., & Lee, S.I.** (2017). A unified approach to interpreting model predictions. *Advances in Neural Information Processing Systems*, 30, 4765-4774.

5. **Zou, H., & Hastie, T.** (2005). Regularization and variable selection via the Elastic Net. *Journal of the Royal Statistical Society Series B*, 67(2), 301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x

6. **Shireby, G.L., et al.** (2020). Recalibrating the epigenetic clock: Implications for assessing biological age in the human cortex. *Brain*, 143(12), 3763-3775. https://doi.org/10.1093/brain/awaa334

### Advanced Epigenetic Clocks

7. **McCrory, C., et al.** (2021). GrimAge outperforms other epigenetic clocks in the prediction of age-related clinical phenotypes and all-cause mortality. *The Journals of Gerontology: Series A*, 76(5), 741-749. https://doi.org/10.1093/gerona/glaa286

8. **Belsky, D.W., et al.** (2022). DunedinPACE, a DNA methylation biomarker of the pace of aging. *eLife*, 11, e73420. https://doi.org/10.7554/eLife.73420

### ELOVL2 as Aging Biomarker

9. **Garagnani, P., et al.** (2012). Methylation of ELOVL2 gene as a new epigenetic marker of age. *Aging Cell*, 11(6), 1132-1134. https://doi.org/10.1111/acel.12005

10. **El-Shishtawy, N.M., et al.** (2024). DNA methylation of ELOVL2 gene as an epigenetic marker of age among Egyptian population. *Egyptian Journal of Medical Human Genetics*, 25(1), 1-8. https://doi.org/10.1186/s43042-024-00477-7

### Comprehensive Reviews

11. **Teschendorff, A.E., & Horvath, S.** (2025). Epigenetic ageing clocks: Statistical methods and emerging computational challenges. *Nature Reviews Genetics*, 26(5), 350-368. https://doi.org/10.1038/s41576-024-00807-w

12. **López-Otín, C., et al.** (2023). Hallmarks of aging: An expanding universe. *Cell*, 186(2), 243-278. https://doi.org/10.1016/j.cell.2022.11.001

### Additional Supporting Literature

- **Pathway Enrichment:** Xie, Z., et al. (2021). Gene set knowledge discovery with Enrichr. *Current Protocols*, 1(3), e90.
- **Causality in Aging:** Ying, K., et al. (2024). Causality-enriched epigenetic age uncouples damage and adaptation. *Nature Aging*, 4(2), 231-246.
- **Sex Differences:** Yusipov, I., et al. (2020). Age-related DNA methylation changes are sex-specific: A comprehensive assessment. *Aging*, 12(24), 24057-24080.
- **Cellular Senescence:** Ajoolabady, A., et al. (2025). Hallmarks and mechanisms of cellular senescence in aging and disease. *Cell Death Discovery*, 11(1), 1-15.

---

## Citation

**Academic Citation:**
```
Singh, N. (2025). Development and Interpretation of an Epigenetic Clock 
Using DNA Methylation Data from GSE87571. Master's Capstone Project, 
MS in Bioinformatics, College of Science, Northeastern University. 
GitHub: https://github.com/humangene/BINF7700_Capstone_NehaSingh
```

**BibTeX:**
```bibtex
@mastersthesis{Singh2025EpigeneticClock,
  author = {Singh, Neha},
  title = {Development and Interpretation of an Epigenetic Clock Using 
           DNA Methylation Data from GSE87571},
  school = {Northeastern University},
  year = {2025},
  type = {Master's Capstone Project},
  url = {https://github.com/humangene/BINF7700_Capstone_NehaSingh}
}
```

---

## Acknowledgments

**Academic Advisors:**
- **Dr. Nabil Atallah** - Project conceptualization, methodological guidance, and critical feedback throughout development
- **Dr. Oyeronke Ayansola, PhD** - Technical mentorship, statistical consultation, and manuscript review

**Computational Resources:**
- **Northeastern University Research Computing** - Access to Discovery HPC cluster enabling large-scale methylation data processing

**Collaborators:**
- **Yash** - Collaborative discussions on research design, grant proposal development, and methodological refinements

**Data Providers:**
- **Johansson et al.** - For making GSE87571 dataset publicly available through GEO
- **Gene Expression Omnibus (GEO)** - Open-access genomic data repository

**Open-Source Community:**
- R and Python developer communities
- Authors of glmnet, SHAP, and bioinformatics packages used in this work

---

## Contact

**Neha Singh**  
MS Bioinformatics Candidate (Class of 2025)  
College of Science  
Northeastern University  
Toronto

**GitHub:** https://github.com/humangene  
**Repository:** https://github.com/humangene/BINF7700_Capstone_NehaSingh

---

## Version History

**v1.0.0** (December 2025)
- Initial release
- Complete analysis pipeline
- Full documentation
- Published results

---

## Related Materials

**Course Information:**
- **Course:** BINF7700 - Bioinformatics Capstone
- **Institution:** Northeastern University
- **Semester:** Fall 2025

**Project Components:**
- Written Report (31 pages): `final_draft_capstone.pdf`
- Research Poster: `poster.pdf`
- Code Repository: This GitHub repository
- Presentation: [Available upon request]

---

**Last Updated:** December 10, 2025  
**Document Version:** 1.0  
**README Maintained by:** Neha Singh

---

*This project represents the culmination of graduate training in bioinformatics, integrating computational methods, statistical analysis, and biological interpretation to address a significant challenge in aging research. It demonstrates mastery of data science, machine learning, and domain expertise in epigenetics.*
