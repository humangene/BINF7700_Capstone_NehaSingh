# R Analysis Code

Complete R analysis workflow implemented in RMarkdown notebook.

---

## Overview

This folder contains the comprehensive R analysis pipeline for developing and validating the epigenetic clock using Elastic Net regression with SHAP interpretability.

---

## File

### `analysis_code.Rmd`

**RMarkdown notebook containing the complete analysis workflow:**

1. **Data Loading & Preprocessing**
   - Import final dataset with methylation + metadata
   - Quality control checks
   - Data formatting for analysis

2. **Feature Selection**
   - Pearson correlation analysis (p < 1×10⁻⁵)
   - Ranking by absolute correlation coefficient
   - Selection of top 500 age-associated CpG sites

3. **Model Development**
   - Elastic Net regression implementation
   - 80/20 train-test split (stratified by age)
   - 10-fold cross-validation for hyperparameter optimization
   - Optimal parameters: α = 0.5, λ = 0.1498
   - Final feature selection: 153 CpG sites

4. **Model Evaluation**
   - Performance metrics: MAE, RMSE, R², Pearson correlation
   - Comparison with Horvath and Hannum clocks
   - Age-stratified analysis
   - Residual diagnostics

5. **SHAP Interpretability Analysis**
   - SHAP value calculation for all 153 features
   - Feature importance validation
   - Dependence plots for top predictors
   - Correlation between coefficients and SHAP values

6. **Gene Annotation**
   - Mapping CpG sites to genes (hg19 assembly)
   - Genomic context analysis (islands, shores, gene bodies)
   - Identification of 108 unique genes

7. **Pathway Enrichment**
   - KEGG pathway analysis
   - Gene Ontology Biological Process enrichment
   - Statistical testing (Fisher's exact, FDR correction)
   - Biological interpretation

8. **Visualization**
   - Predicted vs. actual age scatter plots
   - Heatmaps of methylation patterns
   - SHAP summary and dependence plots
   - Feature importance bar charts
   - Pathway enrichment dot plots

---

## Requirements

### R Version
- R >= 4.4.1

### Required Packages

**Core Analysis:**
```r
glmnet           # Elastic Net regression
caret            # Machine learning framework
```

**Interpretability:**
```r
shapviz          # SHAP analysis
```

**Annotation:**
```r
IlluminaHumanMethylation450kanno.ilmn12.hg19  # CpG annotations
```

**Pathway Enrichment:**
```r
enrichR          # Pathway analysis interface
```

**Data Manipulation:**
```r
dplyr            # Data wrangling
tidyr            # Data tidying
```

**Visualization:**
```r
ggplot2          # Publication-quality plots
pheatmap         # Heatmaps
RColorBrewer     # Color palettes
```

### Installation
```r
# Install required packages
install.packages(c("glmnet", "caret", "shapviz", "dplyr", 
                   "tidyr", "ggplot2", "pheatmap", "RColorBrewer"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# enrichR from CRAN
install.packages("enrichR")
```

---

## Usage

### Opening the Notebook

**In RStudio:**
1. Open `analysis_code.Rmd`
2. Click "Knit" to generate HTML report
3. Or run code chunks interactively

**In R:**
```r
rmarkdown::render("analysis_code.Rmd")
```

### Expected Outputs

The notebook generates:
- All publication figures (saved to `figures/` directory)
- Statistical results tables
- Model objects (.rds files)
- HTML report with embedded results

---

## Reproducibility

**To reproduce the complete analysis:**

1. **Ensure data availability:**
   - Final preprocessed data from Python pipeline
   - `final_data_with_metadata.csv` (729 samples × 450,282 CpGs + age/sex)

2. **Set working directory:**
```r
   setwd("path/to/analysis/")
```

3. **Run the notebook:**
   - Execute all chunks sequentially
   - Or knit entire document

4. **Expected runtime:**
   - Complete analysis: ~30-60 minutes
   - Feature selection: ~5-10 minutes
   - Model training: ~10-20 minutes
   - SHAP analysis: ~15-30 minutes

---

## Key Results

**Model Performance (Test Set, n=146):**
- MAE: 2.55 years
- RMSE: 3.21 years
- R²: 0.976
- Pearson r: 0.988 (p < 2.2×10⁻¹⁶)

**Feature Selection:**
- 153 CpG sites selected
- 81 hypermethylated (52.9%)
- 72 hypomethylated (47.1%)

**Top Predictor:**
- cg16867657 (ELOVL2)
- Coefficient: β = +28.15
- Correlation: r = 0.9464
- SHAP linearity: R² = 1.0

**Biological Pathways:**
- CNS development (p = 2.6×10⁻⁵)
- Cellular senescence (p = 0.010)
- Hippo signaling (p = 0.002)

---

## Notes

- Notebook includes detailed comments explaining each step
- All random seeds set for reproducibility
- Figures optimized for publication quality 

---

## Contact

For questions about the R analysis code, please refer to the main project README or contact the project author.
