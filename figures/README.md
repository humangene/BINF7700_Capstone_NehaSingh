# Project Figures

Publication-quality visualizations generated from the epigenetic clock analysis.

---

## Figure Descriptions

### **Figure1_PredictedVsActual.png**
**Elastic Net Model Performance**
- Shows predicted age vs. actual chronological age on test set (n=146)
- Demonstrates strong correlation (r = 0.988, R² = 0.976)
- Mean Absolute Error: 2.55 years
- Perfect diagonal indicates accurate age prediction across lifespan (14-94 years)

### **Figure2_Heatmap_Top20CpGs.png**
**Age-Associated Methylation Patterns**
- Heatmap of top 20 CpG sites ranked by coefficient magnitude
- Rows: CpG sites (with gene annotations)
- Columns: 729 samples ordered by age
- Color scale: Blue (hypomethylation) to Red (hypermethylation)
- Shows coordinated age-related methylation changes

### **Figure3_ELOVL2_Correlation.png**
**Premier Aging Biomarker Validation**
- Scatter plot: cg16867657 (ELOVL2) methylation vs. chronological age
- Strong positive correlation (r = 0.9464)
- Color-coded by sex (demonstrates sex-independent pattern)
- Linear trend confirms ELOVL2 as robust aging biomarker

### **Figure4_Top20Genes.png**
**Feature Importance Ranking**
- Horizontal bar chart of top 20 genes by absolute Elastic Net coefficient
- Red bars: Hypermethylated with age (positive coefficients)
- Blue bars: Hypomethylated with age (negative coefficients)
- ELOVL2 shows strongest effect (β = +28.15 years)
- Balanced distribution demonstrates coordinated epigenetic remodeling

### **Figure5_PathwayEnrichment.pdf**
**Biological Pathway Analysis**
- Dot plot showing enriched KEGG pathways and GO Biological Processes
- X-axis: -log10(p-value) (significance)
- Red dashed line: p = 0.05 threshold
- Top pathways: CNS development, cellular senescence, Hippo signaling
- Connects statistical predictions to mechanistic aging biology

---

## File Formats

- **PNG files:** High-resolution (300 DPI) for publication quality
- **PDF file:** Vector graphics for scalability

## Usage

These figures are referenced throughout the project report and README documentation.
```

---

## **Then:**

**Scroll down** and in the commit message box, write:
```
Create figures folder with descriptions
