# Analysis Scripts

Python and SBATCH scripts for processing GSE87571 methylation data on Northeastern's Discovery HPC cluster.

---

## Data Acquisition & Preprocessing

### `download_data_files.py`
Downloads GSE87571 dataset from Gene Expression Omnibus (GEO).
- **Input:** GEO accession GSE87571
- **Output:** Raw methylation data files
- **Dependencies:** GEOparse, pandas

### `organized_metadata.py`
Organizes sample metadata from SOFT file.
- **Input:** GSE87571 SOFT file
- **Output:** Structured metadata with sample IDs, age, sex
- **Purpose:** Maps non-descriptive sample IDs to GSM identifiers

---

## HPC Processing Pipeline (SBATCH)

All SBATCH scripts run on Northeastern University's Discovery HPC cluster using Slurm job scheduler.

### `script1_extract_combine.sbatch`
Extracts and combines methylation beta value matrices.
- **Platform:** Northeastern Discovery HPC
- **Input:** Two raw beta value matrices (729 samples each)
- **Process:** 
  - Removes duplicate samples (suffix .1)
  - Chunked processing (50,000 CpGs per chunk) for memory efficiency
  - Combines into single matrix
- **Output:** Combined dataset (729 samples × 485,512 CpG sites)
- **Resources:** 128GB RAM, 8 CPUs, ~2 hours runtime

### `script2_rename_transpose.sbatch`
Renames sample IDs and transposes data matrix.
- **Input:** Combined beta values
- **Process:**
  - Maps generic IDs to GSM identifiers
  - Transposes matrix (CpGs as columns, samples as rows)
- **Output:** 729 rows (samples) × 485,512 columns (CpGs)
- **Purpose:** Prepares data format for downstream statistical analysis

### `script3_combine_metadata.sbatch`
Merges methylation data with sample metadata.
- **Input:** 
  - Transposed beta values
  - Organized metadata (age, sex)
- **Process:** Inner join on sample IDs
- **Output:** Final dataset ready for feature selection (729 complete samples)
- **Quality Control:** Ensures all samples have complete age/sex information

---

## Workflow Summary
```
1. download_data_files.py          → Download from GEO
2. organized_metadata.py           → Extract metadata
3. script1_extract_combine.sbatch  → Combine matrices
4. script2_rename_transpose.sbatch → ID mapping + transpose
5. script3_combine_metadata.sbatch → Add phenotype data
   ↓
Ready for R analysis (feature selection, modeling)
```

---

## Requirements

**Python:**
- Python 3.13.5
- pandas >= 2.0.0
- numpy >= 1.24.0
- GEOparse >= 2.0.3

**HPC Environment:**
- Slurm job scheduler
- 128GB RAM (for large matrix operations)
- Multi-core processing capability

---

## Usage

**On Northeastern Discovery HPC:**
```bash
# Submit jobs sequentially (each depends on previous)
sbatch script1_extract_combine.sbatch
# Wait for job completion, then:
sbatch script2_rename_transpose.sbatch
# Wait for job completion, then:
sbatch script3_combine_metadata.sbatch
```

**Local execution (Python scripts):**
```bash
python download_data_files.py
python organized_metadata.py
```

---

## Output Files

After pipeline completion:
- `combined_methylation.csv` - Full methylation matrix
- `transposed_methylation.csv` - Analysis-ready format
- `final_data_with_metadata.csv` - Complete dataset for R analysis

---

## Notes

- SBATCH scripts contain hard-coded file paths specific to Discovery cluster
- Modify paths if running in different environment
- Chunked processing prevents memory overflow with large arrays
- All intermediate files saved for debugging/verification
```
