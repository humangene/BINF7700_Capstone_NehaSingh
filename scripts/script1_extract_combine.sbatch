#!/bin/bash

#SBATCH --job-name=script1_extract

#SBATCH --output=../logs/script1_extract_%j.out

#SBATCH --error=../logs/script1_extract_%j.err

#SBATCH --partition=short

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem=64G

#SBATCH --time=06:00:00



module load python/3.13.5

source ~/myenv/bin/activate



python3 << 'EOF'

import pandas as pd

import numpy as np

import os

import gc



print("="*80)

print("SCRIPT 1: EXTRACT & COMBINE BOTH MATRICES")

print("="*80)



# ============================================================================

# CONFIGURATION

# ============================================================================



CHUNK_SIZE = 50000

TOTAL_CPGS = 485512

NUM_CHUNKS = (TOTAL_CPGS // CHUNK_SIZE) + 1  # 10 chunks



MATRIX1_PATH = '../data_files/GSE87571_matrix1of2.txt.gz'

MATRIX2_PATH = '../data_files/GSE87571%5Fmatrix2of2.txt.gz'



# Check if Matrix 2 path exists, try alternative

if not os.path.exists(MATRIX2_PATH):

    MATRIX2_PATH = '../data_files/GSE87571_matrix2of2.txt.gz'



OUTPUT_DIR = '../data_files/chunks'



# Create output directory

os.makedirs(OUTPUT_DIR, exist_ok=True)



print(f"\nConfiguration:")

print(f"  Chunk size: {CHUNK_SIZE:,} CpGs")

print(f"  Total CpGs: {TOTAL_CPGS:,}")

print(f"  Number of chunks: {NUM_CHUNKS}")

print(f"  Matrix 1: {MATRIX1_PATH}")

print(f"  Matrix 2: {MATRIX2_PATH}")

print(f"  Output directory: {OUTPUT_DIR}")



# ============================================================================

# FUNCTION: Remove .1 columns

# ============================================================================



def remove_replicate_columns(df):

    """Keep only non-.1 columns (remove bad replicates)"""

    cols_to_keep = ['ID_REF']

    for col in df.columns:

        if col != 'ID_REF' and not col.endswith('.1'):

            cols_to_keep.append(col)

    return df[cols_to_keep]



# ============================================================================

# PROCESS EACH CHUNK

# ============================================================================



print("\n" + "="*80)

print("PROCESSING CHUNKS")

print("="*80)



for chunk_num in range(NUM_CHUNKS):

    

    # Calculate row range for this chunk

    start_row = chunk_num * CHUNK_SIZE

    end_row = min((chunk_num + 1) * CHUNK_SIZE, TOTAL_CPGS)

    n_rows = end_row - start_row

    

    # Skip if no rows to process

    if n_rows <= 0:

        continue

    

    print(f"\n--- CHUNK {chunk_num} ---")

    print(f"  CpG range: {start_row:,} to {end_row:,} ({n_rows:,} CpGs)")

    

    # ========================================================================

    # STEP 1: Load chunk from Matrix 1

    # ========================================================================

    

    print(f"  Loading Matrix 1...")

    

    df1 = pd.read_csv(

        MATRIX1_PATH,

        compression='gzip',

        sep='\t',

        skiprows=range(1, start_row + 1) if start_row > 0 else None,

        nrows=n_rows

    )

    

    print(f"    Loaded: {df1.shape[0]} rows × {df1.shape[1]} columns")

    

    # ========================================================================

    # STEP 2: Remove .1 columns from Matrix 1

    # ========================================================================

    

    print(f"  Removing .1 columns from Matrix 1...")

    df1_clean = remove_replicate_columns(df1)

    print(f"    After cleaning: {df1_clean.shape[0]} rows × {df1_clean.shape[1]} columns")

    

    # Free memory

    del df1

    gc.collect()

    

    # ========================================================================

    # STEP 3: Load chunk from Matrix 2

    # ========================================================================

    

    print(f"  Loading Matrix 2...")

    

    df2 = pd.read_csv(

        MATRIX2_PATH,

        compression='gzip',

        sep='\t',

        skiprows=range(1, start_row + 1) if start_row > 0 else None,

        nrows=n_rows

    )

    

    print(f"    Loaded: {df2.shape[0]} rows × {df2.shape[1]} columns")

    

    # ========================================================================

    # STEP 4: Remove .1 columns from Matrix 2

    # ========================================================================

    

    print(f"  Removing .1 columns from Matrix 2...")

    df2_clean = remove_replicate_columns(df2)

    print(f"    After cleaning: {df2_clean.shape[0]} rows × {df2_clean.shape[1]} columns")

    

    # Free memory

    del df2

    gc.collect()

    

    # ========================================================================

    # STEP 5: Verify CpG IDs match

    # ========================================================================

    

    print(f"  Verifying CpG IDs match...")

    

    if not df1_clean['ID_REF'].equals(df2_clean['ID_REF']):

        print(f"    WARNING: CpG IDs don't match! Attempting to align...")

        # Merge on ID_REF to align

        merged = pd.merge(df1_clean, df2_clean, on='ID_REF', how='inner')

    else:

        print(f"    CpG IDs match perfectly!")

        # ====================================================================

        # STEP 6: Merge horizontally (combine columns)

        # ====================================================================

        

        print(f"  Merging Matrix 1 and Matrix 2...")

        

        # Drop ID_REF from df2_clean before merging (it's duplicate)

        df2_no_id = df2_clean.drop(columns=['ID_REF'])

        

        # Combine horizontally

        merged = pd.concat([df1_clean, df2_no_id], axis=1)

    

    print(f"    Merged: {merged.shape[0]} rows × {merged.shape[1]} columns")

    

    # Free memory

    del df1_clean, df2_clean

    gc.collect()

    

    # ========================================================================

    # STEP 7: Verify we have 732 samples + 1 ID column = 733 columns

    # ========================================================================

    

    expected_cols = 733  # ID_REF + 732 samples

    if merged.shape[1] != expected_cols:

        print(f"    WARNING: Expected {expected_cols} columns, got {merged.shape[1]}")

    else:

        print(f"    Correct number of columns: {merged.shape[1]}")

    

    # ========================================================================

    # STEP 8: Save chunk

    # ========================================================================

    

    output_file = f"{OUTPUT_DIR}/chunk_{chunk_num}.csv.gz"

    print(f"  Saving to {output_file}...")

    

    merged.to_csv(output_file, index=False, compression='gzip')

    

    file_size = os.path.getsize(output_file) / (1024**2)

    print(f"    Saved! File size: {file_size:.1f} MB")

    

    # Free memory

    del merged

    gc.collect()

    

    print(f"  Chunk {chunk_num} complete!")



# ============================================================================

# SUMMARY

# ============================================================================



print("\n" + "="*80)

print("SCRIPT 1 COMPLETE!")

print("="*80)



# List output files

print(f"\nOutput files in {OUTPUT_DIR}:")

for f in sorted(os.listdir(OUTPUT_DIR)):

    if f.startswith('chunk_') and f.endswith('.csv.gz'):

        fpath = f"{OUTPUT_DIR}/{f}"

        fsize = os.path.getsize(fpath) / (1024**2)

        print(f"  {f}: {fsize:.1f} MB")



print("\n Ready for Script 2: Rename X → GSM & Transpose")



EOF



echo ""

echo "Script 1 complete: $(date)"
