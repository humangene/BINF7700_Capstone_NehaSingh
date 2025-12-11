#!/bin/bash

#SBATCH --job-name=script2_rename

#SBATCH --output=../logs/script2_rename_%j.out

#SBATCH --error=../logs/script2_rename_%j.err

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

import re



print("="*80)

print("SCRIPT 2: RENAME X → GSM & TRANSPOSE")

print("="*80)



# ============================================================================

# CONFIGURATION

# ============================================================================



INPUT_DIR = '../data_files/chunks'

OUTPUT_DIR = '../data_files/chunks_transposed'

SOFT_FILE = '../data_files/GSE87571_family.soft'

NUM_CHUNKS = 10



# Create output directory

os.makedirs(OUTPUT_DIR, exist_ok=True)



print(f"\nConfiguration:")

print(f"  Input directory: {INPUT_DIR}")

print(f"  Output directory: {OUTPUT_DIR}")

print(f"  SOFT file: {SOFT_FILE}")

print(f"  Number of chunks: {NUM_CHUNKS}")



# ============================================================================

# STEP 1: READ ACTUAL X → GSM MAPPING FROM SOFT FILE

# ============================================================================



print("\n" + "="*80)

print("READING X → GSM MAPPING FROM SOFT FILE")

print("="*80)



x_to_gsm = {}

current_gsm = None



print(f"\nReading {SOFT_FILE}...")



with open(SOFT_FILE, 'r', encoding='utf-8', errors='ignore') as f:

    for line in f:

        line = line.strip()

        

        # Find GSM ID

        if line.startswith('^SAMPLE'):

            match = re.search(r'GSM\d+', line)

            if match:

                current_gsm = match.group(0)

        

        # Find X number in title

        elif line.startswith('!Sample_title') and current_gsm:

            # Title format: "X1 genomic DNA from whole blood"

            match = re.search(r'X(\d+)\s+genomic', line)

            if match:

                x_num = int(match.group(1))

                x_to_gsm[f"X{x_num}"] = current_gsm



print(f"✓ Found {len(x_to_gsm)} X → GSM mappings")



# Verify mapping

print(f"\n--- MAPPING VERIFICATION (first 10) ---")

sorted_keys = sorted(x_to_gsm.keys(), key=lambda x: int(x.replace('X', '')))

for key in sorted_keys[:10]:

    print(f"  {key} → {x_to_gsm[key]}")



print(f"\n--- MAPPING VERIFICATION (last 10) ---")

for key in sorted_keys[-10:]:

    print(f"  {key} → {x_to_gsm[key]}")



# ============================================================================

# PROCESS EACH CHUNK

# ============================================================================



print("\n" + "="*80)

print("PROCESSING CHUNKS")

print("="*80)



for chunk_num in range(NUM_CHUNKS):

    

    input_file = f"{INPUT_DIR}/chunk_{chunk_num}.csv.gz"

    output_file = f"{OUTPUT_DIR}/chunk_{chunk_num}_transposed.csv.gz"

    

    # Check if input file exists

    if not os.path.exists(input_file):

        print(f"\n chunk {chunk_num}: Input file not found, skipping...")

        continue

    

    print(f"\n--- CHUNK {chunk_num} ---")

    

    # ========================================================================

    # STEP 2: Load chunk

    # ========================================================================

    

    print(f"  Loading {input_file}...")

    df = pd.read_csv(input_file, compression='gzip')

    print(f"    Loaded: {df.shape[0]} CpGs × {df.shape[1]} columns")

    

    # ========================================================================

    # STEP 3: Rename X columns to GSM IDs using SOFT file mapping

    # ========================================================================

    

    print(f"  Renaming X → GSM using SOFT file mapping...")

    

    # Get current columns

    current_cols = df.columns.tolist()

    

    # Create rename dictionary for this chunk

    rename_dict = {}

    renamed_count = 0

    not_found = []

    

    for col in current_cols:

        if col in x_to_gsm:

            rename_dict[col] = x_to_gsm[col]

            renamed_count += 1

        elif col != 'ID_REF' and col.startswith('X'):

            not_found.append(col)

    

    # Rename columns

    df = df.rename(columns=rename_dict)

    

    print(f"    Renamed {renamed_count} columns")

    if len(not_found) > 0:

        print(f"    Not found in mapping: {len(not_found)} columns")

        print(f"       Examples: {not_found[:5]}")

    

    print(f"    First 5 columns now: {df.columns.tolist()[:5]}")

    

    # ========================================================================

    # STEP 4: Set ID_REF as index (CpG IDs become row names)

    # ========================================================================

    

    print(f"  Setting ID_REF as index...")

    df = df.set_index('ID_REF')

    print(f"    Shape after setting index: {df.shape}")

    

    # ========================================================================

    # STEP 5: Transpose (CpGs become columns, Samples become rows)

    # ========================================================================

    

    print(f"  Transposing...")

    df_transposed = df.T

    print(f"    Before transpose: {df.shape[0]} CpGs × {df.shape[1]} samples")

    print(f"    After transpose: {df_transposed.shape[0]} samples × {df_transposed.shape[1]} CpGs")

    

    # Free memory

    del df

    gc.collect()

    

    # ========================================================================

    # STEP 6: Reset index to make SampleID a column

    # ========================================================================

    

    print(f"  Resetting index...")

    df_transposed = df_transposed.reset_index()

    df_transposed = df_transposed.rename(columns={'index': 'SampleID'})

    print(f"    First column: {df_transposed.columns[0]}")

    print(f"    First 5 SampleIDs: {df_transposed['SampleID'].head().tolist()}")

    

    # ========================================================================

    # STEP 7: Save transposed chunk

    # ========================================================================

    

    print(f"  Saving to {output_file}...")

    df_transposed.to_csv(output_file, index=False, compression='gzip')

    

    file_size = os.path.getsize(output_file) / (1024**2)

    print(f"    Saved! File size: {file_size:.1f} MB")

    

    # Free memory

    del df_transposed

    gc.collect()

    

    print(f"  Chunk {chunk_num} complete!")



# ============================================================================

# SUMMARY

# ============================================================================



print("\n" + "="*80)

print("SCRIPT 2 COMPLETE!")

print("="*80)



# List output files

print(f"\nOutput files in {OUTPUT_DIR}:")

for f in sorted(os.listdir(OUTPUT_DIR)):

    if f.endswith('.csv.gz'):

        fpath = f"{OUTPUT_DIR}/{f}"

        fsize = os.path.getsize(fpath) / (1024**2)

        print(f"  {f}: {fsize:.1f} MB")



print("\n Ready for Script 3: Combine Chunks + Add Metadata")



EOF



echo ""

echo "Script 2 complete: $(date)"
