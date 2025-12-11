#!/bin/bash

#SBATCH --job-name=script3_combine

#SBATCH --output=../logs/script3_combine_%j.out

#SBATCH --error=../logs/script3_combine_%j.err

#SBATCH --partition=short

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem=100G

#SBATCH --time=06:00:00



module load python/3.13.5

source ~/myenv/bin/activate



python3 << 'EOF'

import pandas as pd

import numpy as np

import os

import gc



print("="*80)

print("SCRIPT 3: COMBINE CHUNKS + ADD METADATA")

print("="*80)



# ============================================================================

# CONFIGURATION

# ============================================================================



CHUNKS_DIR = '../data_files/chunks_transposed'

METADATA_FILE = '../data_files/GSE87571_metadata.csv'

OUTPUT_FILE = '../data_files/combined_beta_clean.csv.gz'

NUM_CHUNKS = 10



print(f"\nConfiguration:")

print(f"  Chunks directory: {CHUNKS_DIR}")

print(f"  Metadata file: {METADATA_FILE}")

print(f"  Output file: {OUTPUT_FILE}")



# ============================================================================

# SAFETY CHECK 1: VERIFY ALL CHUNKS HAVE SAME SAMPLES IN SAME ORDER

# ============================================================================



print("\n" + "="*80)

print("SAFETY CHECK 1: VERIFYING SAMPLE ORDER ACROSS ALL CHUNKS")

print("="*80)



all_sample_ids = {}



for i in range(NUM_CHUNKS):

    chunk_file = f'{CHUNKS_DIR}/chunk_{i}_transposed.csv.gz'

    df_temp = pd.read_csv(chunk_file, compression='gzip', usecols=['SampleID'])

    all_sample_ids[i] = df_temp['SampleID'].tolist()

    print(f"  Chunk {i}: {len(all_sample_ids[i])} samples")

    del df_temp

    gc.collect()



# Compare all chunks to chunk 0

reference = all_sample_ids[0]

all_match = True



for i in range(1, NUM_CHUNKS):

    if all_sample_ids[i] != reference:

        print(f"  Chunk {i} has DIFFERENT sample order!")

        all_match = False



if all_match:

    print(f"\n SAFETY CHECK 1 PASSED: All {NUM_CHUNKS} chunks have SAME {len(reference)} samples in SAME order!")

else:

    print(f"\n SAFETY CHECK 1 FAILED: Sample order differs between chunks!")

    print("   ABORTING SCRIPT!")

    exit(1)



# Clear memory

del all_sample_ids

gc.collect()



# ============================================================================

# STEP 1: LOAD AND COMBINE ALL CHUNKS

# ============================================================================



print("\n" + "="*80)

print("STEP 1: LOADING AND COMBINING ALL CHUNKS")

print("="*80)



combined_df = None



for i in range(NUM_CHUNKS):

    chunk_file = f'{CHUNKS_DIR}/chunk_{i}_transposed.csv.gz'

    print(f"\n  Loading chunk {i}...")

    

    chunk_df = pd.read_csv(chunk_file, compression='gzip')

    print(f"    Shape: {chunk_df.shape[0]} samples × {chunk_df.shape[1]} columns")

    

    if combined_df is None:

        # First chunk - use as base

        combined_df = chunk_df

        print(f"     Base chunk loaded")

    else:

        # Merge on SampleID (add CpG columns)

        # Drop SampleID from chunk_df to avoid duplicate

        cpg_cols = [col for col in chunk_df.columns if col != 'SampleID']

        chunk_df_cpgs = chunk_df[cpg_cols]

        

        # Verify row order matches before concatenating

        if not chunk_df['SampleID'].equals(combined_df['SampleID']):

            print(f"     ERROR: SampleID order mismatch in chunk {i}!")

            exit(1)

        

        # Concatenate columns (horizontal merge)

        combined_df = pd.concat([combined_df, chunk_df_cpgs], axis=1)

        print(f"     Merged. Combined shape: {combined_df.shape}")

    

    # Free memory

    del chunk_df

    gc.collect()



print(f"\n All chunks combined!")

print(f"   Final shape: {combined_df.shape[0]} samples × {combined_df.shape[1]} columns")



# ============================================================================

# SAFETY CHECK 2: VERIFY COMBINED DATA DIMENSIONS

# ============================================================================



print("\n" + "="*80)

print("SAFETY CHECK 2: VERIFYING COMBINED DATA DIMENSIONS")

print("="*80)



expected_samples = 732

expected_cpgs = 485512

expected_cols = expected_cpgs + 1  # CpGs + SampleID



actual_samples = combined_df.shape[0]

actual_cols = combined_df.shape[1]



print(f"  Expected: {expected_samples} samples × {expected_cols} columns")

print(f"  Actual: {actual_samples} samples × {actual_cols} columns")



if actual_samples == expected_samples and actual_cols == expected_cols:

    print(f"\n SAFETY CHECK 2 PASSED: Dimensions are correct!")

else:

    print(f"\n SAFETY CHECK 2 WARNING: Dimensions differ from expected!")

    print(f"   Continuing anyway...")



# ============================================================================

# STEP 2: LOAD METADATA

# ============================================================================



print("\n" + "="*80)

print("STEP 2: LOADING METADATA")

print("="*80)



metadata = pd.read_csv(METADATA_FILE)

print(f"  Metadata shape: {metadata.shape[0]} samples × {metadata.shape[1]} columns")

print(f"  Metadata columns: {metadata.columns.tolist()}")

print(f"  First 5 samples in metadata:")

print(metadata.head())



# ============================================================================

# STEP 3: MERGE WITH METADATA

# ============================================================================



print("\n" + "="*80)

print("STEP 3: MERGING WITH METADATA")

print("="*80)



print(f"\n  Before merge:")

print(f"    Combined data: {combined_df.shape[0]} samples")

print(f"    Metadata: {metadata.shape[0]} samples")



# Inner merge - keeps only samples in BOTH files

final_df = pd.merge(metadata, combined_df, on='SampleID', how='inner')



print(f"\n  After merge:")

print(f"    Final data: {final_df.shape[0]} samples × {final_df.shape[1]} columns")



# Free memory

del combined_df

gc.collect()



# ============================================================================

# SAFETY CHECK 3: VERIFY MERGE RESULTS

# ============================================================================



print("\n" + "="*80)

print("SAFETY CHECK 3: VERIFYING MERGE RESULTS")

print("="*80)



expected_final_samples = 729

actual_final_samples = final_df.shape[0]



print(f"  Expected samples after merge: {expected_final_samples}")

print(f"  Actual samples after merge: {actual_final_samples}")



if actual_final_samples == expected_final_samples:

    print(f"\n SAFETY CHECK 3 PASSED: Correct number of samples!")

else:

    print(f"\n SAFETY CHECK 3 WARNING: Sample count differs!")

    print(f"   Expected {expected_final_samples}, got {actual_final_samples}")



# Verify columns order

print(f"\n  First 10 columns: {final_df.columns.tolist()[:10]}")

print(f"  Total columns: {len(final_df.columns)}")



# ============================================================================

# SAFETY CHECK 4: SPOT CHECK SPECIFIC SAMPLES

# ============================================================================



print("\n" + "="*80)

print("SAFETY CHECK 4: SPOT CHECK - VERIFY SAMPLE DATA MATCHES")

print("="*80)



# Pick 3 random samples to verify

spot_check_samples = ['GSM2333901', 'GSM2333950', 'GSM2334500']



for sample_id in spot_check_samples:

    if sample_id in final_df['SampleID'].values:

        row = final_df[final_df['SampleID'] == sample_id].iloc[0]

        print(f"\n  Sample: {sample_id}")

        print(f"    Sex: {row['Sex']}")

        print(f"    Age: {row['Age']}")

        

        # Get first 5 CpG columns (after SampleID, Sex, Age)

        cpg_cols = [col for col in final_df.columns if col.startswith('cg')][:5]

        print(f"    First 5 CpG values:")

        for cpg in cpg_cols:

            print(f"      {cpg}: {row[cpg]}")

    else:

        print(f"\n  Sample {sample_id}: Not found in final data")



# ============================================================================

# STEP 4: SAVE FINAL CLEAN DATASET

# ============================================================================



print("\n" + "="*80)

print("STEP 4: SAVING FINAL CLEAN DATASET")

print("="*80)



print(f"\n  Saving to {OUTPUT_FILE}...")

final_df.to_csv(OUTPUT_FILE, index=False, compression='gzip')



file_size = os.path.getsize(OUTPUT_FILE) / (1024**3)

print(f"   Saved! File size: {file_size:.2f} GB")



# ============================================================================

# FINAL SUMMARY

# ============================================================================



print("\n" + "="*80)

print("SCRIPT 3 COMPLETE!")

print("="*80)



print(f"\n FINAL DATASET SUMMARY:")

print(f"   File: {OUTPUT_FILE}")

print(f"   Samples (rows): {final_df.shape[0]}")

print(f"   Columns: {final_df.shape[1]}")

print(f"   Column breakdown:")

print(f"     - SampleID: 1")

print(f"     - Sex: 1")

print(f"     - Age: 1")

print(f"     - CpGs: {final_df.shape[1] - 3}")



print(f"\n READY FOR FEATURE SELECTION!")



EOF



echo ""

echo "Script 3 complete: $(date)"
