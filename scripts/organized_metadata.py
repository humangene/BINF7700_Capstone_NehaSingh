#!/usr/bin/env python3

# imported the libraries that are needed

import GEOparse

import pandas as pd

from pathlib import Path



# Paths for the input and output files are declared

data_dir = Path(__file__).resolve().parent.parent / "data_files"

soft_file = data_dir / "GSE87571_family.soft"

output_csv = data_dir / "GSE87571_metadata.csv"



print(f"Loading SOFT file: {soft_file}")

gse = GEOparse.get_GEO(filepath=str(soft_file))



#  Extract phenotype metadata directly from DataFrame

metadata = gse.phenotype_data.copy()

print(f"\nTotal samples found: {metadata.shape[0]}")



# Keep only SampleID, Age, Sex

columns_needed = ['geo_accession', 'characteristics_ch1.0.gender', 'characteristics_ch1.1.age']

for col in columns_needed:

    if col not in metadata.columns:

        raise ValueError(f"Column not found in phenotype_data!")



metadata = metadata[columns_needed]


# Rename columns

metadata.rename(columns={

    'geo_accession': 'SampleID',

    'characteristics_ch1.0.gender': 'Sex',

    'characteristics_ch1.1.age': 'Age'

}, inplace=True)



# Clean the columns and Remove prefix like "gender: " or "age: "

metadata['Sex'] = metadata['Sex'].astype(str).str.split(':').str[-1].str.strip().str.capitalize()

metadata['Age'] = metadata['Age'].astype(str).str.split(':').str[-1].str.strip()


# Convert Age to numeric

metadata['Age'] = pd.to_numeric(metadata['Age'], errors='coerce')


# Drop rows with missing Age or Sex

metadata.dropna(subset=['Age', 'Sex'], inplace=True)



# Save CSV

metadata.to_csv(output_csv, index=False)

print(f"\n Organized metadata saved as: {output_csv}")

print(f"Final metadata shape: {metadata.shape}")



# Summary

print("\n Summary statistics:")

print(f"Age range: {metadata['Age'].min()} - {metadata['Age'].max()} years")

print(metadata['Sex'].value_counts())

print("\n Preview:")

print(metadata.head())


