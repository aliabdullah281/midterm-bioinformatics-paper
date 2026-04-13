"""
Step 2.7 — Inspect and extract diabetes-relevant data from biokg-source.tar.gz.
Extracts the archive, scans for diabetes-relevant content, saves any useful subset.
Output: output/extracted/biokg_diabetes.csv (if relevant data found)
"""
import sys, tarfile, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    BIOKG_TAR, BIOKG_CONTENTS_DIR, BIOKG_CSV,
    EXTRACTED_DIR, DIABETES_SEED_GENES, DIABETES_UMLS_IDS,
)


def extract():
    print("\n[Step 2.7] Extracting and scanning biokg-source.tar.gz")
    print(f"  Source: {BIOKG_TAR}")

    if not BIOKG_TAR.exists():
        print(f"  ERROR: {BIOKG_TAR} not found. Skipping.")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    BIOKG_CONTENTS_DIR.mkdir(parents=True, exist_ok=True)

    # Extract archive
    print(f"  Extracting to {BIOKG_CONTENTS_DIR} ...")
    with tarfile.open(BIOKG_TAR, 'r:gz') as tar:
        members = tar.getmembers()
        print(f"  Archive contains {len(members):,} files:")
        for m in members[:50]:
            print(f"    {m.name}  ({m.size:,} bytes)")
        if len(members) > 50:
            print(f"    ... and {len(members)-50} more")
        tar.extractall(path=BIOKG_CONTENTS_DIR)

    # Scan extracted files for diabetes-relevant content
    diabetes_keywords = (
        {g.lower() for g in DIABETES_SEED_GENES} |
        DIABETES_UMLS_IDS |
        {'diabetes', 'insulin', 'glucose', 'pancreas', 'beta cell'}
    )

    all_relevant = []
    for fpath in BIOKG_CONTENTS_DIR.rglob('*'):
        if not fpath.is_file():
            continue
        suffix = fpath.suffix.lower()
        if suffix not in {'.tsv', '.csv', '.txt', '.tab'}:
            continue

        try:
            sep = '\t' if suffix in {'.tsv', '.tab', '.txt'} else ','
            df = pd.read_csv(fpath, sep=sep, nrows=5, header='infer', low_memory=False)
            print(f"  Scanning {fpath.name}: columns = {list(df.columns)[:6]}")

            # Check if any column contains diabetes keywords
            found = False
            for col in df.columns:
                sample_vals = df[col].dropna().astype(str).str.lower().tolist()
                if any(kw in val for kw in diabetes_keywords for val in sample_vals):
                    found = True
                    break

            if found:
                # Read full file and filter
                df_full = pd.read_csv(fpath, sep=sep, low_memory=False)
                # Filter rows that mention any keyword anywhere
                mask = df_full.apply(
                    lambda col: col.astype(str).str.lower().str.contains(
                        '|'.join(['diabetes', 'insulin', 'glucose'] + list(DIABETES_SEED_GENES)),
                        na=False
                    )
                ).any(axis=1)
                df_filt = df_full[mask].copy()
                df_filt['source_file'] = fpath.name
                print(f"    DIABETES DATA FOUND: {len(df_filt):,} rows from {fpath.name}")
                all_relevant.append(df_filt)
        except Exception as e:
            print(f"  Could not read {fpath.name}: {e}")

    if all_relevant:
        combined = pd.concat(all_relevant, ignore_index=True)
        combined.to_csv(BIOKG_CSV, index=False)
        print(f"  Saved {len(combined):,} bioKG diabetes rows → {BIOKG_CSV}")
        return combined
    else:
        print("  No diabetes-relevant data found in biokg-source.tar.gz — skipping.")
        return None


if __name__ == "__main__":
    extract()
