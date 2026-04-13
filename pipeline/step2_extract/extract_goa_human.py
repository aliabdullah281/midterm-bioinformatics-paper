"""
Step 2.4 — Extract human GO annotations for diabetes proteins from GOA GAF file.
Depends on: diabetes_gene_set.txt (from extract_biogrid.py)
Output: output/extracted/goa_human_diabetes.csv
"""
import sys, gzip
from pathlib import Path
from io import StringIO
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    GOA_GAF, GOA_CSV, DIABETES_GENE_SET_FILE, GAF_COLS,
    EXTRACTED_DIR, DIABETES_SEED_GENES,
)


def extract():
    print("\n[Step 2.4] Extracting human GO annotations (diabetes proteins)")
    print(f"  Source: {GOA_GAF}")

    if not GOA_GAF.exists():
        print(f"  ERROR: {GOA_GAF} not found. Skipping.")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    # Load expanded diabetes gene set
    if DIABETES_GENE_SET_FILE.exists():
        with open(DIABETES_GENE_SET_FILE) as f:
            diabetes_genes = set(line.strip() for line in f if line.strip())
        print(f"  Loaded {len(diabetes_genes):,} genes from diabetes_gene_set.txt")
    else:
        diabetes_genes = DIABETES_SEED_GENES.copy()
        print(f"  WARNING: diabetes_gene_set.txt not found — using {len(diabetes_genes)} seed genes only")

    # Read GAF, skipping '!' comment lines
    print("  Reading GAF file (skipping comment lines) ...")
    with gzip.open(GOA_GAF, 'rt', encoding='utf-8') as f:
        data_lines = [line for line in f if not line.startswith('!')]
    print(f"  Data lines found   : {len(data_lines):,}")

    df = pd.read_csv(
        StringIO(''.join(data_lines)),
        sep='\t',
        header=None,
        names=GAF_COLS[:17],   # GAF 2.2: up to 17 columns
        low_memory=False,
        on_bad_lines='skip',
    )
    print(f"  Total annotations  : {len(df):,}")

    # Filter: human (taxon field contains 'taxon:9606')
    mask_taxon = df['Taxon'].str.contains('taxon:9606', na=False)
    df = df[mask_taxon].copy()
    print(f"  Human rows         : {len(df):,}")

    # Filter: diabetes gene set
    mask_gene = df['DB_Object_Symbol'].isin(diabetes_genes)
    df = df[mask_gene].copy()
    print(f"  Diabetes protein rows : {len(df):,}")

    if df.empty:
        print("  WARNING: No diabetes annotations found. Check gene set.")
        return df

    # Keep relevant columns
    keep = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier',
            'GO_ID', 'Evidence_Code', 'Aspect', 'Object_Name', 'Taxon']
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    # Normalise: strip whitespace
    for col in ['DB_Object_ID', 'DB_Object_Symbol', 'GO_ID', 'Aspect']:
        df[col] = df[col].str.strip()

    df.to_csv(GOA_CSV, index=False)
    print(f"  Saved {len(df):,} annotations → {GOA_CSV}")
    print(f"  Aspect breakdown:\n{df['Aspect'].value_counts().to_string()}")
    print(f"  Top evidence codes:\n{df['Evidence_Code'].value_counts().head(6).to_string()}\n")
    return df


if __name__ == "__main__":
    extract()
