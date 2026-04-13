"""
Step 2.6 — Extract NCBI gene_info for human genes.
Run this FIRST — it provides the authoritative Gene Symbol ↔ GeneID mapping.
Output: output/extracted/gene_info_human.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import GENE_INFO, GENE_INFO_CSV, EXTRACTED_DIR


def extract():
    print("\n[Step 2.6] Extracting NCBI gene_info (human only)")
    print(f"  Source: {GENE_INFO}")

    if not GENE_INFO.exists():
        print(f"  ERROR: {GENE_INFO} not found. Skipping.")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    # Header line starts with "#tax_id" — read normally then rename
    df = pd.read_csv(
        GENE_INFO,
        sep='\t',
        compression='gzip',
        header=0,
        low_memory=False,
    )
    # Rename the '#tax_id' column
    df.rename(columns={'#tax_id': 'tax_id'}, inplace=True)
    print(f"  Total rows loaded   : {len(df):,}")

    # Filter to human
    df = df[df['tax_id'] == 9606].copy()
    print(f"  Human rows (9606)   : {len(df):,}")

    # Keep relevant columns
    keep = ['tax_id', 'GeneID', 'Symbol', 'chromosome', 'description', 'type_of_gene']
    existing = [c for c in keep if c in df.columns]
    df = df[existing].copy()

    # Clean
    df['Symbol'] = df['Symbol'].str.strip()
    df.dropna(subset=['Symbol', 'GeneID'], inplace=True)
    df.drop_duplicates(subset=['GeneID'], inplace=True)

    df.to_csv(GENE_INFO_CSV, index=False)
    print(f"  Saved {len(df):,} genes → {GENE_INFO_CSV}")
    print(f"  Sample:\n{df.head(3).to_string()}\n")
    return df


if __name__ == "__main__":
    extract()
