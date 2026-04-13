"""
Step 2.8 — Extract diabetes gene-disease associations from DisGeNET.
Requires manual download: curated_gene_disease_associations.tsv
  https://www.disgenet.org/downloads → curated_gene_disease_associations.tsv
Output: output/extracted/disgenet_diabetes.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import DISGENET_FILE, DISGENET_CSV, EXTRACTED_DIR, DIABETES_UMLS_IDS


def extract():
    print("\n[Step 2.8] Extracting DisGeNET diabetes associations")
    print(f"  Source: {DISGENET_FILE}")

    if not DISGENET_FILE.exists():
        print(f"  WARNING: {DISGENET_FILE} not found.")
        print("  Download from: https://www.disgenet.org/downloads")
        print("  → curated_gene_disease_associations.tsv (free registration required)")
        print("  Place file in: d:\\ali_bioinfo\\  then re-run this script.\n")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DISGENET_FILE, sep='\t', header=0, low_memory=False)
    print(f"  Total rows loaded   : {len(df):,}")
    print(f"  Columns: {list(df.columns)}")

    # Normalise column names (DisGeNET versions differ slightly)
    df.columns = [c.strip() for c in df.columns]
    col_map = {
        'geneid': 'geneId', 'genesymbol': 'geneSymbol',
        'diseaseid': 'diseaseId', 'diseasename': 'diseaseName',
    }
    df.rename(columns={c: col_map.get(c.lower(), c) for c in df.columns}, inplace=True)

    if 'diseaseId' not in df.columns:
        print("  ERROR: 'diseaseId' column not found. Check file format.")
        return None

    # Filter: diabetes UMLS IDs or disease name contains 'diabetes'
    mask_id = df['diseaseId'].isin(DIABETES_UMLS_IDS)
    mask_name = df['diseaseName'].str.contains('diabetes', case=False, na=False)
    df = df[mask_id | mask_name].copy()
    print(f"  Diabetes rows       : {len(df):,}")

    if df.empty:
        print("  WARNING: No diabetes associations found. Check UMLS IDs in file.")
        return df

    # Keep useful columns
    keep = ['geneId', 'geneSymbol', 'diseaseId', 'diseaseName', 'score', 'NofPmids', 'source']
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    for col in ['geneSymbol', 'diseaseId', 'diseaseName']:
        if col in df.columns:
            df[col] = df[col].str.strip()
    df.dropna(subset=['geneSymbol', 'diseaseId'], inplace=True)

    df.to_csv(DISGENET_CSV, index=False)
    print(f"  Saved {len(df):,} associations → {DISGENET_CSV}")
    print(f"  Disease breakdown:\n{df['diseaseName'].value_counts().head(6).to_string()}")
    print(f"  Top genes:\n{df['geneSymbol'].value_counts().head(6).to_string()}\n")
    return df


if __name__ == "__main__":
    extract()
