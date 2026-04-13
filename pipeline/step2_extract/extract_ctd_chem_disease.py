"""
Step 2.10 — Extract CTD chemical-disease associations for diabetes.
Requires manual download: CTD_chemicals_diseases.csv.gz
  https://ctdbase.org/downloads/ → CTD_chemicals_diseases.csv.gz
Output: output/extracted/ctd_chemicals_diabetes.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import CTD_FILE, CTD_CSV, EXTRACTED_DIR, CTD_CHEM_DISEASE_COLS, DIABETES_MESH_IDS

# CTD disease IDs to keep (MESH: prefix in file)
DIABETES_CTD_IDS = {f"MESH:{mid}" for mid in DIABETES_MESH_IDS}


def extract():
    print("\n[Step 2.10] Extracting CTD chemical-disease associations (diabetes)")
    print(f"  Source: {CTD_FILE}")

    if not CTD_FILE.exists():
        print(f"  WARNING: {CTD_FILE} not found.")
        print("  Download from: https://ctdbase.org/downloads/")
        print("  → CTD_chemicals_diseases.csv.gz  (free, no registration)")
        print("  Place file in: d:\\ali_bioinfo\\  then re-run this script.\n")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    # CTD files have 29 comment lines starting with '#', no header in data
    df = pd.read_csv(
        CTD_FILE,
        sep=',',
        comment='#',
        header=None,
        names=CTD_CHEM_DISEASE_COLS,
        low_memory=False,
    )
    print(f"  Total rows loaded   : {len(df):,}")

    # Filter: diabetes MeSH disease IDs
    mask = (
        df['DiseaseID'].str.startswith('MESH:D003920', na=False) |
        df['DiseaseID'].str.startswith('MESH:D003924', na=False) |
        df['DiseaseName'].str.contains('diabetes', case=False, na=False)
    )
    df = df[mask].copy()
    print(f"  Diabetes rows       : {len(df):,}")

    if df.empty:
        print("  WARNING: No diabetes chemical-disease rows. Check file format.")
        return df

    # Clean
    for col in ['ChemicalName', 'DiseaseName', 'DiseaseID']:
        df[col] = df[col].str.strip()
    df.dropna(subset=['ChemicalName', 'DiseaseID'], inplace=True)

    # Keep relevant columns
    keep = ['ChemicalName', 'ChemicalID', 'CasRN', 'DiseaseName', 'DiseaseID',
            'DirectEvidence', 'InferenceScore']
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    df.to_csv(CTD_CSV, index=False)
    print(f"  Saved {len(df):,} rows → {CTD_CSV}")
    print(f"  Unique chemicals    : {df['ChemicalName'].nunique():,}")
    print(f"  Disease IDs present : {df['DiseaseID'].value_counts().head(5).to_string()}")
    print(f"  Top chemicals:\n{df['ChemicalName'].value_counts().head(8).to_string()}\n")
    return df


if __name__ == "__main__":
    extract()
