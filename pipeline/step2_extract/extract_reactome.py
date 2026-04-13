"""
Step 2.9 — Extract Reactome pathways for diabetes proteins.
Requires manual download: UniProt2Reactome_All_Levels.txt
  https://reactome.org/download-data → UniProt2Reactome_All_Levels.txt
Output: output/extracted/reactome_diabetes.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    REACTOME_FILE, REACTOME_CSV, GOA_CSV,
    EXTRACTED_DIR, REACTOME_COLS,
)

# Key insulin/glucose Reactome pathway IDs (Homo sapiens)
DIABETES_REACTOME_PATHWAYS = {
    'R-HSA-264876',   # Insulin signalling
    'R-HSA-70326',    # Glucose metabolism
    'R-HSA-422085',   # Synthesis, secretion and inactivation of Glucagon-like Peptide-1
    'R-HSA-6807878',  # COPI-dependent Golgi-to-ER retrograde traffic (beta-cell)
    'R-HSA-392499',   # Metabolism of proteins
    'R-HSA-1430728',  # Metabolism
}


def extract():
    print("\n[Step 2.9] Extracting Reactome pathways (diabetes proteins)")
    print(f"  Source: {REACTOME_FILE}")

    if not REACTOME_FILE.exists():
        print(f"  WARNING: {REACTOME_FILE} not found.")
        print("  Download from: https://reactome.org/download-data")
        print("  → UniProt2Reactome_All_Levels.txt  (open access)")
        print("  Place file in: d:\\ali_bioinfo\\  then re-run this script.\n")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    # Load diabetes protein set from GOA
    diabetes_proteins = set()
    if GOA_CSV.exists():
        goa = pd.read_csv(GOA_CSV, usecols=['DB_Object_ID'])
        diabetes_proteins = set(goa['DB_Object_ID'].dropna().str.strip().tolist())
        print(f"  Diabetes UniProt IDs (from GOA): {len(diabetes_proteins):,}")
    else:
        print("  WARNING: goa_human_diabetes.csv not found — no protein filter applied")

    df = pd.read_csv(
        REACTOME_FILE,
        sep='\t',
        header=None,
        names=REACTOME_COLS,
        low_memory=False,
    )
    print(f"  Total rows loaded   : {len(df):,}")

    # Filter: Homo sapiens only
    df = df[df['species'] == 'Homo sapiens'].copy()
    print(f"  Human rows          : {len(df):,}")

    # Filter: diabetes protein set (if available)
    if diabetes_proteins:
        mask_prot = df['UniProt_ID'].isin(diabetes_proteins)
        df_diabetes = df[mask_prot].copy()
        print(f"  Diabetes protein rows: {len(df_diabetes):,}")
    else:
        df_diabetes = df.copy()

    if df_diabetes.empty:
        print("  WARNING: No diabetes protein rows found.")
        return df_diabetes

    # Clean up
    for col in ['UniProt_ID', 'Reactome_pathway_ID', 'pathway_name']:
        df_diabetes[col] = df_diabetes[col].str.strip()

    df_diabetes.drop_duplicates(subset=['UniProt_ID', 'Reactome_pathway_ID'], inplace=True)

    df_diabetes.to_csv(REACTOME_CSV, index=False)
    print(f"  Saved {len(df_diabetes):,} protein→pathway rows → {REACTOME_CSV}")
    print(f"  Unique pathways     : {df_diabetes['Reactome_pathway_ID'].nunique():,}")
    print(f"  Top pathways:\n{df_diabetes['pathway_name'].value_counts().head(8).to_string()}\n")
    return df_diabetes


if __name__ == "__main__":
    extract()
