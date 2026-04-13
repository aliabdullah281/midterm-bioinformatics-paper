"""
Step 3 — ID Harmonization.
Builds cross-reference mapping tables:
  - gene_map.csv      : GeneSymbol ↔ Entrez GeneID
  - protein_map.csv   : UniProt_ID ↔ GeneSymbol
  - chemical_map.csv  : ChEBI_ID ↔ CTD ChemicalID ↔ ChemicalName
  - disease_map.csv   : UMLS_ID ↔ MeSH_ID ↔ DiseaseName
Output: output/integrated/gene_map.csv, protein_map.csv, chemical_map.csv, disease_map.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    GENE_INFO_CSV, GOA_CSV, CHEBI_CSV, CTD_CSV, DISGENET_CSV,
    GENE_MAP_CSV, PROTEIN_MAP_CSV, CHEMICAL_MAP_CSV, DISEASE_MAP_CSV,
    INTEGRATED_DIR, DIABETES_SEED_GENES, DIABETES_UMLS_IDS, CHEBI_DIABETES_DRUGS,
)


def build_gene_map():
    """GeneSymbol → Entrez GeneID from NCBI gene_info."""
    print("  [3.1] Building gene_map (Symbol ↔ GeneID) ...")
    if not GENE_INFO_CSV.exists():
        print("    WARNING: gene_info_human.csv not found — empty gene map")
        return pd.DataFrame(columns=['Symbol', 'GeneID', 'chromosome', 'description', 'type_of_gene'])

    df = pd.read_csv(GENE_INFO_CSV)
    df = df[['Symbol', 'GeneID', 'chromosome', 'description', 'type_of_gene']].copy()
    df.drop_duplicates(subset=['Symbol'], inplace=True)
    df.to_csv(GENE_MAP_CSV, index=False)
    print(f"    Saved {len(df):,} gene mappings → {GENE_MAP_CSV}")
    return df


def build_protein_map():
    """UniProt_ID ↔ GeneSymbol from GOA annotations."""
    print("  [3.2] Building protein_map (UniProt ↔ Symbol) ...")
    if not GOA_CSV.exists():
        print("    WARNING: goa_human_diabetes.csv not found — empty protein map")
        return pd.DataFrame(columns=['uniprot_id', 'gene_symbol', 'protein_name'])

    goa = pd.read_csv(GOA_CSV, usecols=['DB_Object_ID', 'DB_Object_Symbol', 'Object_Name'])
    goa.rename(columns={
        'DB_Object_ID':     'uniprot_id',
        'DB_Object_Symbol': 'gene_symbol',
        'Object_Name':      'protein_name',
    }, inplace=True)
    goa.dropna(subset=['uniprot_id', 'gene_symbol'], inplace=True)
    goa.drop_duplicates(subset=['uniprot_id'], inplace=True)
    goa.to_csv(PROTEIN_MAP_CSV, index=False)
    print(f"    Saved {len(goa):,} protein mappings → {PROTEIN_MAP_CSV}")
    return goa


def build_chemical_map():
    """ChEBI_ID ↔ ChemicalName, with CTD cross-reference where available."""
    print("  [3.3] Building chemical_map (ChEBI ↔ CTD) ...")
    rows = []

    if CHEBI_CSV.exists():
        chebi = pd.read_csv(CHEBI_CSV)
        for _, r in chebi.iterrows():
            rows.append({
                'chebi_id':      r.get('chebi_id', ''),
                'chebi_name':    r.get('chebi_name', ''),
                'drug_name':     r.get('drug_name', ''),
                'ctd_chem_id':   '',
                'ctd_chem_name': '',
            })

    # Cross-reference to CTD by name matching
    if CTD_CSV.exists() and rows:
        ctd = pd.read_csv(CTD_CSV)
        ctd_name_to_id = dict(zip(
            ctd['ChemicalName'].str.lower().str.strip(),
            ctd['ChemicalID'].str.strip() if 'ChemicalID' in ctd.columns else [''] * len(ctd)
        ))
        for row in rows:
            name_key = row['chebi_name'].lower().strip()
            row['ctd_chem_id'] = ctd_name_to_id.get(name_key, '')
            if row['ctd_chem_id']:
                row['ctd_chem_name'] = row['chebi_name']

    df = pd.DataFrame(rows)
    df.to_csv(CHEMICAL_MAP_CSV, index=False)
    print(f"    Saved {len(df):,} chemical mappings → {CHEMICAL_MAP_CSV}")
    ctd_linked = df[df['ctd_chem_id'] != '']
    print(f"    ChEBI↔CTD links   : {len(ctd_linked):,}")
    return df


def build_disease_map():
    """UMLS_ID ↔ MeSH_ID ↔ DiseaseName from DisGeNET + CTD."""
    print("  [3.4] Building disease_map (UMLS ↔ MeSH ↔ Name) ...")
    rows = {}

    # From DisGeNET
    if DISGENET_CSV.exists():
        dis = pd.read_csv(DISGENET_CSV)
        for _, r in dis[['diseaseId', 'diseaseName']].drop_duplicates().iterrows():
            did = str(r['diseaseId']).strip()
            rows[did] = {
                'disease_id':   did,
                'disease_name': str(r['diseaseName']).strip(),
                'id_type':      'UMLS',
                'mesh_id':      '',
                'source':       'DisGeNET',
            }

    # From CTD
    if CTD_CSV.exists():
        ctd = pd.read_csv(CTD_CSV)
        for _, r in ctd[['DiseaseID', 'DiseaseName']].drop_duplicates().iterrows():
            did = str(r['DiseaseID']).strip()
            if did not in rows:
                rows[did] = {
                    'disease_id':   did,
                    'disease_name': str(r['DiseaseName']).strip(),
                    'id_type':      'MeSH',
                    'mesh_id':      did.replace('MESH:', ''),
                    'source':       'CTD',
                }

    if not rows:
        print("    WARNING: No disease data found (DisGeNET and CTD not available).")
        return pd.DataFrame(columns=['disease_id', 'disease_name', 'id_type', 'mesh_id', 'source'])

    df = pd.DataFrame(list(rows.values()))
    df.to_csv(DISEASE_MAP_CSV, index=False)
    print(f"    Saved {len(df):,} disease mappings → {DISEASE_MAP_CSV}")
    print(f"    ID types:\n{df['id_type'].value_counts().to_string()}")
    return df


def harmonize():
    print("\n[Step 3] ID Harmonization")
    INTEGRATED_DIR.mkdir(parents=True, exist_ok=True)
    gene_map    = build_gene_map()
    protein_map = build_protein_map()
    chem_map    = build_chemical_map()
    dis_map     = build_disease_map()
    print(f"\n  Summary:")
    print(f"    Genes    : {len(gene_map):,}")
    print(f"    Proteins : {len(protein_map):,}")
    print(f"    Chemicals: {len(chem_map):,}")
    print(f"    Diseases : {len(dis_map):,}")
    return gene_map, protein_map, chem_map, dis_map


if __name__ == "__main__":
    harmonize()
