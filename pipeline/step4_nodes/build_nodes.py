"""
Step 4 — Build unified node table.
Creates: output/integrated/diabeteskg_nodes.csv
Columns: node_id, node_type, node_name, attributes_json
Node ID format: "Type:SOURCE_ID"
  Gene:3630, Protein:P01308, GO:GO:0006006,
  Phenotype:HP:0000819, Chemical:CHEBI:6801, Disease:C0011847, Pathway:R-HSA-264876
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    GENE_INFO_CSV, GENE_MAP_CSV, GOA_CSV, PROTEIN_MAP_CSV,
    GO_TERMS_CSV, HP_TERMS_CSV, CHEBI_CSV, DISGENET_CSV, REACTOME_CSV, CTD_CSV,
    NODES_CSV, INTEGRATED_DIR, GO_NAMESPACE_TO_NODE_TYPE,
    DIABETES_SEED_GENES, DIABETES_GENE_SET_FILE,
)


def _attrs(d: dict) -> str:
    """Convert dict to JSON string, dropping NaN/None values."""
    return json.dumps({k: str(v) for k, v in d.items()
                       if v is not None and str(v) not in ('nan', '', 'None')})


def build_gene_nodes() -> list:
    print("  Building Gene nodes ...")
    rows = []
    if not GENE_INFO_CSV.exists():
        print("    WARNING: gene_info_human.csv missing — using seed genes only")
        for sym in DIABETES_SEED_GENES:
            rows.append({'node_id': f"Gene:{sym}", 'node_type': 'Gene',
                         'node_name': sym, 'attributes_json': _attrs({'gene_symbol': sym})})
        return rows

    # Load expanded diabetes gene set
    diabetes_genes = DIABETES_SEED_GENES.copy()
    if DIABETES_GENE_SET_FILE.exists():
        with open(DIABETES_GENE_SET_FILE) as f:
            diabetes_genes = set(line.strip() for line in f if line.strip())

    df = pd.read_csv(GENE_INFO_CSV)
    df = df[df['Symbol'].isin(diabetes_genes)].copy()
    # Also ensure all seed genes are present (may be missing from NCBI file)
    have = set(df['Symbol'].tolist())
    for sym in DIABETES_SEED_GENES - have:
        df = pd.concat([df, pd.DataFrame([{'Symbol': sym, 'GeneID': '',
                                            'chromosome': '', 'description': '',
                                            'type_of_gene': 'protein-coding'}])],
                       ignore_index=True)
    for _, r in df.iterrows():
        rows.append({
            'node_id':   f"Gene:{r['GeneID']}" if str(r.get('GeneID', '')) not in ('', 'nan') else f"Gene:{r['Symbol']}",
            'node_type': 'Gene',
            'node_name': str(r['Symbol']),
            'attributes_json': _attrs({
                'gene_symbol':  r['Symbol'],
                'entrez_id':    r.get('GeneID', ''),
                'chromosome':   r.get('chromosome', ''),
                'description':  r.get('description', ''),
                'type_of_gene': r.get('type_of_gene', ''),
            }),
        })
    print(f"    Gene nodes : {len(rows):,}")
    return rows


def build_protein_nodes() -> list:
    print("  Building Protein nodes ...")
    rows = []
    if not GOA_CSV.exists():
        print("    WARNING: goa_human_diabetes.csv missing")
        return rows

    goa = pd.read_csv(GOA_CSV, usecols=['DB_Object_ID', 'DB_Object_Symbol', 'Object_Name'])
    goa.drop_duplicates(subset=['DB_Object_ID'], inplace=True)
    for _, r in goa.iterrows():
        rows.append({
            'node_id':   f"Protein:{r['DB_Object_ID']}",
            'node_type': 'Protein',
            'node_name': str(r['DB_Object_Symbol']),
            'attributes_json': _attrs({
                'uniprot_id':   r['DB_Object_ID'],
                'gene_symbol':  r['DB_Object_Symbol'],
                'protein_name': r.get('Object_Name', ''),
                'taxon_id':     '9606',
            }),
        })
    print(f"    Protein nodes : {len(rows):,}")
    return rows


def build_go_nodes() -> list:
    print("  Building GO nodes (BiologicalProcess / MolecularFunction / CellularComponent) ...")
    rows = []
    if not GO_TERMS_CSV.exists():
        print("    WARNING: go_terms_diabetes.csv missing")
        return rows

    df = pd.read_csv(GO_TERMS_CSV)
    for _, r in df.iterrows():
        ns = str(r.get('namespace', ''))
        node_type = GO_NAMESPACE_TO_NODE_TYPE.get(ns, 'BiologicalProcess')
        rows.append({
            'node_id':   f"GO:{r['go_id']}",
            'node_type': node_type,
            'node_name': str(r['go_name']),
            'attributes_json': _attrs({
                'go_id':     r['go_id'],
                'go_name':   r['go_name'],
                'namespace': ns,
            }),
        })
    print(f"    GO nodes : {len(rows):,}")
    return rows


def build_phenotype_nodes() -> list:
    print("  Building Phenotype nodes ...")
    rows = []
    if not HP_TERMS_CSV.exists():
        print("    WARNING: hp_diabetes_terms.csv missing")
        return rows

    df = pd.read_csv(HP_TERMS_CSV)
    for _, r in df.iterrows():
        rows.append({
            'node_id':   f"Phenotype:{r['hp_id']}",
            'node_type': 'Phenotype',
            'node_name': str(r['hp_name']),
            'attributes_json': _attrs({
                'hp_id':      r['hp_id'],
                'hp_name':    r['hp_name'],
                'definition': str(r.get('definition', ''))[:300],
                'synonyms':   str(r.get('synonyms', '')),
            }),
        })
    print(f"    Phenotype nodes : {len(rows):,}")
    return rows


def build_chemical_nodes() -> list:
    print("  Building Chemical nodes ...")
    rows = []
    if not CHEBI_CSV.exists():
        print("    WARNING: chebi_diabetes.csv missing")
        return rows

    df = pd.read_csv(CHEBI_CSV)
    for _, r in df.iterrows():
        rows.append({
            'node_id':   f"Chemical:{r['chebi_id']}",
            'node_type': 'Chemical',
            'node_name': str(r['chebi_name']) if str(r.get('drug_name', '')) == '' else str(r['drug_name']),
            'attributes_json': _attrs({
                'chebi_id':       r['chebi_id'],
                'chebi_name':     r['chebi_name'],
                'drug_name':      r.get('drug_name', ''),
                'definition':     str(r.get('definition', ''))[:200],
                'parent_chebi_id': r.get('parent_chebi_id', ''),
            }),
        })
    print(f"    Chemical nodes : {len(rows):,}")
    return rows


def build_disease_nodes() -> list:
    print("  Building Disease nodes ...")
    rows = []
    seen = set()

    for src_csv, id_col, name_col in [
        (DISGENET_CSV, 'diseaseId', 'diseaseName'),
        (CTD_CSV,      'DiseaseID', 'DiseaseName'),
    ]:
        if not src_csv.exists():
            continue
        df = pd.read_csv(src_csv)
        if id_col not in df.columns:
            continue
        for _, r in df[[id_col, name_col]].drop_duplicates().iterrows():
            did  = str(r[id_col]).strip()
            name = str(r[name_col]).strip()
            nid  = f"Disease:{did}"
            if nid in seen:
                continue
            seen.add(nid)
            rows.append({
                'node_id':   nid,
                'node_type': 'Disease',
                'node_name': name,
                'attributes_json': _attrs({
                    'disease_id':   did,
                    'disease_name': name,
                }),
            })
    print(f"    Disease nodes : {len(rows):,}")
    return rows


def build_pathway_nodes() -> list:
    print("  Building Pathway nodes ...")
    rows = []
    if not REACTOME_CSV.exists():
        print("    WARNING: reactome_diabetes.csv missing")
        return rows

    df = pd.read_csv(REACTOME_CSV)
    df.drop_duplicates(subset=['Reactome_pathway_ID'], inplace=True)
    for _, r in df.iterrows():
        rows.append({
            'node_id':   f"Pathway:{r['Reactome_pathway_ID']}",
            'node_type': 'Pathway',
            'node_name': str(r['pathway_name']),
            'attributes_json': _attrs({
                'reactome_id':  r['Reactome_pathway_ID'],
                'pathway_name': r['pathway_name'],
                'species':      r.get('species', 'Homo sapiens'),
            }),
        })
    print(f"    Pathway nodes : {len(rows):,}")
    return rows


def build_nodes():
    print("\n[Step 4] Building unified node table")
    INTEGRATED_DIR.mkdir(parents=True, exist_ok=True)

    all_nodes = (
        build_gene_nodes() +
        build_protein_nodes() +
        build_go_nodes() +
        build_phenotype_nodes() +
        build_chemical_nodes() +
        build_disease_nodes() +
        build_pathway_nodes()
    )

    df = pd.DataFrame(all_nodes)
    df.drop_duplicates(subset=['node_id'], inplace=True)

    df.to_csv(NODES_CSV, index=False)
    print(f"\n  Total nodes saved : {len(df):,}")
    print(f"  By type:\n{df['node_type'].value_counts().to_string()}")
    print(f"  Saved → {NODES_CSV}\n")
    return df


if __name__ == "__main__":
    build_nodes()
