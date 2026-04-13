"""
Step 5 — Build unified edge table.
Creates: output/integrated/diabeteskg_edges.csv
Columns: source_id, target_id, edge_type, data_source, confidence_score, evidence_count
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    GOA_CSV, BIOGRID_CSV, DISGENET_CSV, REACTOME_CSV, CTD_CSV,
    GO_HIER_CSV, HP_HIER_CSV, CHEBI_CSV,
    PROTEIN_MAP_CSV, GENE_MAP_CSV, NODES_CSV,
    EDGES_CSV, INTEGRATED_DIR,
    GO_NAMESPACE_TO_NODE_TYPE,
)


def _make_edge(src, tgt, etype, source, score=None, evidence=None):
    return {
        'source_id':       str(src).strip(),
        'target_id':       str(tgt).strip(),
        'edge_type':       etype,
        'data_source':     source,
        'confidence_score': float(score) if score is not None else None,
        'evidence_count':  int(evidence) if evidence is not None else None,
    }


def _load_protein_map():
    """Returns dict: gene_symbol → uniprot_id"""
    if PROTEIN_MAP_CSV.exists():
        pm = pd.read_csv(PROTEIN_MAP_CSV)
        return dict(zip(pm['gene_symbol'].str.strip(), pm['uniprot_id'].str.strip()))
    if GOA_CSV.exists():
        goa = pd.read_csv(GOA_CSV, usecols=['DB_Object_ID', 'DB_Object_Symbol'])
        goa.drop_duplicates(subset=['DB_Object_Symbol'], inplace=True)
        return dict(zip(goa['DB_Object_Symbol'].str.strip(), goa['DB_Object_ID'].str.strip()))
    return {}


def _load_gene_map():
    """Returns dict: gene_symbol → node_id (Gene:GeneID or Gene:Symbol)"""
    if GENE_MAP_CSV.exists():
        gm = pd.read_csv(GENE_MAP_CSV)
        result = {}
        for _, r in gm.iterrows():
            gid = str(r.get('GeneID', '')).strip()
            sym = str(r['Symbol']).strip()
            result[sym] = f"Gene:{gid}" if gid and gid != 'nan' else f"Gene:{sym}"
        return result
    return {sym: f"Gene:{sym}" for sym in []}


def build_goa_edges() -> list:
    """protein_has_biological_process, protein_has_molecular_function,
       protein_located_in_component, gene_encodes_protein"""
    print("  [5.1] GOA annotation edges ...")
    if not GOA_CSV.exists():
        print("    WARNING: goa_human_diabetes.csv missing")
        return []

    goa = pd.read_csv(GOA_CSV)
    rows = []
    aspect_to_etype = {
        'P': 'protein_has_biological_process',
        'F': 'protein_has_molecular_function',
        'C': 'protein_located_in_component',
    }
    # Protein → GO edges
    for _, r in goa.iterrows():
        etype = aspect_to_etype.get(str(r.get('Aspect', '')).strip())
        if etype is None:
            continue
        rows.append(_make_edge(
            f"Protein:{r['DB_Object_ID']}",
            f"GO:{r['GO_ID']}",
            etype, 'GOA',
        ))

    # gene_encodes_protein edges (unique gene→protein pairs)
    gene_prot = goa[['DB_Object_Symbol', 'DB_Object_ID']].drop_duplicates()
    gene_map = _load_gene_map()
    for _, r in gene_prot.iterrows():
        sym = str(r['DB_Object_Symbol']).strip()
        gene_nid = gene_map.get(sym, f"Gene:{sym}")
        rows.append(_make_edge(
            gene_nid,
            f"Protein:{r['DB_Object_ID']}",
            'gene_encodes_protein', 'GOA',
        ))

    print(f"    GOA edges : {len(rows):,}")
    return rows


def build_biogrid_edges() -> list:
    """protein_interacts_with_protein"""
    print("  [5.2] BioGRID PPI edges ...")
    if not BIOGRID_CSV.exists():
        print("    WARNING: biogrid_human_diabetes.csv missing")
        return []

    df = pd.read_csv(BIOGRID_CSV)
    prot_map = _load_protein_map()
    rows = []
    missing = 0
    for _, r in df.iterrows():
        sym_a = str(r.get('symbol_a', '')).strip()
        sym_b = str(r.get('symbol_b', '')).strip()
        pid_a = prot_map.get(sym_a)
        pid_b = prot_map.get(sym_b)
        if not pid_a:
            missing += 1
            pid_a = f"Gene:{sym_a}"   # fall back to Gene node
        if not pid_b:
            pid_b = f"Gene:{sym_b}"

        score_raw = r.get('Score', None)
        try:
            score = float(score_raw)
        except (TypeError, ValueError):
            score = None

        rows.append(_make_edge(
            f"Protein:{pid_a}" if pid_a.startswith('P') else pid_a,
            f"Protein:{pid_b}" if pid_b.startswith('P') else pid_b,
            'protein_interacts_with_protein', 'BioGRID', score,
        ))

    print(f"    BioGRID edges  : {len(rows):,}  (proteins not in GOA: {missing:,})")
    return rows


def build_disgenet_edges() -> list:
    """gene_associated_with_disease"""
    print("  [5.3] DisGeNET gene-disease edges ...")
    if not DISGENET_CSV.exists():
        print("    WARNING: disgenet_diabetes.csv missing")
        return []

    df = pd.read_csv(DISGENET_CSV)
    gene_map = _load_gene_map()
    rows = []
    for _, r in df.iterrows():
        sym = str(r.get('geneSymbol', '')).strip()
        did = str(r.get('diseaseId', '')).strip()
        gene_nid = gene_map.get(sym, f"Gene:{sym}")
        score = r.get('score', None)
        pmids = r.get('NofPmids', None)
        rows.append(_make_edge(gene_nid, f"Disease:{did}",
                               'gene_associated_with_disease', 'DisGeNET',
                               score, pmids))
    print(f"    DisGeNET edges : {len(rows):,}")
    return rows


def build_reactome_edges() -> list:
    """protein_participates_in_pathway"""
    print("  [5.4] Reactome protein-pathway edges ...")
    if not REACTOME_CSV.exists():
        print("    WARNING: reactome_diabetes.csv missing")
        return []

    df = pd.read_csv(REACTOME_CSV)
    rows = []
    for _, r in df.iterrows():
        rows.append(_make_edge(
            f"Protein:{r['UniProt_ID']}",
            f"Pathway:{r['Reactome_pathway_ID']}",
            'protein_participates_in_pathway', 'Reactome',
        ))
    print(f"    Reactome edges : {len(rows):,}")
    return rows


def build_ctd_edges() -> list:
    """chemical_associated_with_disease"""
    print("  [5.5] CTD chemical-disease edges ...")
    if not CTD_CSV.exists():
        print("    WARNING: ctd_chemicals_diabetes.csv missing")
        return []

    df = pd.read_csv(CTD_CSV)
    rows = []
    for _, r in df.iterrows():
        cname = str(r.get('ChemicalName', '')).strip()
        did   = str(r.get('DiseaseID', '')).strip()
        score = r.get('InferenceScore', None)
        # Use ChemicalID (CTD mesh-style) or name as chemical node id
        cid_raw = str(r.get('ChemicalID', '')).strip()
        cid = f"Chemical:CTD:{cid_raw}" if cid_raw and cid_raw != 'nan' else f"Chemical:CTD:{cname}"
        rows.append(_make_edge(cid, f"Disease:{did}",
                               'chemical_associated_with_disease', 'CTD', score))
    print(f"    CTD edges      : {len(rows):,}")
    return rows


def build_go_hierarchy_edges() -> list:
    """go_term_is_a, go_term_part_of"""
    print("  [5.6] GO hierarchy edges ...")
    if not GO_HIER_CSV.exists():
        print("    WARNING: go_hierarchy_diabetes.csv missing")
        return []

    df = pd.read_csv(GO_HIER_CSV)
    rows = []
    for _, r in df.iterrows():
        etype = 'go_term_is_a' if r.get('relation') == 'is_a' else 'go_term_part_of'
        rows.append(_make_edge(
            f"GO:{r['child_id']}", f"GO:{r['parent_id']}",
            etype, 'GeneOntology',
        ))
    print(f"    GO hierarchy edges : {len(rows):,}")
    return rows


def build_hp_hierarchy_edges() -> list:
    """phenotype_is_a"""
    print("  [5.7] HPO hierarchy edges ...")
    if not HP_HIER_CSV.exists():
        print("    WARNING: hp_diabetes_hierarchy.csv missing")
        return []

    df = pd.read_csv(HP_HIER_CSV)
    rows = []
    for _, r in df.iterrows():
        rows.append(_make_edge(
            f"Phenotype:{r['child_id']}", f"Phenotype:{r['parent_id']}",
            'phenotype_is_a', 'HPO',
        ))
    print(f"    HP hierarchy edges : {len(rows):,}")
    return rows


def build_chebi_hierarchy_edges() -> list:
    """chemical_is_a_subclass_of"""
    print("  [5.8] ChEBI hierarchy edges ...")
    if not CHEBI_CSV.exists():
        print("    WARNING: chebi_diabetes.csv missing")
        return []

    df = pd.read_csv(CHEBI_CSV)
    if 'parent_chebi_id' not in df.columns:
        return []
    rows = []
    for _, r in df.iterrows():
        pid = str(r.get('parent_chebi_id', '')).strip()
        if pid and pid != 'nan':
            rows.append(_make_edge(
                f"Chemical:{r['chebi_id']}", f"Chemical:{pid}",
                'chemical_is_a_subclass_of', 'ChEBI',
            ))
    print(f"    ChEBI hierarchy edges : {len(rows):,}")
    return rows


def build_edges():
    print("\n[Step 5] Building unified edge table")
    INTEGRATED_DIR.mkdir(parents=True, exist_ok=True)

    all_edges = (
        build_goa_edges() +
        build_biogrid_edges() +
        build_disgenet_edges() +
        build_reactome_edges() +
        build_ctd_edges() +
        build_go_hierarchy_edges() +
        build_hp_hierarchy_edges() +
        build_chebi_hierarchy_edges()
    )

    df = pd.DataFrame(all_edges)
    before = len(df)

    # Deduplicate: same (src, tgt, edge_type) → keep max confidence, sum evidence
    if not df.empty:
        df['confidence_score'] = pd.to_numeric(df['confidence_score'], errors='coerce')
        df['evidence_count']   = pd.to_numeric(df['evidence_count'],   errors='coerce')
        agg = {
            'data_source':     lambda x: '|'.join(sorted(set(x.dropna().astype(str)))),
            'confidence_score': 'max',
            'evidence_count':  'sum',
        }
        df = (
            df.groupby(['source_id', 'target_id', 'edge_type'], as_index=False)
            .agg(agg)
        )
    after = len(df)
    print(f"\n  Edges before dedup : {before:,}")
    print(f"  Edges after dedup  : {after:,}")

    df.to_csv(EDGES_CSV, index=False)
    print(f"  By type:\n{df['edge_type'].value_counts().to_string()}")
    print(f"  Saved → {EDGES_CSV}\n")
    return df


if __name__ == "__main__":
    build_edges()
