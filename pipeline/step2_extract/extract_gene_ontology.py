"""
Step 2.2 — Parse Gene Ontology OWL and extract GO terms + hierarchy for diabetes proteins.
Uses go-basic.owl (preferred for speed) or falls back to go.owl.
Depends on: goa_human_diabetes.csv (for filtering to relevant GO terms only)
Output: output/extracted/go_terms_diabetes.csv
         output/extracted/go_hierarchy_diabetes.csv
"""
import sys, re, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from rdflib import Graph, URIRef, BNode
from rdflib.namespace import RDF, RDFS, OWL
from config import (
    GO_BASIC_OWL, GO_OWL, GOA_CSV,
    GO_TERMS_CSV, GO_HIER_CSV, EXTRACTED_DIR,
    GO_NAMESPACE_TO_NODE_TYPE, DIABETES_GO_SEEDS,
)

OBO         = "http://purl.obolibrary.org/obo/"
OBOANN      = "http://www.geneontology.org/formats/oboInOwl#"
BFO_PART_OF = URIRef(f"{OBO}BFO_0000050")


def _go_iri_to_id(iri: str):
    """'http://…/GO_0006006' → 'GO:0006006'"""
    m = re.search(r'GO_(\d+)', iri)
    return f"GO:{m.group(1).zfill(7)}" if m else None


def extract():
    print("\n[Step 2.2] Parsing Gene Ontology OWL")

    # Choose OWL file
    owl_path = GO_BASIC_OWL if GO_BASIC_OWL.exists() else GO_OWL
    if not owl_path.exists():
        print(f"  ERROR: Neither go-basic.owl nor go.owl found. Skipping.")
        return None
    if owl_path == GO_OWL:
        print("  WARNING: go-basic.owl not found — using go.owl (may be slow).")
    print(f"  Source: {owl_path}")

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    # Determine relevant GO IDs from GOA annotations
    relevant_go_ids = set(DIABETES_GO_SEEDS)
    if GOA_CSV.exists():
        goa = pd.read_csv(GOA_CSV, usecols=['GO_ID'])
        relevant_go_ids |= set(goa['GO_ID'].dropna().str.strip().tolist())
        print(f"  GO IDs from GOA annotations : {len(relevant_go_ids):,}")
    else:
        print("  WARNING: goa_human_diabetes.csv not found — keeping all GO terms")
        relevant_go_ids = None  # keep all

    # Parse OWL with rdflib
    print("  Parsing OWL (rdflib) ...")
    t0 = time.time()
    g = Graph()
    g.parse(str(owl_path), format='xml')
    print(f"  Parsed {len(g):,} triples in {time.time()-t0:.1f}s")

    # Helper: get annotation value
    def ann(cls, prop_local):
        vals = list(g.objects(cls, URIRef(f"{OBOANN}{prop_local}")))
        return str(vals[0]) if vals else ""

    # Extract terms + hierarchy
    terms = []
    hier  = []

    for cls in g.subjects(RDF.type, OWL.Class):
        if not isinstance(cls, URIRef):
            continue
        go_id = _go_iri_to_id(str(cls))
        if go_id is None:
            continue

        # Get label
        labels = list(g.objects(cls, RDFS.label))
        go_name = str(labels[0]) if labels else ""

        # Get namespace
        ns = ann(cls, 'hasOBONamespace')
        node_type = GO_NAMESPACE_TO_NODE_TYPE.get(ns, "GO")

        is_obsolete = ann(cls, 'is_obsolete')
        if is_obsolete == 'true':
            continue

        terms.append({
            'go_id':     go_id,
            'go_name':   go_name,
            'namespace': ns,
            'node_type': node_type,
        })

        # Hierarchy edges
        for sc in g.objects(cls, RDFS.subClassOf):
            if isinstance(sc, URIRef):
                parent_id = _go_iri_to_id(str(sc))
                if parent_id:
                    hier.append({'child_id': go_id, 'parent_id': parent_id,
                                 'relation': 'is_a'})
            elif isinstance(sc, BNode):
                # OWL restriction: check if part_of
                for prop in g.objects(sc, OWL.onProperty):
                    if prop == BFO_PART_OF:
                        for val in g.objects(sc, OWL.someValuesFrom):
                            if isinstance(val, URIRef):
                                parent_id = _go_iri_to_id(str(val))
                                if parent_id:
                                    hier.append({
                                        'child_id': go_id,
                                        'parent_id': parent_id,
                                        'relation': 'part_of',
                                    })

    df_terms = pd.DataFrame(terms)
    df_hier  = pd.DataFrame(hier)
    print(f"  Total GO terms   : {len(df_terms):,}")
    print(f"  Total GO edges   : {len(df_hier):,}")

    # Filter to relevant GO IDs (terms + their parents collected above)
    if relevant_go_ids is not None:
        # Also collect all ancestors of relevant IDs from hierarchy
        child_to_parents = {}
        for _, row in df_hier.iterrows():
            child_to_parents.setdefault(row['child_id'], []).append(row['parent_id'])

        to_keep = set(relevant_go_ids)
        frontier = list(relevant_go_ids)
        while frontier:
            node = frontier.pop()
            for p in child_to_parents.get(node, []):
                if p not in to_keep:
                    to_keep.add(p)
                    frontier.append(p)

        df_terms = df_terms[df_terms['go_id'].isin(to_keep)].copy()
        df_hier  = df_hier[
            df_hier['child_id'].isin(to_keep) &
            df_hier['parent_id'].isin(to_keep)
        ].copy()
        print(f"  After filtering  : {len(df_terms):,} terms, {len(df_hier):,} edges")

    df_terms.to_csv(GO_TERMS_CSV, index=False)
    df_hier.to_csv(GO_HIER_CSV,   index=False)
    print(f"  Saved → {GO_TERMS_CSV}")
    print(f"  Saved → {GO_HIER_CSV}")
    print(f"  Namespace breakdown:\n{df_terms['namespace'].value_counts().to_string()}\n")
    return df_terms, df_hier


if __name__ == "__main__":
    extract()
