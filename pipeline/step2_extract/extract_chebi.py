"""
Step 2.3 — Parse ChEBI OWL and extract diabetes-relevant chemical terms + hierarchy.
Starts from CHEBI_DIABETES_DRUGS seeds and walks subClassOf ancestors.
Output: output/extracted/chebi_diabetes.csv
"""
import sys, re, time
from pathlib import Path
from collections import deque
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from rdflib import Graph, URIRef
from rdflib.namespace import RDF, RDFS, OWL
from config import CHEBI_OWL, CHEBI_CSV, EXTRACTED_DIR, CHEBI_DIABETES_DRUGS

OBO = "http://purl.obolibrary.org/obo/"
IAO_DEF = URIRef(f"http://purl.obolibrary.org/obo/IAO_0000115")


def _iri_to_chebi(iri: str):
    """'http://…/CHEBI_6801' → 'CHEBI:6801'"""
    m = re.search(r'CHEBI_(\d+)', str(iri))
    return f"CHEBI:{m.group(1)}" if m else None


def extract():
    print("\n[Step 2.3] Parsing ChEBI OWL (diabetes drug terms + ancestors)")
    print(f"  Source: {CHEBI_OWL}")

    if not CHEBI_OWL.exists():
        print(f"  ERROR: {CHEBI_OWL} not found. Skipping.")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    print("  Parsing ChEBI OWL with rdflib ...")
    t0 = time.time()
    g = Graph()
    g.parse(str(CHEBI_OWL), format='xml')
    print(f"  Parsed {len(g):,} triples in {time.time()-t0:.1f}s")

    # Build IRI → ChEBI_ID map for all classes
    iri_to_id = {}
    id_to_iri = {}
    for cls in g.subjects(RDF.type, OWL.Class):
        cid = _iri_to_chebi(str(cls))
        if cid:
            iri_to_id[str(cls)] = cid
            id_to_iri[cid] = str(cls)

    print(f"  Total ChEBI classes : {len(iri_to_id):,}")

    # Build parent map: chebi_id → list of parent chebi_ids
    parent_map = {}
    for cls_iri, cid in iri_to_id.items():
        cls = URIRef(cls_iri)
        parents = []
        for sc in g.objects(cls, RDFS.subClassOf):
            if isinstance(sc, URIRef):
                pid = iri_to_id.get(str(sc))
                if pid:
                    parents.append(pid)
        parent_map[cid] = parents

    # BFS from seed drugs → collect all ancestors
    seed_ids = set(CHEBI_DIABETES_DRUGS.keys())
    to_extract = set(seed_ids)
    queue = deque(seed_ids)
    while queue:
        cid = queue.popleft()
        for pid in parent_map.get(cid, []):
            if pid not in to_extract:
                to_extract.add(pid)
                queue.append(pid)

    print(f"  Seeds + ancestors   : {len(to_extract):,}")

    # Extract term attributes
    rows = []
    hier = []
    for cid in to_extract:
        cls_iri_str = id_to_iri.get(cid)
        if not cls_iri_str:
            continue
        cls = URIRef(cls_iri_str)

        labels = list(g.objects(cls, RDFS.label))
        name = str(labels[0]) if labels else ""

        defs = list(g.objects(cls, IAO_DEF))
        definition = str(defs[0]) if defs else ""

        rows.append({
            'chebi_id':   cid,
            'chebi_name': name,
            'definition': definition[:300] if definition else "",
            'is_drug_seed': cid in seed_ids,
        })

        for pid in parent_map.get(cid, []):
            if pid in to_extract:
                hier.append({'child_id': cid, 'parent_id': pid,
                             'relation': 'is_a'})

    df = pd.DataFrame(rows)
    df_hier = pd.DataFrame(hier)

    # Annotate known drug names from seed list
    df['drug_name'] = df['chebi_id'].map(CHEBI_DIABETES_DRUGS).fillna('')

    df.to_csv(CHEBI_CSV, index=False)
    print(f"  Saved {len(df):,} ChEBI terms → {CHEBI_CSV}")
    print(f"  Hierarchy edges     : {len(df_hier):,}")

    # Embed hierarchy as parent_chebi_id (first parent per term, for node attribute)
    parent_col = {row['child_id']: row['parent_id'] for _, row in df_hier.iterrows()}
    df['parent_chebi_id'] = df['chebi_id'].map(parent_col).fillna('')
    df.to_csv(CHEBI_CSV, index=False)  # resave with parent column

    print(f"  Drug seeds present  : {df[df['is_drug_seed']]['chebi_name'].tolist()}")
    print(f"  Sample:\n{df[df['is_drug_seed']].head(5)[['chebi_id','chebi_name','drug_name']].to_string()}\n")
    return df


if __name__ == "__main__":
    extract()
