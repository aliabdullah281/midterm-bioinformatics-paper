"""
Step 2.5 — Parse Human Phenotype Ontology (OBO format) and extract diabetes subtree.
Uses pronto library to parse hp.obo.
Output: output/extracted/hp_diabetes_terms.csv
         output/extracted/hp_diabetes_hierarchy.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import HP_OBO, HP_TERMS_CSV, HP_HIER_CSV, EXTRACTED_DIR, HP_DIABETES_ROOT


def extract():
    print("\n[Step 2.5] Parsing Human Phenotype Ontology (OBO)")
    print(f"  Source: {HP_OBO}")

    if not HP_OBO.exists():
        print(f"  ERROR: {HP_OBO} not found. Skipping.")
        return None

    try:
        import pronto
    except ImportError:
        print("  ERROR: pronto not installed. Run: pip install pronto")
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    print("  Loading hp.obo with pronto ...")
    onto = pronto.Ontology(str(HP_OBO))
    print(f"  Total HPO terms     : {len(onto):,}")

    root_id = HP_DIABETES_ROOT
    if root_id not in onto:
        print(f"  ERROR: Root term {root_id} not found in ontology.")
        return None

    root = onto[root_id]
    print(f"  Root term: {root.id} — {root.name}")

    # Collect all descendants (includes root itself)
    subtree_terms = list(root.subclasses(with_self=True))
    print(f"  Diabetes subtree    : {len(subtree_terms):,} terms")

    # Extract term attributes
    rows = []
    hier = []

    for term in subtree_terms:
        if term.obsolete:
            continue

        # Definition
        definition = ""
        if term.definition:
            definition = str(term.definition)[:500]

        # Synonyms (pipe-separated)
        synonyms = "|".join(str(s.description) for s in term.synonyms)

        rows.append({
            'hp_id':       term.id,
            'hp_name':     term.name,
            'definition':  definition,
            'synonyms':    synonyms,
        })

        # is_a parents within the subtree
        subtree_ids = {t.id for t in subtree_terms}
        for parent in term.superclasses(distance=1, with_self=False):
            if parent.id in subtree_ids:
                hier.append({
                    'child_id':  term.id,
                    'parent_id': parent.id,
                    'relation':  'is_a',
                })

    df_terms = pd.DataFrame(rows)
    df_hier  = pd.DataFrame(hier)

    df_terms.to_csv(HP_TERMS_CSV, index=False)
    df_hier.to_csv(HP_HIER_CSV,   index=False)
    print(f"  Saved {len(df_terms):,} phenotype terms → {HP_TERMS_CSV}")
    print(f"  Saved {len(df_hier):,} hierarchy edges → {HP_HIER_CSV}")
    print(f"  Sample:\n{df_terms.head(5)[['hp_id','hp_name']].to_string()}\n")
    return df_terms, df_hier


if __name__ == "__main__":
    extract()
