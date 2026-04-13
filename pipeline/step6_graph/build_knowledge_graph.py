"""
Step 6 — Build DiabetesKG as a NetworkX heterogeneous graph.
Loads nodes.csv + edges.csv → builds MultiDiGraph → exports GraphML.
Output: output/graph/diabeteskg.graphml
         output/graph/diabeteskg_nodes.csv
         output/graph/diabeteskg_edges.csv
"""
import sys, json, shutil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
import networkx as nx
from config import (
    NODES_CSV, EDGES_CSV,
    GRAPHML_FILE, FINAL_NODES_CSV, FINAL_EDGES_CSV,
    GRAPH_DIR,
)


def build_kg():
    print("\n[Step 6] Building DiabetesKG as NetworkX MultiDiGraph")
    GRAPH_DIR.mkdir(parents=True, exist_ok=True)

    if not NODES_CSV.exists():
        print(f"  ERROR: {NODES_CSV} not found. Run Steps 3-5 first.")
        return None
    if not EDGES_CSV.exists():
        print(f"  ERROR: {EDGES_CSV} not found. Run Step 5 first.")
        return None

    nodes = pd.read_csv(NODES_CSV)
    edges = pd.read_csv(EDGES_CSV)
    print(f"  Nodes loaded : {len(nodes):,}")
    print(f"  Edges loaded : {len(edges):,}")

    G = nx.MultiDiGraph()

    # Add nodes
    print("  Adding nodes ...")
    for _, r in nodes.iterrows():
        nid   = str(r['node_id'])
        attrs = {}
        try:
            attrs = json.loads(r.get('attributes_json', '{}') or '{}')
        except Exception:
            pass
        attrs['node_type'] = str(r['node_type'])
        attrs['node_name'] = str(r['node_name'])
        G.add_node(nid, **attrs)

    # Add edges
    print("  Adding edges ...")
    skipped = 0
    for _, r in edges.iterrows():
        src = str(r['source_id'])
        tgt = str(r['target_id'])
        if src not in G or tgt not in G:
            skipped += 1
            continue
        edge_attrs = {
            'edge_type':       str(r.get('edge_type', '')),
            'data_source':     str(r.get('data_source', '')),
            'confidence_score': r.get('confidence_score', None),
            'evidence_count':  r.get('evidence_count', None),
        }
        G.add_edge(src, tgt, **edge_attrs)

    print(f"  Edges added  : {G.number_of_edges():,}  (skipped dangling: {skipped:,})")

    # ── Basic Statistics ────────────────────────────────────────────────────
    print("\n  === Graph Statistics ===")
    print(f"  Nodes total         : {G.number_of_nodes():,}")
    print(f"  Edges total         : {G.number_of_edges():,}")

    # Node counts per type
    type_counts = {}
    for n, attr in G.nodes(data=True):
        t = attr.get('node_type', 'Unknown')
        type_counts[t] = type_counts.get(t, 0) + 1
    print("\n  Nodes by type:")
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"    {t:<25} {c:>6,}")

    # Edge counts per type
    edge_counts = {}
    for u, v, d in G.edges(data=True):
        t = d.get('edge_type', 'Unknown')
        edge_counts[t] = edge_counts.get(t, 0) + 1
    print("\n  Edges by type:")
    for t, c in sorted(edge_counts.items(), key=lambda x: -x[1]):
        print(f"    {t:<45} {c:>6,}")

    # Connectivity
    undirected = G.to_undirected()
    components = list(nx.connected_components(undirected))
    print(f"\n  Connected components: {len(components):,}")
    if components:
        largest = max(components, key=len)
        print(f"  Largest component   : {len(largest):,} nodes "
              f"({100*len(largest)/G.number_of_nodes():.1f}%)")

    # Degree stats
    degrees = [d for _, d in G.degree()]
    if degrees:
        avg_deg = sum(degrees) / len(degrees)
        max_deg = max(degrees)
        print(f"  Average degree      : {avg_deg:.2f}")
        print(f"  Max degree          : {max_deg:,}")

    # ── Export ────────────────────────────────────────────────────────────
    print(f"\n  Exporting GraphML → {GRAPHML_FILE} ...")
    # GraphML cannot store None values; replace with empty strings
    for n, data in G.nodes(data=True):
        for k, v in list(data.items()):
            if v is None:
                data[k] = ''
    for u, v, k, data in G.edges(data=True, keys=True):
        for key, val in list(data.items()):
            if val is None:
                data[key] = ''

    nx.write_graphml(G, str(GRAPHML_FILE))
    print(f"  GraphML saved.")

    # Copy final node/edge CSVs to graph directory
    shutil.copy(NODES_CSV, FINAL_NODES_CSV)
    shutil.copy(EDGES_CSV, FINAL_EDGES_CSV)
    print(f"  Final CSVs  → {FINAL_NODES_CSV.parent}\n")

    return G


if __name__ == "__main__":
    build_kg()
