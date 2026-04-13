"""
Step 7 — Visualizations for DiabetesKG.
Produces 7 outputs in output/visualizations/:
  1. schema_diagram.png        — node types + edge types
  2. diabetes_subgraph.png     — ego graph around INS/INSR/GCK
  3. degree_distribution.png   — log-log scatter
  4. node_type_distribution.png
  5. edge_type_distribution.png
  6. diabetes_network.html     — interactive pyvis (≤500 nodes)
  7. go_enrichment_heatmap.png — top GO terms × top diabetes proteins
"""
import sys, json, warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

warnings.filterwarnings('ignore')

import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from config import (
    GRAPHML_FILE, NODES_CSV, EDGES_CSV,
    VIZ_DIR, DIABETES_SEED_GENES,
)

# Colour palette per node type
NODE_COLORS = {
    'Gene':               '#E41A1C',
    'Protein':            '#FF7F00',
    'BiologicalProcess':  '#4DAF4A',
    'MolecularFunction':  '#984EA3',
    'CellularComponent':  '#377EB8',
    'Phenotype':          '#F781BF',
    'Chemical':           '#A65628',
    'Disease':            '#999999',
    'Pathway':            '#FFFF33',
}


def _load_graph():
    if not GRAPHML_FILE.exists():
        print(f"  ERROR: {GRAPHML_FILE} not found. Run Step 6 first.")
        return None
    print(f"  Loading graph from {GRAPHML_FILE} ...")
    G = nx.read_graphml(str(GRAPHML_FILE))
    print(f"  Nodes: {G.number_of_nodes():,}  Edges: {G.number_of_edges():,}")
    return G


# ── 1. Schema diagram ─────────────────────────────────────────────────────────
def plot_schema(nodes_csv, edges_csv, out_path):
    print("  [7.1] Schema diagram ...")
    nodes = pd.read_csv(nodes_csv)
    edges = pd.read_csv(edges_csv)

    type_counts = nodes['node_type'].value_counts().to_dict()
    edge_types  = edges['edge_type'].value_counts()

    from config import SCHEMA
    # Build a small schema graph
    Gs = nx.DiGraph()
    for nt, cnt in type_counts.items():
        Gs.add_node(nt, count=cnt)

    for etype, pair in SCHEMA.items():
        if pair:
            src_t, tgt_t = pair
            if src_t in Gs and tgt_t in Gs:
                lbl = etype.replace('_', '\n')
                Gs.add_edge(src_t, tgt_t, label=lbl)

    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_facecolor('#F8F9FA')
    ax.axis('off')

    pos = nx.spring_layout(Gs, seed=42, k=3.5)
    node_color_list = [NODE_COLORS.get(n, '#CCCCCC') for n in Gs.nodes()]
    node_sizes = [max(300, type_counts.get(n, 0) * 0.5 + 600) for n in Gs.nodes()]

    nx.draw_networkx_nodes(Gs, pos, node_color=node_color_list,
                           node_size=node_sizes, alpha=0.9, ax=ax)
    nx.draw_networkx_labels(Gs, pos, font_size=9, font_weight='bold', ax=ax)
    nx.draw_networkx_edges(Gs, pos, edge_color='#555555', arrows=True,
                           arrowsize=20, width=1.5,
                           connectionstyle='arc3,rad=0.1', ax=ax)

    # Legend
    legend = [mpatches.Patch(color=c, label=f"{t} ({type_counts.get(t,0):,})")
              for t, c in NODE_COLORS.items() if t in type_counts]
    ax.legend(handles=legend, loc='lower left', fontsize=8, framealpha=0.8)
    ax.set_title('DiabetesKG Schema: Node Types and Relationships', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {out_path}")


# ── 2. Diabetes core subgraph ─────────────────────────────────────────────────
def plot_diabetes_subgraph(G, out_path):
    print("  [7.2] Diabetes core subgraph (seed genes, depth=2) ...")

    # Find Gene/Protein nodes for seed genes
    seed_nodes = set()
    for n, d in G.nodes(data=True):
        name = d.get('node_name', '') or d.get('gene_symbol', '')
        if str(name).upper() in {g.upper() for g in ['INS', 'INSR', 'GCK', 'PPARG', 'TCF7L2']}:
            seed_nodes.add(n)

    if not seed_nodes:
        print("    WARNING: Seed gene nodes not found in graph — skipping subgraph plot.")
        return

    # 1-hop ego graph
    ego_nodes = set(seed_nodes)
    for sn in seed_nodes:
        ego_nodes |= set(G.successors(sn)) | set(G.predecessors(sn))
    # Limit size
    ego_nodes = list(ego_nodes)[:200]

    Gsub = G.subgraph(ego_nodes).copy().to_undirected()
    print(f"    Subgraph: {Gsub.number_of_nodes()} nodes, {Gsub.number_of_edges()} edges")

    fig, ax = plt.subplots(figsize=(16, 14))
    ax.set_facecolor('#F8F9FA')
    ax.axis('off')

    pos = nx.spring_layout(Gsub, seed=42, k=1.5)
    node_colors = [NODE_COLORS.get(Gsub.nodes[n].get('node_type', ''), '#CCCCCC') for n in Gsub.nodes()]
    is_seed = [n in seed_nodes for n in Gsub.nodes()]
    node_sizes = [800 if s else 200 for s in is_seed]

    nx.draw_networkx_nodes(Gsub, pos, node_color=node_colors,
                           node_size=node_sizes, alpha=0.85, ax=ax)
    seed_labels = {n: Gsub.nodes[n].get('node_name', n) for n in Gsub.nodes() if n in seed_nodes}
    nx.draw_networkx_labels(Gsub, pos, labels=seed_labels, font_size=8, font_weight='bold', ax=ax)
    nx.draw_networkx_edges(Gsub, pos, edge_color='#AAAAAA', alpha=0.5, width=0.7, ax=ax)

    legend = [mpatches.Patch(color=c, label=t) for t, c in NODE_COLORS.items()]
    ax.legend(handles=legend, loc='lower left', fontsize=7, framealpha=0.8, ncol=2)
    ax.set_title('DiabetesKG: Core Diabetes Gene Neighbourhood (INS, INSR, GCK, PPARG, TCF7L2)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {out_path}")


# ── 3. Degree distribution ────────────────────────────────────────────────────
def plot_degree_distribution(G, out_path):
    print("  [7.3] Degree distribution (log-log) ...")
    import numpy as np
    degrees = sorted([d for _, d in G.degree()], reverse=True)
    if not degrees:
        return
    from collections import Counter
    counts = Counter(degrees)
    x = sorted(counts.keys())
    y = [counts[k] for k in x]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, s=15, alpha=0.6, color='steelblue')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Degree (k)', fontsize=12)
    ax.set_ylabel('Number of nodes P(k)', fontsize=12)
    ax.set_title('DiabetesKG Degree Distribution (log-log)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')

    # Fit line
    log_x = np.log10([xi for xi in x if xi > 0])
    log_y = np.log10([yi for yi, xi in zip(y, x) if xi > 0])
    if len(log_x) > 5:
        coef = np.polyfit(log_x, log_y, 1)
        fit_x = np.logspace(min(log_x), max(log_x), 50)
        fit_y = 10 ** np.polyval(coef, np.log10(fit_x))
        ax.plot(fit_x, fit_y, 'r--', alpha=0.7, label=f'Fit: γ={-coef[0]:.2f}')
        ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {out_path}")


# ── 4 & 5. Node/Edge type bar charts ─────────────────────────────────────────
def plot_type_distributions(nodes_csv, edges_csv, out_nodes, out_edges):
    print("  [7.4] Node type distribution ...")
    nodes = pd.read_csv(nodes_csv)
    edges = pd.read_csv(edges_csv)

    # Node types
    nt = nodes['node_type'].value_counts()
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = [NODE_COLORS.get(t, '#CCCCCC') for t in nt.index]
    nt.plot(kind='bar', color=colors, edgecolor='black', ax=ax)
    ax.set_title('DiabetesKG: Node Count by Type', fontsize=13, fontweight='bold')
    ax.set_xlabel('Node Type', fontsize=11); ax.set_ylabel('Count', fontsize=11)
    ax.tick_params(axis='x', rotation=30)
    for i, v in enumerate(nt.values):
        ax.text(i, v + 5, f'{v:,}', ha='center', va='bottom', fontsize=8)
    plt.tight_layout()
    plt.savefig(out_nodes, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {out_nodes}")

    # Edge types
    print("  [7.5] Edge type distribution ...")
    et = edges['edge_type'].value_counts()
    fig, ax = plt.subplots(figsize=(14, 6))
    et.plot(kind='bar', color='steelblue', edgecolor='black', ax=ax)
    ax.set_title('DiabetesKG: Edge Count by Type', fontsize=13, fontweight='bold')
    ax.set_xlabel('Edge Type', fontsize=11); ax.set_ylabel('Count', fontsize=11)
    ax.tick_params(axis='x', rotation=40)
    for i, v in enumerate(et.values):
        ax.text(i, v + 2, f'{v:,}', ha='center', va='bottom', fontsize=7)
    plt.tight_layout()
    plt.savefig(out_edges, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {out_edges}")


# ── 6. Interactive pyvis HTML ─────────────────────────────────────────────────
def plot_interactive_html(G, out_path, max_nodes=500):
    print(f"  [7.6] Interactive HTML graph (max {max_nodes} nodes) ...")
    try:
        from pyvis.network import Network
    except ImportError:
        print("    WARNING: pyvis not installed. Run: pip install pyvis")
        return

    # Subgraph: diabetes seed genes + 1 hop, capped at max_nodes
    seed_nodes = set()
    for n, d in G.nodes(data=True):
        nt = d.get('node_type', '')
        nm = d.get('node_name', '') or d.get('gene_symbol', '')
        if nt == 'Gene' and str(nm).upper() in {g.upper() for g in DIABETES_SEED_GENES}:
            seed_nodes.add(n)

    viz_nodes = set(seed_nodes)
    for sn in list(seed_nodes):
        viz_nodes |= set(list(G.successors(sn))[:20]) | set(list(G.predecessors(sn))[:20])
    viz_nodes = list(viz_nodes)[:max_nodes]

    Gsub = G.subgraph(viz_nodes)

    net = Network(height='750px', width='100%', bgcolor='#222222', font_color='white',
                  directed=True)
    net.set_options("""
    var options = {
      "nodes": {"shape": "dot", "scaling": {"min": 10, "max": 30}},
      "edges": {"arrows": {"to": {"enabled": true, "scaleFactor": 0.5}}, "smooth": true},
      "physics": {"enabled": true, "stabilization": {"iterations": 100}}
    }""")

    for n, d in Gsub.nodes(data=True):
        nt = d.get('node_type', 'Unknown')
        color = NODE_COLORS.get(nt, '#CCCCCC')
        label = d.get('node_name', n)[:30]
        title = f"{nt}: {d.get('node_name', n)}"
        size = 20 if n in seed_nodes else 10
        net.add_node(n, label=label, color=color, title=title, size=size)

    for u, v, d in Gsub.edges(data=True):
        net.add_edge(u, v, title=d.get('edge_type', ''), width=1)

    net.write_html(str(out_path))
    print(f"    Saved → {out_path}")


# ── 7. GO Term Enrichment Heatmap ─────────────────────────────────────────────
def plot_go_heatmap(nodes_csv, edges_csv, out_path):
    print("  [7.7] GO term enrichment heatmap ...")
    nodes = pd.read_csv(nodes_csv)
    edges = pd.read_csv(edges_csv)

    # Get protein→GO edges
    go_edges = edges[edges['edge_type'].isin([
        'protein_has_biological_process',
        'protein_has_molecular_function',
        'protein_located_in_component',
    ])].copy()

    if go_edges.empty:
        print("    No GO annotation edges found — skipping heatmap.")
        return

    # Get GO node names
    go_nodes = nodes[nodes['node_type'].isin(
        ['BiologicalProcess', 'MolecularFunction', 'CellularComponent']
    )][['node_id', 'node_name']].copy()
    go_name_map = dict(zip(go_nodes['node_id'], go_nodes['node_name']))

    # Get protein node names
    prot_nodes = nodes[nodes['node_type'] == 'Protein'][['node_id', 'node_name']].copy()
    prot_name_map = dict(zip(prot_nodes['node_id'], prot_nodes['node_name']))

    go_edges['protein_name'] = go_edges['source_id'].map(prot_name_map).fillna(go_edges['source_id'])
    go_edges['go_name']      = go_edges['target_id'].map(go_name_map).fillna(go_edges['target_id'])

    # Top 10 proteins by GO term count; top 20 GO terms by protein count
    top_prots = go_edges['protein_name'].value_counts().head(10).index.tolist()
    top_gos   = go_edges['go_name'].value_counts().head(20).index.tolist()

    sub = go_edges[
        go_edges['protein_name'].isin(top_prots) &
        go_edges['go_name'].isin(top_gos)
    ]
    if sub.empty:
        print("    Insufficient data for heatmap — skipping.")
        return

    matrix = sub.pivot_table(index='protein_name', columns='go_name',
                              values='edge_type', aggfunc='count', fill_value=0)

    fig, ax = plt.subplots(figsize=(max(14, len(top_gos)), max(6, len(top_prots))))
    sns.heatmap(matrix, cmap='Blues', linewidths=0.5, linecolor='grey',
                annot=True, fmt='d', ax=ax, cbar_kws={'label': 'Association Count'})
    ax.set_title('DiabetesKG: GO Term Enrichment (Top Diabetes Proteins × GO Terms)',
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('GO Term', fontsize=10)
    ax.set_ylabel('Protein', fontsize=10)
    ax.tick_params(axis='x', rotation=45, labelsize=7)
    ax.tick_params(axis='y', rotation=0, labelsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=100, bbox_inches='tight')
    plt.close()
    print(f"    Saved → {out_path}")


def visualize():
    print("\n[Step 7] Generating Visualizations")
    VIZ_DIR.mkdir(parents=True, exist_ok=True)

    if not NODES_CSV.exists() or not EDGES_CSV.exists():
        print("  ERROR: nodes.csv or edges.csv missing. Run Steps 4-5 first.")
        return

    G = _load_graph()

    plot_schema(NODES_CSV, EDGES_CSV, VIZ_DIR / 'schema_diagram.png')

    if G is not None:
        plot_diabetes_subgraph(G, VIZ_DIR / 'diabetes_subgraph.png')
        plot_degree_distribution(G, VIZ_DIR / 'degree_distribution.png')
        plot_interactive_html(G, VIZ_DIR / 'diabetes_network.html')

    plot_type_distributions(
        NODES_CSV, EDGES_CSV,
        VIZ_DIR / 'node_type_distribution.png',
        VIZ_DIR / 'edge_type_distribution.png',
    )
    plot_go_heatmap(NODES_CSV, EDGES_CSV, VIZ_DIR / 'go_enrichment_heatmap.png')

    print(f"\n  All visualizations saved to {VIZ_DIR}\n")


if __name__ == "__main__":
    visualize()
