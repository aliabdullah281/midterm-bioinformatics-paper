"""
Step 8 — Validate DiabetesKG.

Checks:
  1. Schema validation   — every edge's node types match SCHEMA
  2. Referential integrity — all edge endpoints present as nodes
  3. Seed gene coverage    — which of 15 seed genes appear
  4. Connectivity          — connected components analysis
  5. Scale-free check      — log-log degree fit R²
  6. Spot checks           — specific known associations
  7. Duplicate detection   — duplicate node_ids or edges
  8. Summary statistics    — complete report

Output: output/validation/validation_report.txt
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
import networkx as nx
from config import (
    NODES_CSV, EDGES_CSV, GRAPHML_FILE, VALIDATION_DIR,
    DIABETES_SEED_GENES, SCHEMA, NODE_TYPES,
)


def _load_data():
    if not NODES_CSV.exists() or not EDGES_CSV.exists():
        print(f"  ERROR: nodes.csv or edges.csv not found. Run Steps 4-5 first.")
        return None, None
    nodes = pd.read_csv(NODES_CSV)
    edges = pd.read_csv(EDGES_CSV)
    return nodes, edges


def _load_graph():
    if not GRAPHML_FILE.exists():
        return None
    G = nx.read_graphml(str(GRAPHML_FILE))
    return G


def validate():
    print("\n[Step 8] Validating DiabetesKG")
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

    nodes, edges = _load_data()
    if nodes is None:
        return

    report_lines = ["=" * 70, "DiabetesKG Validation Report", "=" * 70, ""]

    def log(msg=""):
        print(f"  {msg}")
        report_lines.append(msg)

    # ── Overview ──────────────────────────────────────────────────────────────
    log(f"Nodes total  : {len(nodes):,}")
    log(f"Edges total  : {len(edges):,}")
    log()

    node_id_set = set(nodes['node_id'])

    # ── Check 1: Schema Validation ────────────────────────────────────────────
    log("─" * 60)
    log("CHECK 1: Schema Validation")
    log("─" * 60)

    # Build a lookup: node_id → node_type
    type_map = dict(zip(nodes['node_id'], nodes['node_type']))

    schema_violations = []
    schema_pass = 0
    for _, row in edges.iterrows():
        etype = row.get('edge_type', '')
        src_type = type_map.get(row.get('source_id', ''), 'UNKNOWN')
        tgt_type = type_map.get(row.get('target_id', ''), 'UNKNOWN')

        expected = SCHEMA.get(etype)
        if expected is None:
            # edge type not in schema (or None = flexible)
            continue
        exp_src, exp_tgt = expected
        if src_type != exp_src or tgt_type != exp_tgt:
            schema_violations.append({
                'edge_type': etype,
                'source_id': row['source_id'],
                'target_id': row['target_id'],
                'src_type_found': src_type,
                'tgt_type_found': tgt_type,
                'src_type_expected': exp_src,
                'tgt_type_expected': exp_tgt,
            })
        else:
            schema_pass += 1

    log(f"  Schema-conformant edges : {schema_pass:,}")
    log(f"  Schema violations       : {len(schema_violations):,}")
    if schema_violations:
        log("  First 5 violations:")
        for v in schema_violations[:5]:
            log(f"    {v['edge_type']}: src={v['src_type_found']}≠{v['src_type_expected']}, "
                f"tgt={v['tgt_type_found']}≠{v['tgt_type_expected']}")
    log()

    # ── Check 2: Referential Integrity ────────────────────────────────────────
    log("─" * 60)
    log("CHECK 2: Referential Integrity")
    log("─" * 60)

    dangling_src = edges[~edges['source_id'].isin(node_id_set)]
    dangling_tgt = edges[~edges['target_id'].isin(node_id_set)]
    log(f"  Edges with missing source node : {len(dangling_src):,}")
    log(f"  Edges with missing target node : {len(dangling_tgt):,}")
    if len(dangling_src) > 0:
        log(f"  Example missing source IDs : {dangling_src['source_id'].head(3).tolist()}")
    if len(dangling_tgt) > 0:
        log(f"  Example missing target IDs : {dangling_tgt['target_id'].head(3).tolist()}")
    log()

    # ── Check 3: Seed Gene Coverage ───────────────────────────────────────────
    log("─" * 60)
    log("CHECK 3: Diabetes Seed Gene Coverage")
    log("─" * 60)

    gene_nodes = nodes[nodes['node_type'] == 'Gene']
    gene_names_in_kg = set(gene_nodes['node_name'].dropna().str.upper())
    present  = {g for g in DIABETES_SEED_GENES if g.upper() in gene_names_in_kg}
    missing  = DIABETES_SEED_GENES - present
    log(f"  Seed genes expected  : {len(DIABETES_SEED_GENES)}")
    log(f"  Seed genes found     : {len(present)}  →  {sorted(present)}")
    log(f"  Seed genes missing   : {len(missing)}  →  {sorted(missing)}")
    log()

    # ── Check 4: Connectivity ─────────────────────────────────────────────────
    log("─" * 60)
    log("CHECK 4: Connectivity")
    log("─" * 60)

    G = _load_graph()
    if G is None:
        log("  GraphML not found — skipping connectivity checks.")
    else:
        G_und = G.to_undirected()
        comps = list(nx.connected_components(G_und))
        largest = max(comps, key=len)
        log(f"  Connected components : {len(comps):,}")
        log(f"  Largest component    : {len(largest):,} nodes "
            f"({100 * len(largest) / G_und.number_of_nodes():.1f}%)")
        log(f"  Isolated nodes       : {sum(1 for c in comps if len(c) == 1):,}")
    log()

    # ── Check 5: Scale-free Check ─────────────────────────────────────────────
    log("─" * 60)
    log("CHECK 5: Scale-free Degree Distribution (log-log R²)")
    log("─" * 60)

    if G is not None:
        import numpy as np
        from collections import Counter
        degrees = [d for _, d in G.degree() if d > 0]
        if degrees:
            cnt = Counter(degrees)
            x = sorted(cnt.keys())
            y = [cnt[k] for k in x]
            log_x = np.log10(x)
            log_y = np.log10(y)
            coef = np.polyfit(log_x, log_y, 1)
            y_pred = np.polyval(coef, log_x)
            ss_res = np.sum((log_y - y_pred) ** 2)
            ss_tot = np.sum((log_y - log_y.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
            log(f"  Power-law exponent γ : {-coef[0]:.3f}  (expected ≈ 2-3 for scale-free)")
            log(f"  Log-log R²           : {r2:.3f}  (≥0.80 suggests scale-free)")
    else:
        log("  GraphML not found — skipping.")
    log()

    # ── Check 6: Spot Checks ──────────────────────────────────────────────────
    log("─" * 60)
    log("CHECK 6: Known Association Spot Checks")
    log("─" * 60)

    def _check_association(src_name_substr, tgt_name_substr, etype_substr, label):
        hits = edges[
            edges['edge_type'].str.contains(etype_substr, na=False) &
            (
                edges['source_id'].str.contains(src_name_substr, case=False, na=False) |
                edges['source_id'].isin(
                    nodes[nodes['node_name'].str.contains(src_name_substr, case=False, na=False)]['node_id']
                )
            ) &
            (
                edges['target_id'].str.contains(tgt_name_substr, case=False, na=False) |
                edges['target_id'].isin(
                    nodes[nodes['node_name'].str.contains(tgt_name_substr, case=False, na=False)]['node_id']
                )
            )
        ]
        status = "PASS" if len(hits) > 0 else "FAIL"
        log(f"  [{status}] {label} ({len(hits)} edges found)")

    _check_association('CHEBI:6801', 'D003920', 'chemical',   "Metformin → Diabetes")
    _check_association('Metformin',  'iabetes',  'chemical',   "Metformin → Diabetes (by name)")
    _check_association('INS',        'disease',  'disease',    "INS gene → Disease")
    _check_association('PPARG',      'GO:',      'process',    "PPARG → GO:biological_process")
    _check_association('HP:',        'HP:',      'phenotype',  "HPO phenotype hierarchy (is_a)")
    log()

    # ── Check 7: Duplicate Detection ─────────────────────────────────────────
    log("─" * 60)
    log("CHECK 7: Duplicate Detection")
    log("─" * 60)

    dup_nodes = nodes[nodes.duplicated(subset=['node_id'])]['node_id']
    dup_edges = edges[edges.duplicated(subset=['source_id', 'target_id', 'edge_type'])]
    log(f"  Duplicate node_ids   : {len(dup_nodes):,}")
    log(f"  Duplicate (src, tgt, type) edges : {len(dup_edges):,}")
    if len(dup_nodes) > 0:
        log(f"  First duplicate node_ids : {dup_nodes.head(3).tolist()}")
    log()

    # ── Check 8: Summary Statistics ───────────────────────────────────────────
    log("─" * 60)
    log("CHECK 8: Full Summary Statistics")
    log("─" * 60)

    log("  Node type breakdown:")
    for nt, cnt in nodes['node_type'].value_counts().items():
        log(f"    {nt:<30} {cnt:>7,}")
    log()
    log("  Edge type breakdown:")
    for et, cnt in edges['edge_type'].value_counts().items():
        log(f"    {et:<40} {cnt:>7,}")
    log()

    if G is not None:
        log(f"  Graph density : {nx.density(G):.6f}")
        if G.number_of_nodes() > 0:
            avg_deg = sum(d for _, d in G.degree()) / G.number_of_nodes()
            log(f"  Average degree : {avg_deg:.2f}")
    log()

    # ── Write report ──────────────────────────────────────────────────────────
    report_path = VALIDATION_DIR / 'validation_report.txt'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    print(f"\n  Validation report saved → {report_path}\n")


if __name__ == "__main__":
    validate()
