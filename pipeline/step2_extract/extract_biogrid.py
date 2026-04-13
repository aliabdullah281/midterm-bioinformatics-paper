"""
Step 2.1 — Extract BioGRID human physical interactions involving diabetes genes.
Strategy: chunked read (100k rows/chunk), filter human-human + physical + seed genes,
then expand 1-hop to collect all interaction partners.
Output: output/extracted/biogrid_human_diabetes.csv
         output/extracted/diabetes_gene_set.txt
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd
from config import (
    BIOGRID_FILE, BIOGRID_CSV, DIABETES_GENE_SET_FILE,
    DIABETES_SEED_GENES, BIOGRID_CHUNK_SIZE, MAX_GENE_SET_SIZE, EXTRACTED_DIR,
)


def extract():
    print("\n[Step 2.1] Extracting BioGRID human physical interactions (diabetes seed genes)")
    print(f"  Source: {BIOGRID_FILE}")

    if not BIOGRID_FILE.exists():
        print(f"  ERROR: {BIOGRID_FILE} not found. Skipping.")
        # Still save the seed gene set so downstream steps can proceed
        _save_gene_set(DIABETES_SEED_GENES)
        return None

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    kept_chunks = []
    total_rows = 0
    kept_rows = 0
    header_cols = None

    reader = pd.read_csv(
        BIOGRID_FILE,
        sep='\t',
        header=0,
        low_memory=False,
        chunksize=BIOGRID_CHUNK_SIZE,
        dtype=str,
    )

    for chunk_num, chunk in enumerate(reader):
        # Normalise column names (may have BOM '#' prefix on first column)
        chunk.columns = [c.strip().lstrip('\ufeff').lstrip('#').strip() for c in chunk.columns]
        if header_cols is None:
            header_cols = list(chunk.columns)

        total_rows += len(chunk)

        # Require these columns to exist
        req = ['Organism ID Interactor A', 'Organism ID Interactor B',
               'Experimental System Type', 'Official Symbol Interactor A',
               'Official Symbol Interactor B']
        if not all(c in chunk.columns for c in req):
            if chunk_num == 0:
                print(f"  WARNING: Expected columns not found. Found: {list(chunk.columns)[:8]}")
            continue

        # Filter: human-human
        mask_human = (
            (chunk['Organism ID Interactor A'].str.strip() == '9606') &
            (chunk['Organism ID Interactor B'].str.strip() == '9606')
        )
        # Filter: physical interactions only
        mask_phys = (
            chunk['Experimental System Type'].str.strip().str.lower() == 'physical'
        )
        # Filter: at least one interactor is a diabetes seed gene
        sym_a = chunk['Official Symbol Interactor A'].str.strip()
        sym_b = chunk['Official Symbol Interactor B'].str.strip()
        mask_seed = sym_a.isin(DIABETES_SEED_GENES) | sym_b.isin(DIABETES_SEED_GENES)

        filtered = chunk[mask_human & mask_phys & mask_seed].copy()
        kept_rows += len(filtered)
        kept_chunks.append(filtered)

        if (chunk_num + 1) % 5 == 0:
            print(f"    Processed {total_rows:,} rows, kept {kept_rows:,} ...")

    print(f"  Total rows scanned                : {total_rows:,}")
    print(f"  Rows with diabetes seed (physical): {kept_rows:,}")

    if kept_rows == 0:
        print("  WARNING: No interactions found. Saving seed genes only.")
        _save_gene_set(DIABETES_SEED_GENES)
        return None

    df = pd.concat(kept_chunks, ignore_index=True)

    # Keep useful columns
    keep = [
        'Official Symbol Interactor A', 'Official Symbol Interactor B',
        'Experimental System', 'Experimental System Type',
        'Throughput', 'Score', 'Source Database',
    ]
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()
    df.rename(columns={
        'Official Symbol Interactor A': 'symbol_a',
        'Official Symbol Interactor B': 'symbol_b',
        'Experimental System':          'exp_system',
        'Experimental System Type':     'exp_type',
        'Source Database':              'source_db',
    }, inplace=True)

    df.dropna(subset=['symbol_a', 'symbol_b'], inplace=True)
    df.drop_duplicates(subset=['symbol_a', 'symbol_b'], inplace=True)

    # Build expanded gene set = seeds ∪ 1-hop partners
    all_syms = set(df['symbol_a'].tolist()) | set(df['symbol_b'].tolist())
    diabetes_gene_set = DIABETES_SEED_GENES | all_syms

    if len(diabetes_gene_set) > MAX_GENE_SET_SIZE:
        print(f"  Gene set ({len(diabetes_gene_set):,}) > cap ({MAX_GENE_SET_SIZE}). Trimming...")
        counts = (
            pd.concat([df['symbol_a'], df['symbol_b']])
            .value_counts()
        )
        non_seeds = [g for g in counts.index if g not in DIABETES_SEED_GENES]
        top_partners = set(non_seeds[: MAX_GENE_SET_SIZE - len(DIABETES_SEED_GENES)])
        diabetes_gene_set = DIABETES_SEED_GENES | top_partners
        df = df[
            df['symbol_a'].isin(diabetes_gene_set) &
            df['symbol_b'].isin(diabetes_gene_set)
        ].copy()

    print(f"  Expanded diabetes gene set size   : {len(diabetes_gene_set):,}")
    print(f"  Final interaction rows            : {len(df):,}")

    _save_gene_set(diabetes_gene_set)
    df.to_csv(BIOGRID_CSV, index=False)
    print(f"  Saved interactions → {BIOGRID_CSV}")
    print(f"  Sample:\n{df.head(3).to_string()}\n")
    return df


def _save_gene_set(gene_set):
    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)
    with open(DIABETES_GENE_SET_FILE, 'w') as f:
        f.write('\n'.join(sorted(gene_set)))
    print(f"  Saved gene set ({len(gene_set):,} genes) → {DIABETES_GENE_SET_FILE}")


if __name__ == "__main__":
    extract()
