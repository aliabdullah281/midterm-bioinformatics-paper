"""
DiabetesKG Master Pipeline Runner
===================================
Usage:
    cd d:\\ali_bioinfo
    python pipeline/run_pipeline.py                # run all steps
    python pipeline/run_pipeline.py --from-step 4  # start from step 4
    python pipeline/run_pipeline.py --from-step 2 --to-step 5

Steps:
    1  — Setup directories
    2  — Extract raw data (10 sub-scripts)
    3  — Harmonise IDs
    4  — Build nodes
    5  — Build edges
    6  — Assemble knowledge graph (NetworkX + GraphML)
    7  — Visualizations
    8  — Validate
"""
import sys, argparse, time, traceback
from pathlib import Path

# Allow imports from pipeline/ and pipeline/step*/
PIPELINE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PIPELINE_DIR))

from config import setup_dirs

# ── Step 2: extract sub-modules ───────────────────────────────────────────────
from step2_extract.extract_gene_info          import extract as extract_gene_info
from step2_extract.extract_biogrid            import extract as extract_biogrid
from step2_extract.extract_goa_human          import extract as extract_goa_human
from step2_extract.extract_gene_ontology      import extract as extract_gene_ontology
from step2_extract.extract_chebi              import extract as extract_chebi
from step2_extract.extract_hp_obo             import extract as extract_hp_obo
from step2_extract.extract_biokg              import extract as extract_biokg
from step2_extract.extract_disgenet           import extract as extract_disgenet
from step2_extract.extract_reactome           import extract as extract_reactome
from step2_extract.extract_ctd_chem_disease   import extract as extract_ctd

from step3_harmonize.harmonize_ids            import harmonize
from step4_nodes.build_nodes                  import build_nodes
from step5_edges.build_edges                  import build_edges
from step6_graph.build_knowledge_graph        import build_kg
from step7_viz.visualize_kg                   import visualize
from step8_validate.validate_kg               import validate


def _run_step(name, fn, errors, skip=False):
    if skip:
        print(f"\n  [SKIP] {name}")
        return True
    print(f"\n{'=' * 60}")
    print(f"  RUNNING: {name}")
    print(f"{'=' * 60}")
    t0 = time.time()
    try:
        fn()
        elapsed = time.time() - t0
        print(f"  ✔ {name} completed in {elapsed:.1f}s")
        return True
    except KeyboardInterrupt:
        raise
    except Exception as e:
        elapsed = time.time() - t0
        print(f"\n  ✖ {name} FAILED after {elapsed:.1f}s")
        print(f"  Error: {e}")
        traceback.print_exc()
        errors.append(name)
        return False


def run_pipeline(from_step: int = 1, to_step: int = 8, skip_on_error: bool = True):
    t_start = time.time()
    errors = []

    def should_run(step_num: int) -> bool:
        return from_step <= step_num <= to_step

    # ── Step 1: Setup ─────────────────────────────────────────────────────────
    if should_run(1):
        print("\n[Step 1] Setting up output directories ...")
        setup_dirs()
        print("  ✔ Directories ready")

    # ── Step 2: Extraction ────────────────────────────────────────────────────
    if should_run(2):
        print("\n[Step 2] Extracting raw data sources ...")
        # gene_info must run first (provides gene map for subsequent steps)
        _run_step("2.1 extract_gene_info",        extract_gene_info,  errors)
        # BioGRID builds the expanded gene set needed by GOA extraction
        _run_step("2.2 extract_biogrid",           extract_biogrid,    errors)
        # GOA needs the gene set;  GO OWL needs the GOA GO IDs
        _run_step("2.3 extract_goa_human",         extract_goa_human,  errors)
        _run_step("2.4 extract_gene_ontology",     extract_gene_ontology, errors)
        # These are independent
        _run_step("2.5 extract_chebi",             extract_chebi,      errors)
        _run_step("2.6 extract_hp_obo",            extract_hp_obo,     errors)
        _run_step("2.7 extract_biokg",             extract_biokg,      errors)
        # Optional downloads (graceful if file not present)
        _run_step("2.8 extract_disgenet",          extract_disgenet,   errors)
        _run_step("2.9 extract_reactome",          extract_reactome,   errors)
        _run_step("2.10 extract_ctd_chem_disease", extract_ctd,        errors)

    # ── Step 3: Harmonise ─────────────────────────────────────────────────────
    if should_run(3):
        _run_step("3. harmonize_ids", harmonize, errors)

    # ── Step 4: Nodes ─────────────────────────────────────────────────────────
    if should_run(4):
        _run_step("4. build_nodes", build_nodes, errors)

    # ── Step 5: Edges ─────────────────────────────────────────────────────────
    if should_run(5):
        _run_step("5. build_edges", build_edges, errors)

    # ── Step 6: Graph ─────────────────────────────────────────────────────────
    if should_run(6):
        _run_step("6. build_knowledge_graph", build_kg, errors)

    # ── Step 7: Visualize ─────────────────────────────────────────────────────
    if should_run(7):
        _run_step("7. visualize_kg", visualize, errors)

    # ── Step 8: Validate ──────────────────────────────────────────────────────
    if should_run(8):
        _run_step("8. validate_kg", validate, errors)

    # ── Summary ───────────────────────────────────────────────────────────────
    total = time.time() - t_start
    print("\n" + "=" * 60)
    print(f"  DiabetesKG Pipeline Complete  ({total:.1f}s total)")
    if errors:
        print(f"  Steps with errors : {errors}")
        print("  Check logs above for details. Outputs may be partial.")
    else:
        print("  All steps succeeded. Check output/ for results.")
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='DiabetesKG master pipeline runner',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('--from-step', type=int, default=1, dest='from_step',
                        help='Start from this step number (1-8). Default: 1')
    parser.add_argument('--to-step',   type=int, default=8, dest='to_step',
                        help='End at this step number (1-8).   Default: 8')
    parser.add_argument('--list-steps', action='store_true',
                        help='List all steps and exit.')
    args = parser.parse_args()

    if args.list_steps:
        print("""
Steps:
  1  Setup directories
  2  Extract raw data  (gene_info, BioGRID, GOA, GO OWL,
                         ChEBI, HP OBO, BioKG, DisGeNET,
                         Reactome, CTD)
  3  Harmonise IDs     (gene/protein/chemical/disease maps)
  4  Build nodes       (diabeteskg_nodes.csv)
  5  Build edges       (diabeteskg_edges.csv)
  6  Assemble KG       (GraphML + summary stats)
  7  Visualise         (7 plots + interactive HTML)
  8  Validate          (8 checks + validation_report.txt)
""")
        return

    run_pipeline(from_step=args.from_step, to_step=args.to_step)


if __name__ == "__main__":
    main()
