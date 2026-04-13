"""
DiabetesKG Pipeline — Shared Configuration
All paths, constants, seed genes, disease IDs, and column definitions live here.
"""
from pathlib import Path

# ── Base Directories ──────────────────────────────────────────────────────────
DATA_DIR   = Path(r"d:\ali_bioinfo")
OUTPUT_DIR = DATA_DIR / "output"

# ── Output Subdirectories ─────────────────────────────────────────────────────
EXTRACTED_DIR  = OUTPUT_DIR / "extracted"
INTEGRATED_DIR = OUTPUT_DIR / "integrated"
GRAPH_DIR      = OUTPUT_DIR / "graph"
VIZ_DIR        = OUTPUT_DIR / "visualizations"
VALIDATION_DIR = OUTPUT_DIR / "validation"

ALL_OUTPUT_DIRS = [OUTPUT_DIR, EXTRACTED_DIR, INTEGRATED_DIR, GRAPH_DIR, VIZ_DIR, VALIDATION_DIR]

# ── Raw Data Files (already in d:\ali_bioinfo) ───────────────────────────────
BIOGRID_FILE  = DATA_DIR / "BIOGRID-ALL-5.0.256.tab2.txt"
CHEBI_OWL     = DATA_DIR / "chebi_1.owl"
GO_OWL        = DATA_DIR / "go.owl"
GO_BASIC_OWL  = DATA_DIR / "go-basic.owl"
GOA_GAF       = DATA_DIR / "goa_human.gaf.gz"
HP_OBO        = DATA_DIR / "hp.obo"
GENE_INFO     = DATA_DIR / "gene_info.gz"
BIOKG_TAR     = DATA_DIR / "biokg-source.tar.gz"

# ── Downloaded Datasets (may not exist yet) ───────────────────────────────────
DISGENET_FILE  = DATA_DIR / "curated_gene_disease_associations.tsv"
REACTOME_FILE  = DATA_DIR / "UniProt2Reactome_All_Levels.txt"
CTD_FILE       = DATA_DIR / "CTD_chemicals_diseases.csv.gz"

# ── Extracted Intermediate Files ──────────────────────────────────────────────
GENE_INFO_CSV          = EXTRACTED_DIR / "gene_info_human.csv"
DIABETES_GENE_SET_FILE = EXTRACTED_DIR / "diabetes_gene_set.txt"
BIOGRID_CSV            = EXTRACTED_DIR / "biogrid_human_diabetes.csv"
GOA_CSV                = EXTRACTED_DIR / "goa_human_diabetes.csv"
GO_TERMS_CSV           = EXTRACTED_DIR / "go_terms_diabetes.csv"
GO_HIER_CSV            = EXTRACTED_DIR / "go_hierarchy_diabetes.csv"
CHEBI_CSV              = EXTRACTED_DIR / "chebi_diabetes.csv"
HP_TERMS_CSV           = EXTRACTED_DIR / "hp_diabetes_terms.csv"
HP_HIER_CSV            = EXTRACTED_DIR / "hp_diabetes_hierarchy.csv"
BIOKG_CSV              = EXTRACTED_DIR / "biokg_diabetes.csv"
BIOKG_CONTENTS_DIR     = EXTRACTED_DIR / "biokg_contents"
DISGENET_CSV           = EXTRACTED_DIR / "disgenet_diabetes.csv"
REACTOME_CSV           = EXTRACTED_DIR / "reactome_diabetes.csv"
CTD_CSV                = EXTRACTED_DIR / "ctd_chemicals_diabetes.csv"

# ── Integrated Files ──────────────────────────────────────────────────────────
GENE_MAP_CSV     = INTEGRATED_DIR / "gene_map.csv"
PROTEIN_MAP_CSV  = INTEGRATED_DIR / "protein_map.csv"
CHEMICAL_MAP_CSV = INTEGRATED_DIR / "chemical_map.csv"
DISEASE_MAP_CSV  = INTEGRATED_DIR / "disease_map.csv"
NODES_CSV        = INTEGRATED_DIR / "diabeteskg_nodes.csv"
EDGES_CSV        = INTEGRATED_DIR / "diabeteskg_edges.csv"

# ── Final Graph Output Files ──────────────────────────────────────────────────
GRAPHML_FILE    = GRAPH_DIR / "diabeteskg.graphml"
FINAL_NODES_CSV = GRAPH_DIR / "diabeteskg_nodes.csv"
FINAL_EDGES_CSV = GRAPH_DIR / "diabeteskg_edges.csv"

# ── Disease Focus: Diabetes Mellitus ─────────────────────────────────────────
HP_DIABETES_ROOT = "HP:0000819"

DIABETES_UMLS_IDS = {
    "C0011847",  # Diabetes Mellitus (general)
    "C0011854",  # Type 2 DM
    "C0011860",  # Type 1 DM
    "C0011882",  # Gestational DM
}

DIABETES_MESH_IDS = {
    "D003920",   # Diabetes Mellitus
    "D003924",   # Type 2 Diabetes Mellitus
}

# ── Seed Gene Set (always included regardless of BioGRID) ────────────────────
DIABETES_SEED_GENES = {
    "INS", "INSR", "GCK", "PPARG", "TCF7L2", "KCNJ11",
    "HNF1A", "HNF4A", "SLC2A2", "PDX1", "ABCC8",
    "IRS1", "IRS2", "PIK3R1", "FOXO1",
}

# Maximum gene set size after 1-hop BioGRID expansion
MAX_GENE_SET_SIZE = 2000

# ── ChEBI Diabetes Drug Seeds ─────────────────────────────────────────────────
CHEBI_DIABETES_DRUGS = {
    "CHEBI:5931":  "Insulin",
    "CHEBI:6801":  "Metformin",
    "CHEBI:5384":  "Glipizide",
    "CHEBI:50122": "Rosiglitazone",
    "CHEBI:10070": "Sitagliptin",
    "CHEBI:5391":  "Glucagon",
    "CHEBI:17234": "Glucose",
}

# ── Key GO Term Seeds for Diabetes ───────────────────────────────────────────
DIABETES_GO_SEEDS = {
    "GO:0006006",  # glucose metabolic process
    "GO:0008286",  # insulin receptor signaling pathway
    "GO:0005979",  # regulation of glycogen biosynthetic process
    "GO:0046323",  # glucose import
    "GO:0042593",  # glucose homeostasis
}

# ── Node and Edge Type Definitions ───────────────────────────────────────────
NODE_TYPES = {
    "Gene", "Protein", "BiologicalProcess", "MolecularFunction",
    "CellularComponent", "Phenotype", "Chemical", "Disease", "Pathway",
}

GO_NAMESPACE_TO_NODE_TYPE = {
    "biological_process": "BiologicalProcess",
    "molecular_function": "MolecularFunction",
    "cellular_component": "CellularComponent",
}

# Allowed (source_type, target_type) per edge type; None = validate per-row
SCHEMA = {
    "protein_has_biological_process":    ("Protein",   "BiologicalProcess"),
    "protein_has_molecular_function":    ("Protein",   "MolecularFunction"),
    "protein_located_in_component":      ("Protein",   "CellularComponent"),
    "protein_interacts_with_protein":    ("Protein",   "Protein"),
    "gene_encodes_protein":              ("Gene",      "Protein"),
    "gene_associated_with_disease":      ("Gene",      "Disease"),
    "protein_participates_in_pathway":   ("Protein",   "Pathway"),
    "chemical_associated_with_disease":  ("Chemical",  "Disease"),
    "phenotype_is_a":                    ("Phenotype", "Phenotype"),
    "go_term_is_a":                      None,  # same namespace source→target
    "go_term_part_of":                   None,  # same namespace source→target
    "chemical_is_a_subclass_of":         ("Chemical",  "Chemical"),
    "phenotype_associated_with_disease": ("Phenotype", "Disease"),
}

# ── Column Definitions ───────────────────────────────────────────────────────
GAF_COLS = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
    "DB_Reference", "Evidence_Code", "With_From", "Aspect", "Object_Name",
    "DB_Object_Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By",
    "Annotation_Extension", "Gene_Product_Form_ID",
]

CTD_CHEM_DISEASE_COLS = [
    "ChemicalName", "ChemicalID", "CasRN", "DiseaseName", "DiseaseID",
    "DirectEvidence", "InferenceGeneSymbol", "InferenceScore", "OmimIDs", "PubMedIDs",
]

REACTOME_COLS = [
    "UniProt_ID", "Reactome_pathway_ID", "URL", "pathway_name",
    "evidence_code", "species",
]

DISGENET_COLS = [
    "geneId", "geneSymbol", "DSI", "DPI", "diseaseId", "diseaseName",
    "diseaseType", "diseaseClass", "diseaseSemanticType", "score",
    "EI", "YI", "NofPmids", "NofSnps", "source",
]

# BioGRID chunked-read size (rows per chunk)
BIOGRID_CHUNK_SIZE = 100_000


def setup_dirs():
    """Create all output directories."""
    for d in ALL_OUTPUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
    print(f"  Output directories ready under {OUTPUT_DIR}")


if __name__ == "__main__":
    print("=== DiabetesKG Configuration ===")
    print(f"  DATA_DIR    : {DATA_DIR}")
    print(f"  OUTPUT_DIR  : {OUTPUT_DIR}")
    print(f"  Seed genes  : {len(DIABETES_SEED_GENES)}")
    print(f"  UMLS IDs    : {DIABETES_UMLS_IDS}")
    print(f"  MeSH IDs    : {DIABETES_MESH_IDS}")
    print(f"  ChEBI drugs : {len(CHEBI_DIABETES_DRUGS)}")
    print(f"  Node types  : {NODE_TYPES}")
    setup_dirs()
    print("Configuration OK.")
