# DiabetesKG: Biomedical Knowledge Graph Construction Pipeline (Midterm Paper) 

Live Demo: Bioinformatics Knowledge Grap  ([URL](https://bioinformatics-knowledge-graph.vercel.app/)).

This repository contains a modular, step-by-step pipeline designed to construct a comprehensive Biomedical Knowledge Graph (KG) focused on
Diabetes. The project integrates diverse data sources—including genomic, chemical, and ontological data—to create an interoperable graph 
structure for clinical reasoning and research.

## 📁 Repository Structure

## 1 Core Pipeline (`/pipeline`) 
The heart of the project. It contains the execution scripts for each stage of the KG construction.

- `run_pipeline.py`: The main entry point to execute the entire workflow.
- `config.py`: Centralized configuration for file paths and parameters.

Step-by-step Modules:
  - `step2_extract/`: Scripts for parsing raw biomedical datasets.
  - `step3_harmonize/`: Logic for entity alignment and data normalization.
  - `step4_nodes/` &  `step5_edges/`: Generators for graph components.
  - `step6_graph/`: Construction of the final GraphML and CSV structures.
  - `step7_viz/`: Automated visualization generation.
  - `step8_validate/`: Scripts for generating data integrity and quality reports.

## 2 Generated Outputs (`/output`)
Contains all data and visual artifacts produced by the pipeline.

- `extracted/`: Intermediary processed files like `chebi_diabetes.csv, go_terms_diabetes.csv` and gene sets.
- `graph/`: The final Knowledge Graph files, including `diabeteskg.graphml` for use in tools like Gephi or Cytoscape.
- `integrated/`: Harmonized node and edge lists (CSV) ready for database ingestion.
- `validation/`: Contains `validation_report.txt` which summarizes graph statistics.
- `visualizations/`: PNG and HTML exports showing degree distributions, schema diagrams, and network subgraphs.

## 3 Web Application (`/output/extracted/biokg_contents/biokg`)
A Next.js-based frontend application designed to interactively explore the generated Knowledge Graph ([URL](https://bioinformatics-knowledge-graph.vercel.app/)).

- Uses `kg_output.json` to render interactive network graphs.
- Built with TypeScript and Tailwind CSS for a modern, responsive research interface.


## 4 External Tools & Libraries (`/lib` & `/Protege-5.6.9`) 


- `lib/`: Client-side dependencies for graph rendering (e.g., `vis-9.1.2)`.
- `Protege-5.6.9/`: Local Protegé instance used for managing the Neuro Behavior Ontology (NBO) and clinical reasoning.


## 🚀 Getting Started
<b> Prerequisites </b> 

  Python 3.8+
  Node.js (for the web dashboard)
  Protegé (for ontology editing)


<b> Running the Pipeline </b>
To generate the knowledge graph from scratch, ensure you have Python 3.8+ installed, then run

```bash
python pipeline/run_pipeline.py
```

<b> Local Web Development <b>
To run the visualization dashboard locally:
  1. Navigate to the app directory:
     ```bash
     cd output/extracted/biokg_contents/biokg
     ```
  2. Install dependencies and start:
     ```bash
     npm install && npm run dev
     ```

