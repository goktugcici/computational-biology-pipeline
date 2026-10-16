# Computational Biology & Multi-Omics Science Pipeline 🧬🔬

[![CI Pipeline](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml)
[![Pages Deployment](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/deploy-pages.yml/badge.svg)](https://goktugcici.github.io/computational-biology-pipeline/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg)](https://www.docker.com/)

A production-grade, modular, and reproducible computational biology platform integrating **3D Structural Biology**, **Genomic Variant Annotation**, **Pan-Cancer Oncology Profiling**, **Pharmacological Drug Discovery**, **Biological Pathway Mapping**, **Biomedical Literature Intelligence**, and an **Interactive WebGL Dashboard**.

---

## 🏛️ Architecture Overview

```
                      +-----------------------------------------+
                      |     Computational Biology Platform      |
                      +--------------------+--------------------+
                                           |
       +--------------------+--------------+--------------+--------------------+
       |                    |                             |                    |
       v                    v                             v                    v
+--------------+    +---------------+             +---------------+    +---------------+
| 3D Structure |    | Genomic State |             |  Pan-Cancer   |    |  Literature   |
| (PDB/Mol*)   |    | (ClinVar/SNP) |             | (EGFR/KRAS)   |    | (PubMed/PMC)  |
+--------------+    +---------------+             +---------------+    +---------------+
       |                    |                             |                    |
       | 4HHB Structure     | HBB rs334                   | EGFR L858R/T790M   | CRISPR Exa-cel
       | Active Histidines  | Pathogenic Missense         | KRAS G12C Sotorasib| NEJM Landmark
       v                    v                             v                    v
       +--------------------+--------------+--------------+--------------------+
                                           |
                                           v
       +-----------------------------------+-----------------------------------+
       |                                                                       |
       v                                                                       v
+-------------------------------+                       +-------------------------------+
|  Interactive WebGL Dashboard  |                       |   Executive HTML Reports      |
|  (3Dmol.js / GitHub Pages)    |                       |   (Multi-Omics Charts / PDF)  |
+-------------------------------+                       +-------------------------------+
```

---

## 🌟 Core Modules

### 1. 🔬 3D Structural Biology (`src/structural_biology.py`)
- Analyzes macromolecular complexes and PDB/mmCIF coordinates (`4HHB` - Human Deoxyhemoglobin).
- Identifies critical coordination sites:
  - **Proximal Histidine (`Chain A: His87`)**: Covalent coordination to heme iron ($Fe^{2+}$).
  - **Distal Histidine (`Chain A: His58`)**: Oxygen binding stabilization and steric resistance to CO.
  - **Sickle Cell Mutation Site (`Chain B/D: Glu6`)**: Glutamic acid to Valine substitution ($E6V$) causing pathological polymerization.

### 2. 🎯 Pan-Cancer Precision Oncology (`src/oncology_profiler.py`)
- Curated knowledgebase for high-impact oncology targets:
  - **EGFR**: Hotspots `p.L858R`, `p.T790M`, `Exon 19 del` $\rightarrow$ Targeted by *Osimertinib*, *Gefitinib*.
  - **KRAS**: Hotspots `p.G12C`, `p.G12D`, `p.G12V` $\rightarrow$ Targeted by *Sotorasib*, *Adagrasib*.
  - **BRAF**: Hotspots `p.V600E`, `p.V600K` $\rightarrow$ Targeted by *Dabrafenib + Trametinib*.
  - **TP53**: Hotspots `p.R175H`, `p.R248Q`, `p.R273H` $\rightarrow$ WEE1/ATR synthetic lethality.

### 3. ⚡ Dynamic Live Search Engine (`src/live_fetcher.py`)
- Real-time querying against live public APIs with an automatic local caching layer in `data/cache/`:
  - **RCSB PDB REST API**: Structural metadata and resolution.
  - **PubChem PUG-REST**: Real-time molecular weights, formulas, and lipophilicity.
  - **PubMed E-utilities**: Instant article search and abstract resolution.

### 4. 🧬 Genomic Variant & Clinical Significance (`src/variant_analyzer.py`)
- Integrates with **ClinVar** and **dbSNP** (`rs334`, `chr11:5227002:T>A`).
- Screens clinical interpretations:
  - `HBB:c.20A>T (p.Glu7Val / Glu6Val)` classified as **Pathogenic**.
  - Negative control comparison against `BRCA2:c.198A>G (rs28897700)` classified as **Benign** by expert panels.

### 5. 💊 Pharmacology & Cheminformatics (`src/pharmacology.py`)
- Retrieves chemical properties, SMILES structures, and mechanisms of action from **PubChem** and **ChEMBL**:
  - **Voxelotor (Oxbryta)** (PubChem CID: `71602803`): Binds reversibly to the N-terminal valine of hemoglobin's $\alpha$-chain to allosterically increase oxygen affinity, stabilizing non-polymerizing oxy-HbS.
  - **Metformin** (`CHEMBL1431`): Standard profiling for metabolic regulators.

### 6. 📊 Biological Pathway Mapping (`src/pathway_mapper.py`)
- Maps human genes to curated **Reactome** pathways:
  - `R-HSA-1247673`: *Erythrocytes take up oxygen and release carbon dioxide*
  - `R-HSA-1480926`: *O2/CO2 exchange in erythrocytes*
  - `R-HSA-9927020`: *Heme assimilation*

### 7. 📚 Scientific Literature Mining (`src/literature_miner.py`)
- Searches **PubMed / PMC** for clinical trial breakthroughs:
  - Landmark Phase 3 trial on **Exagamglogene Autotemcel (Casgevy)** (*NEJM*, PMID: `38661449`).
  - Demonstrates CRISPR-Cas9 editing of the *BCL11A* erythroid enhancer to reactivate fetal hemoglobin (HbF).

### 8. 📑 Publication-Grade Report Generator (`src/report_generator.py`)
- Generates a standalone, responsive, interactive HTML scientific report (`reports/executive_summary.html`) with embedded visual charts, executive metrics, and multi-omics findings.

---

## 🚀 Quick Start

### 1. Installation
Clone the repository and install dependencies:
```bash
git clone https://github.com/goktugcici/computational-biology-pipeline.git
cd computational-biology-pipeline
pip install -r requirements.txt
```

### 2. Run the Full Multi-Omics Pipeline
Execute the standard pipeline:
```bash
python src/main.py --data_dir data
```

### 3. Advanced CLI Operations
```bash
# Profile an oncogene and generate publication HTML report
python src/main.py --pan_cancer EGFR --generate_report

# Query live biological APIs on the fly
python src/main.py --live_pdb 1TSR --live_drug Imatinib --live_pubmed "CAR-T cell therapy"
```

### 4. Interactive Web Dashboard
Open `dashboard/index.html` directly in any web browser to access the 3D WebGL protein visualizer and oncology consoles, or view it live via GitHub Pages:
👉 **Live Web Dashboard:** [https://goktugcici.github.io/computational-biology-pipeline/](https://goktugcici.github.io/computational-biology-pipeline/)

---

## 🐳 Docker Deployment

Run the complete pipeline with zero local dependencies using Docker:
```bash
# Build and run with Docker Compose
docker-compose up --build
```

---

## 🧪 Testing Suite

Run the automated test matrix:
```bash
python -m unittest discover -s tests
```

---

## 📂 Project Structure

```
computational-biology-pipeline/
├── src/
│   ├── main.py                 # Master pipeline orchestrator & CLI
│   ├── structural_biology.py   # PDB/CIF coordinate parser
│   ├── oncology_profiler.py    # Pan-cancer hotspot & targeted therapy profiler
│   ├── live_fetcher.py         # Dynamic REST API search engine with cache
│   ├── variant_analyzer.py     # ClinVar/dbSNP variant classifier
│   ├── pharmacology.py         # PubChem/ChEMBL chemical profiler
│   ├── pathway_mapper.py       # Reactome biological pathway mapper
│   ├── literature_miner.py     # PubMed literature intelligence
│   └── report_generator.py     # Interactive HTML report generator
├── data/                       # Curated scientific data files
│   ├── cache/                  # Persistent cache for live API queries
│   ├── pdb_00004hhb.cif        # Human Deoxyhemoglobin coordinates
│   ├── summary_all.json        # ClinVar variant classifications
│   └── ...
├── dashboard/
│   └── index.html              # Standalone WebGL 3D molecular dashboard
├── reports/
│   ├── executive_summary.html  # Generated publication-grade HTML report
│   └── scientific_investigation_report.md
├── tests/
│   └── test_pipeline.py        # Automated unit test suite (100% pass)
├── .github/workflows/
│   ├── ci.yml                  # Automated test matrix (Python 3.9 - 3.12)
│   └── deploy-pages.yml        # Automatic GitHub Pages deployment
├── Dockerfile                  # Container definition
├── docker-compose.yml          # Container orchestration
├── requirements.txt            # Python dependencies
├── .gitignore                  # Git ignore rules
└── README.md                   # Platform documentation
```

---

## 📄 License
This project is open-source under the MIT License.
