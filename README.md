# Computational Biology & Multi-Omics Science Pipeline 🧬🔬

[![CI Pipeline](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml)
[![Pages Deployment](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/deploy-pages.yml/badge.svg)](https://goktugcici.github.io/computational-biology-pipeline/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version: 2.0 Enterprise](https://img.shields.io/badge/Release-v2.0_Enterprise-emerald.svg)](https://github.com/goktugcici/computational-biology-pipeline/releases)
[![Tests Passing](https://img.shields.io/badge/Tests-11%2F11%20Passed-brightgreen.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions)

A production-grade, modular, and reproducible computational biology platform integrating **3D Structural Biology**, **Genomic Variant Annotation**, **Pan-Cancer Precision Oncology**, **Needleman-Wunsch Sequence Alignment**, **In Silico Drug Screening (Lipinski Ro5 & Binding Scoring)**, **Protein-Protein Interaction (PPI) Interactome Networks**, **Biological Pathways**, **Literature Intelligence**, and an **Interactive WebGL Dashboard**.

---

## 🏛️ Multi-Omics Architecture (v2.0)

```
                            +-----------------------------------------------+
                            |   Computational Biology Platform v2.0         |
                            +-----------------------+-----------------------+
                                                    |
       +--------------------+-----------------------+-----------------------+--------------------+
       |                    |                       |                       |                    |
       v                    v                       v                       v                    v
+--------------+    +---------------+       +---------------+       +---------------+    +---------------+
| 3D Structure |    | Genomic State |       | Needleman-    |       | In Silico     |    | PPI Network   |
| (PDB/Mol*)   |    | (ClinVar/SNP) |       | Wunsch Align  |       | Drug Screen   |    | (STRING DB)   |
+--------------+    +---------------+       +---------------+       +---------------+    +---------------+
| 4HHB / 2J6M  |    | HBB rs334     |       | BLOSUM62      |       | Lipinski Ro5  |    | Hub Centrality|
| Active Sites |    | Pathogenicity |       | Mutation Diff |       | Estimated ΔG  |    | Signalosomes  |
+--------------+    +---------------+       +---------------+       +---------------+    +---------------+
       |                    |                       |                       |                    |
       +--------------------+-----------------------+-----------------------+--------------------+
                                                    |
                                                    v
       +--------------------------------------------+--------------------------------------------+
       |                                                                                         |
       v                                                                                         v
+---------------------------------------+                       +---------------------------------------+
|  Interactive WebGL Dashboard v2.0     |                       |   Executive Publication Report v2.0   |
|  (3Dmol.js, SVG Interactome, Radars)  |                       |   (HTML / Printable PDF / JSON)       |
+---------------------------------------+                       +---------------------------------------+
```

---

## 🌟 Core Modules

### 1. 🔬 3D Structural Biology (`src/structural_biology.py`)
- Analyzes macromolecular complexes and PDB/mmCIF coordinates (`4HHB` - Human Deoxyhemoglobin).
- Identifies critical coordination sites:
  - **Proximal Histidine (`Chain A: His87`)**: Covalent coordination to heme iron ($Fe^{2+}$).
  - **Distal Histidine (`Chain A: His58`)**: Oxygen binding stabilization and steric resistance to CO.
  - **Sickle Cell Mutation Site (`Chain B/D: Glu6`)**: Glutamic acid to Valine substitution ($E6V$) causing pathological polymerization.

### 2. 🧬 Pairwise Sequence Alignment & Mutation Consequence (`src/sequence_aligner.py`)
- Needleman-Wunsch dynamic programming algorithm for global sequence alignment.
- Scored with the canonical **BLOSUM62** substitution matrix and affine gap penalties.
- Pinpoints precise amino acid substitutions, insertions, deletions, percent identity, and similarity:
  - `HBB_SICKLE`: $p.E7V$ substitution ($99.32\%$ identity, score $768$).
  - `EGFR_L858R`: $p.L858R$ kinase domain activation ($99.24\%$ identity, score $697$).
  - `KRAS_G12C`: $p.G12C$ codon 12 oncogenic switch ($99.41\%$ identity, score $842$).
  - `BRAF_V600E`: $p.V600E$ activation loop constitutive signaling ($99.42\%$ identity, score $885$).

### 3. 💊 In Silico Drug Screening & Lipinski Scoring (`src/docking_scorer.py`)
- Comprehensive **Lipinski Rule of 5** and **Veber Criteria** evaluation (MW $\le 500$, LogP $\le 5$, HBD $\le 5$, HBA $\le 10$, Rotatable Bonds $\le 10$, TPSA $\le 140\ \text{\AA}^2$).
- Computes Quantitative Estimate of Drug-likeness (QED proxy) from 0 to 100.
- Calculates an in silico binding affinity proxy ($\Delta G$ in $\text{kcal/mol}$) against target pockets:
  - **EGFR Kinase Cleft**: Rank #1 *Osimertinib* ($\Delta G = -12.52\ \text{kcal/mol}$, $97.3\%$ fitness).
  - **KRAS Switch-II Pocket**: Rank #1 *Sotorasib* ($\Delta G = -12.63\ \text{kcal/mol}$, $98.9\%$ fitness).
  - **BRAF V600E Pocket**: Rank #1 *Dabrafenib* ($\Delta G = -12.87\ \text{kcal/mol}$, $99.5\%$ fitness).
  - **Hemoglobin Tetramer Cleft**: Rank #1 *Voxelotor* ($\Delta G = -12.94\ \text{kcal/mol}$, $99.2\%$ fitness).

### 4. 🕸️ Protein-Protein Interaction (PPI) Interactome (`src/interaction_network.py`)
- Models high-confidence interactome topologies based on **STRING DB** standards.
- Computes degree centrality, edge weight connectivity, and network density.
- Identifies critical bottleneck hub proteins (e.g., $EGFR$, $KRAS$, $TP53$).

### 5. 🎯 Pan-Cancer Precision Oncology (`src/oncology_profiler.py`)
- Curated knowledgebase for high-impact oncology targets:
  - **EGFR**: Hotspots `p.L858R`, `p.T790M`, `Exon 19 del` $\rightarrow$ Targeted by *Osimertinib*, *Gefitinib*.
  - **KRAS**: Hotspots `p.G12C`, `p.G12D`, `p.G12V` $\rightarrow$ Targeted by *Sotorasib*, *Adagrasib*.
  - **BRAF**: Hotspots `p.V600E`, `p.V600K` $\rightarrow$ Targeted by *Dabrafenib + Trametinib*.
  - **TP53**: Hotspots `p.R175H`, `p.R248Q`, `p.R273H` $\rightarrow$ Structural & contact mutants.

### 6. ⚡ Dynamic Live Search Engine (`src/live_fetcher.py`)
- Real-time querying against live public APIs with an automatic local caching layer in `data/cache/`:
  - **RCSB PDB REST API**: Structural metadata and resolution.
  - **PubChem PUG-REST**: Real-time molecular weights, formulas, and lipophilicity.
  - **PubMed E-utilities**: Instant article search and abstract resolution.

### 7. 🧬 Genomic Variant & Clinical Significance (`src/variant_analyzer.py`)
- Integrates with **ClinVar** and **dbSNP** (`rs334`, `chr11:5227002:T>A`).
- Screens clinical interpretations:
  - `HBB:c.20A>T (p.Glu7Val / Glu6Val)` classified as **Pathogenic**.
  - Negative control comparison against `BRCA2:c.198A>G (rs28897700)` classified as **Benign**.

### 8. 📊 Biological Pathway Mapping (`src/pathway_mapper.py`)
- Maps human genes to curated **Reactome** pathways:
  - `R-HSA-1247673`: *Erythrocytes take up oxygen and release carbon dioxide*
  - `R-HSA-1480926`: *O2/CO2 exchange in erythrocytes*
  - `R-HSA-9927020`: *Heme assimilation*

### 9. 📚 Scientific Literature Mining (`src/literature_miner.py`)
- Searches **PubMed / PMC** for clinical trial breakthroughs:
  - Landmark Phase 3 trial on **Exagamglogene Autotemcel (Casgevy)** (*NEJM*, PMID: `38661449`).
  - Demonstrates CRISPR-Cas9 editing of the *BCL11A* erythroid enhancer to reactivate fetal hemoglobin (HbF).

### 10. 📑 Publication Report Generator (`src/report_generator.py`)
- Generates a standalone, responsive, interactive HTML scientific report (`reports/executive_summary.html`) with embedded visual charts, executive metrics, sequence alignment tracks, and multi-omics findings.

---

## 🚀 Quick Start

### 1. Installation
```bash
git clone https://github.com/goktugcici/computational-biology-pipeline.git
cd computational-biology-pipeline
pip install -r requirements.txt
```

### 2. Running the Pipeline

```bash
# 1. Run full pipeline with Pan-Cancer profiling and HTML report generation:
python src/main.py --pan_cancer EGFR --generate_report

# 2. Run Sequence Alignment for Sickle Cell or EGFR:
python src/main.py --align EGFR_L858R

# 3. Screen drug library against target binding pocket:
python src/main.py --screen_drugs KRAS

# 4. Analyze Protein-Protein Interaction (PPI) Interactome:
python src/main.py --ppi TP53

# 5. Export comprehensive analysis to JSON:
python src/main.py --export_json reports/pipeline_output.json

# 6. Live API querying (PDB, PubChem, PubMed):
python src/main.py --live_pdb 1TSR --live_drug imatinib --live_pubmed "immunotherapy NSCLC"
```

### 3. Launching the Interactive Web Dashboard

To inspect the 3D molecular viewer and multi-omics console locally:
```bash
python -m http.server 8000
```
Open **[http://localhost:8000/dashboard/index.html](http://localhost:8000/dashboard/index.html)** in your browser.

---

## 🧪 Automated Testing

Run the test suite covering all 11 multi-omics test cases:
```bash
python tests/test_pipeline.py
# or
pytest tests/
```
Output:
```
...........
----------------------------------------------------------------------
Ran 11 tests in 0.127s

OK
```

---

## 🐳 Docker Deployment

Run the containerized pipeline with a single command:
```bash
docker compose up
```

---

## 📜 License
This project is open-source and released under the **MIT License**.
