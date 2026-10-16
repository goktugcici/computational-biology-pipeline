# Computational Biology & Multi-Omics Science Pipeline 🧬🔬

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Scientific Suite: Antigravity Science](https://img.shields.io/badge/Suite-Antigravity%20Science-emerald.svg)](#)

A modular, reproducible bioinformatics pipeline that integrates **3D Structural Biology**, **Genomic Variant Annotation**, **Pharmacological Drug Discovery**, **Biological Pathway Mapping**, and **Biomedical Literature Intelligence**.

---

## 🏛️ Architecture Overview

```
                      +-----------------------------------+
                      |   Computational Biology Pipeline  |
                      +-----------------+-----------------+
                                        |
       +-------------------+------------+------------+--------------------+
       |                   |                         |                    |
       v                   v                         v                    v
+--------------+   +---------------+         +---------------+    +---------------+
| 3D Structure |   | Genomic State |         | Pharmacology  |    |  Literature   |
| (PDB/Mol*)   |   | (ClinVar/SNP) |         | (PubChem/ChEM)|    | (PubMed/PMC)  |
+--------------+   +---------------+         +---------------+    +---------------+
       |                   |                         |                    |
       | 4HHB Structure    | HBB rs334               | Voxelotor (Oxbryta)| CRISPR Exa-cel
       | Active Histidines | Pathogenic Missense     | Metformin (ChEMBL) | NEJM Landmark
       v                   v                         v                    v
       +-------------------+-------------------------+--------------------+
                                        |
                                        v
                      +-----------------------------------+
                      |      Multi-Omics Unified Model    |
                      +-----------------------------------+
```

---

## 🌟 Core Modules

### 1. 3D Structural Biology (`src/structural_biology.py`)
- Analyzes macromolecular complexes and PDB/mmCIF coordinates (`4HHB` - Human Deoxyhemoglobin).
- Identifies critical coordination sites:
  - **Proximal Histidine (`Chain A: His87`)**: Covalent coordination to heme iron ($Fe^{2+}$).
  - **Distal Histidine (`Chain A: His58`)**: Oxygen binding stabilization and steric resistance to CO.
  - **Sickle Cell Mutation Site (`Chain B/D: Glu6`)**: Glutamic acid to Valine substitution ($E6V$) causing pathological polymerization.

### 2. Genomic Variant & Clinical Significance (`src/variant_analyzer.py`)
- Integrates with **ClinVar** and **dbSNP** (`rs334`, `chr11:5227002:T>A`).
- Screens clinical interpretations:
  - `HBB:c.20A>T (p.Glu7Val / Glu6Val)` classified as **Pathogenic**.
  - Negative control comparison against `BRCA2:c.198A>G (rs28897700)` classified as **Benign** by expert panels.

### 3. Pharmacology & Cheminformatics (`src/pharmacology.py`)
- Retrieves chemical properties, SMILES structures, and mechanisms of action from **PubChem** and **ChEMBL**:
  - **Voxelotor (Oxbryta)** (PubChem CID: `71602803`): Binds reversibly to the N-terminal valine of hemoglobin's $\alpha$-chain to allosterically increase oxygen affinity, stabilizing non-polymerizing oxy-HbS.
  - **Metformin** (`CHEMBL1431`): Standard profiling for metabolic regulators.

### 4. Biological Pathway Mapping (`src/pathway_mapper.py`)
- Maps human genes to curated **Reactome** pathways:
  - `R-HSA-1247673`: *Erythrocytes take up oxygen and release carbon dioxide*
  - `R-HSA-1480926`: *O2/CO2 exchange in erythrocytes*
  - `R-HSA-9927020`: *Heme assimilation*

### 5. Scientific Literature Mining (`src/literature_miner.py`)
- Searches **PubMed / PMC** for clinical trial breakthroughs:
  - Landmark Phase 3 trial on **Exagamglogene Autotemcel (Casgevy)** (*NEJM*, PMID: `38661449`).
  - Demonstrates CRISPR-Cas9 editing of the *BCL11A* erythroid enhancer to reactivate fetal hemoglobin (HbF).

---

## 🚀 Quick Start

### 1. Installation
Clone the repository and install dependencies:
```bash
git clone https://github.com/<your-username>/computational-biology-pipeline.git
cd computational-biology-pipeline
pip install -r requirements.txt
```

### 2. Run the Full Pipeline
Execute the end-to-end multi-omics analysis:
```bash
python src/main.py --data_dir data
```

---

## 📂 Project Structure

```
computational-biology-pipeline/
├── src/
│   ├── main.py                 # Pipeline orchestrator
│   ├── structural_biology.py   # PDB/CIF coordinate parser
│   ├── variant_analyzer.py     # ClinVar/dbSNP variant classifier
│   ├── pharmacology.py         # PubChem/ChEMBL chemical profiler
│   ├── pathway_mapper.py       # Reactome biological pathway mapper
│   └── literature_miner.py     # PubMed literature intelligence
├── data/                       # Curated scientific data files
│   ├── pdb_00004hhb.cif        # Human Deoxyhemoglobin coordinates
│   ├── summary_all.json        # ClinVar variant classifications
│   ├── dbsnp_rs334.json        # dbSNP genomic coordinates
│   ├── voxelotor_props.json    # PubChem properties for Voxelotor
│   ├── voxelotor_pharm.json    # Pharmacology & mechanism of action
│   ├── chembl_metformin.json   # ChEMBL drug records
│   ├── reactome_hbb.json       # Reactome pathways for HBB
│   └── pubmed_abstract.json    # NEJM CRISPR paper metadata
├── reports/
│   └── scientific_investigation_report.md # Comprehensive scientific report
├── requirements.txt            # Python dependencies
├── .gitignore                  # Git ignore rules
└── README.md                   # Documentation
```

---

## 📄 License
This project is open-source under the MIT License.
