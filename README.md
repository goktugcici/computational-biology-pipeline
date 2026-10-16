# Computational Biology & Multi-Omics Science Pipeline 🧬🔬

[![CI Pipeline](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml)
[![Pages Deployment](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/deploy-pages.yml/badge.svg)](https://goktugcici.github.io/computational-biology-pipeline/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version: 3.2 Early Detection](https://img.shields.io/badge/Release-v3.2_Early_Detection-emerald.svg)](https://github.com/goktugcici/computational-biology-pipeline/releases)
[![Tests Passing](https://img.shields.io/badge/Tests-20%2F20%20Passed-brightgreen.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions)

A state-of-the-art, production-grade computational biology and precision oncology platform integrating **Early Cancer Detection & 100,000+ Patient Statistical Risk Engine**, **Unified Relational Multi-Omics Knowledgebase (SQLite 3)**, **3D Structural Crystallography**, **Genomic Pathogenicity Annotation**, **Pan-Cancer Targeted Oncology**, **Needleman-Wunsch Sequence Alignment**, **In Silico Drug Screening (Lipinski Ro5 & Binding Free Energy)**, **Protein-Protein Interaction (PPI) Interactome Networks**, **Drug Combination Synergy Modeling (Bliss / Loewe CI)**, **Immunotherapy Biomarker Response Prediction (TMB & MSI-H)**, **AlphaFold Structural Quality (pLDDT & Disordered Regions)**, and an **Interactive WebGL Dashboard with Live 100k Cancer Risk Simulator**.

---

## 🏛️ Comprehensive Architecture (v3.2 Early Detection Edition)

```
                            +-----------------------------------------------+
                            |   Computational Biology Platform v3.2         |
                            +-----------------------+-----------------------+
                                                    |
        +--------------------+----------------------+-----------------------+--------------------+
        |                    |                      |                       |                    |
        v                    v                      v                       v                    v
+--------------+    +---------------+      +-----------------+      +---------------+    +---------------+
| 3D Structure |    | Genomic State |      | Needleman-Wunsch|      | In Silico     |    | Drug Synergy  |
| & AlphaFold  |    | (ClinVar/SNP) |      | Alignment       |      | Drug Screen   |    | (Loewe / Bliss|
+--------------+    +---------------+      +-----------------+      +---------------+    +---------------+
| PDB & pLDDT  |    | HBB rs334     |      | BLOSUM62        |      | Lipinski Ro5  |    | Combination CI|
| Disordered   |    | Pathogenicity |      | Mutation Diff   |      | Estimated ΔG  |    | Bypass Break  |
+--------------+    +---------------+      +-----------------+      +---------------+    +---------------+
        |                    |                      |                       |                    |
        +--------------------+----------------------+-----------------------+--------------------+
                                                    |
                                                    v
                    +-------------------------------+-------------------------------+
                    |                                                               |
                    v                                                               v
    +---------------------------------------+       +-----------------------------------------------+
    |  Unified SQLite 3 Bio-Knowledgebase   |       |  100,000+ Patient Cancer Epidemiology DB     |
    |  (8 Relational Multi-Omics Tables)    |       |  (NICE NG12 Red Flags & Bayesian Risk Model)  |
    +-------------------+-------------------+       +-----------------------+-----------------------+
                        |                                                   |
                        +---------------------------+-----------------------+
                                                    |
                                                    v
    +-----------------------------------------------+       +---------------------------------------+
    |  Interactive WebGL Dashboard & Risk Simulator |       |   Executive Publication Report v3.2   |
    |  (3Dmol.js, 100k Risk Gauge, DB Explorer, etc)|       |   (HTML / Printable PDF / JSON)       |
    +-----------------------------------------------+       +---------------------------------------+
```

---

## 🩺 Early Cancer Detection & 100,000+ Patient Risk Engine (`src/cancer_risk_predictor.py`)

The pipeline implements an advanced, clinically grounded statistical estimation engine for early cancer detection:

1. **100,000+ Patient Population Cohort (`data/cancer_epidemiology_100k.db`)**:
   - A fully indexed SQLite database holding 100,000 population-calibrated patient records modeled on **CDC WONDER**, **SEER Cancer Statistics Review**, and **NICE NG12** primary care data.
   - High-speed bitwise masking indexes allow sub-millisecond query matching for combinations of symptoms, age brackets, smoking pack-years, and family history.

2. **14 Clinical Alarm Symptoms ("Red Flags" - NICE NG12 & QCancer Guidelines)**:
   - **Hemoptysis (`hemoptysis`)**: Coughing blood (OR $8.5$, PPV $5.5\%$, Lung Cancer primary driver).
   - **Persistent Cough (`chronic_cough`)**: Cough $>3$ weeks / hoarseness (OR $3.2$, PPV $2.4\%$, Lung/Larynx).
   - **Unexplained Weight Loss (`unexplained_weight_loss`)**: $>5-10\%$ loss in 6 months (OR $4.8$, PPV $3.8\%$, Pan-cancer red flag).
   - **Rectal Bleeding (`rectal_bleeding`)**: Blood in stool / hematochezia (OR $7.5$, PPV $4.8\%$, Colorectal).
   - **Bowel Habit Change (`bowel_habit_change`)**: Persistent alteration $>4-6$ weeks (OR $4.2$, PPV $3.2\%$, Colorectal).
   - **Iron Deficiency Anemia (`iron_deficiency_anemia`)**: Microcytic anemia (OR $5.2$, PPV $6.2\%$, Colorectal & Gastric).
   - **Palpable Breast Lump (`breast_lump`)**: Firm, painless mass (OR $12.0$, PPV $11.5\%$, Breast Cancer).
   - **Nipple Abnormality (`nipple_abnormality`)**: Spontaneous discharge / inversion (OR $6.0$, PPV $4.5\%$, Breast).
   - **Gross Hematuria (`hematuria`)**: Blood in urine (OR $7.2$, PPV $5.8\%$, Bladder / Renal / Prostate).
   - **Urinary Hesitancy (`urinary_hesitancy`)**: Weak stream / nocturia (OR $2.8$, PPV $2.1\%$, Prostate).
   - **Painless Jaundice (`jaundice`)**: Skin/scleral yellowing, dark urine (OR $18.5$, PPV $22.0\%$, Pancreatic Head).
   - **Progressive Dysphagia (`dysphagia`)**: Difficulty swallowing solids/liquids (OR $8.4$, PPV $7.5\%$, Gastric / Esophageal).
   - **Persistent Dyspepsia (`persistent_dyspepsia`)**: Epigastric pain, early satiety (OR $3.1$, PPV $2.5\%$, Gastric).
   - **Suspicious Mole Changes (`skin_mole_changes`)**: ABCDE criteria alteration (OR $14.2$, PPV $12.5\%$, Cutaneous Melanoma).

3. **NICE NG12 Urgent Suspected Cancer Pathway (2-Week-Wait / 2WW)**:
   - Automatically evaluates if symptoms meet or exceed the **$\ge 3\%$ Positive Predictive Value (PPV)** threshold established by the UK National Institute for Health and Care Excellence (NICE).
   - Triggers actionable referral guidance and prioritizes specialized diagnostic imaging (Low-dose CT, Colonoscopy, Mammography, Upper GI Endoscopy, etc.).

---

## 🗄️ Unified Multi-Omics Relational Database (`src/database.py`)

The pipeline features an embedded, zero-configuration relational knowledgebase (`SQLite 3`) stored in `data/bio_knowledgebase.db` with 8 normalized tables:
- **`genes`**: Target oncogenes & tumor suppressors (*EGFR, KRAS, BRAF, TP53, HBB, BRCA1, BRCA2, MET, ALK, PIK3CA*).
- **`variants`**: ClinVar & dbSNP clinical records, HGVS cDNA/protein notation, chromosomal loci, allele frequencies, and clinical pathogenicity.
- **`structures`**: Macromolecular crystallography coordinates (PDB IDs) & AlphaFold predictive models with mean pLDDT confidence scores.
- **`compounds`**: Targeted therapies & small molecules (*Osimertinib, Sotorasib, Dabrafenib, Trametinib, Voxelotor, Metformin*), Lipinski Ro5 compliance, and QED scores.
- **`oncology_hotspots`**: Actionable cancer mutations across NSCLC, CRC, Melanoma, and PDAC, linked to targeted regimens.
- **`pathways`**: Reactome cellular events and signal cascades (RTK, MAPK, p53 checkpoints, oxygen transport).
- **`drug_synergy`**: Clinically validated dual-agent regimens, Loewe Combination Index ($\text{CI}$), and Bliss excess scores.
- **`literature_trials`**: Landmark Phase 2/3 oncology and hematology studies published in *The New England Journal of Medicine (NEJM)*.

---

## 🌟 12 Integrated Multi-Omics Modules

1. **🩺 Early Cancer Detection & 100k Risk Engine (`src/cancer_risk_predictor.py`)**: Population cohort matching, NICE 2WW criteria, and multi-organ Bayesian log-odds prediction.
2. **🗄️ Relational Bio-Knowledgebase (`src/database.py`)**: SQLite 3 engine for relational data queries and JSON exports.
3. **🔬 3D Structural Biology (`src/structural_biology.py`)**: Macromolecular coordinate analysis, catalytic histidines, and heme pockets (`4HHB`).
4. **🧬 Needleman-Wunsch Sequence Alignment (`src/sequence_aligner.py`)**: Dynamic programming with BLOSUM62 matrix for global amino acid alignments.
5. **💊 In Silico Drug Screening & Lipinski Scoring (`src/docking_scorer.py`)**: Lipinski Rule of 5, Veber criteria, QED proxy, and binding free energy estimation ($\Delta G$).
6. **💥 Combination Therapy & Drug Synergy Predictor (`src/synergy_predictor.py`)**: Loewe Combination Index (CI) and Bliss Independence models.
7. **🛡️ Immunotherapy Biomarker & TMB Response Predictor (`src/immunotherapy_response.py`)**: TMB ($\ge 10\ \text{mut/Mb}$), MSI-H/dMMR, and PD-L1 TPS evaluation.
8. **📐 AlphaFold Structural Quality & Disorder (`src/alphafold_quality.py`)**: Per-residue pLDDT confidence bands and Intrinsically Disordered Regions (IDRs).
9. **🕸️ Protein-Protein Interaction (PPI) Interactome (`src/interaction_network.py`)**: STRING-calibrated network topology, centrality hubs, and bottleneck proteins.
10. **🎯 Pan-Cancer Precision Oncology (`src/oncology_profiler.py`)**: Clinically actionable hotspot catalogs for *EGFR, KRAS, BRAF, TP53*.
11. **⚡ Dynamic Live Scientific API Search (`src/live_fetcher.py`)**: Real-time REST integrations for RCSB PDB, PubChem, and PubMed with smart caching.
12. **📑 Publication Report Generator (`src/report_generator.py`)**: Standalone publication-ready HTML reports with embedded charts and metrics.

---

## 🚀 Quick Start

### 1. Installation
```bash
git clone https://github.com/goktugcici/computational-biology-pipeline.git
cd computational-biology-pipeline
pip install -r requirements.txt
```

### 2. Early Cancer Detection CLI Usage

```bash
# 1. Evaluate early cancer risk for a 62 y/o male smoker with hemoptysis and weight loss:
python src/main.py --predict_cancer --patient_age 62 --patient_gender M --patient_symptoms "hemoptysis,chronic_cough,weight_loss" --smoking_pack_years 30.0 --family_history

# 2. Evaluate breast cancer risk for a 52 y/o female with palpable breast lump:
python src/main.py --predict_cancer --patient_age 52 --patient_gender F --patient_symptoms "breast_lump,nipple_abnormality"

# 3. Evaluate colorectal cancer risk with rectal bleeding and iron-deficiency anemia:
python src/main.py --predict_cancer --patient_age 58 --patient_gender M --patient_symptoms "rectal_bleeding,iron_deficiency_anemia"

# 4. Generate or re-seed the 100,000-patient epidemiological cohort:
python src/main.py --build_cohort 100000
```

### 3. Comprehensive Pipeline Run
```bash
python src/main.py --pan_cancer EGFR --synergy Dabrafenib Trametinib --tmb 16.2 --msi MSI-H --db_query EGFR --generate_report
```

### 4. Interactive Web Dashboard & 100k Risk Simulator
```bash
python -m http.server 8000
```
Open **[http://localhost:8000/dashboard/index.html](http://localhost:8000/dashboard/index.html)** in your browser:
- **🩺 Erken Teşhis (100k)**: Real-time risk gauge, 14 interactive NICE alarm symptom pills, 100,000-patient cohort matching card, organ-specific progress bars, and evidence-based clinical investigation recommendations.
- **🗄️ Veri Tabanı (Database Explorer)**: Search, category filtering, and JSON export over the 8 SQLite relational tables.
- **🧬 3B Molekül Görselleştirici**: Interactive 3Dmol.js WebGL crystallographic rendering.
- **🛡️ İmmünoterapi Simülatörü**: Interactive TMB, MSI, and PD-L1 response prediction sliders.

---

## 🧪 Automated Testing

Covering all 20 unit tests across epidemiology cohort generation, NICE guidelines, multi-organ Bayesian risk models, SQLite queries, structural biology, and sequence alignments:
```bash
python -m unittest tests/test_pipeline.py
# or
pytest tests/
```
Output:
```
....................
----------------------------------------------------------------------
Ran 20 tests in 0.715s

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
