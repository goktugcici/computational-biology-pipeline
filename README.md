# Computational Biology & Multi-Omics Science Pipeline 🧬🔬

[![CI Pipeline](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml)
[![Pages Deployment](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/deploy-pages.yml/badge.svg)](https://goktugcici.github.io/computational-biology-pipeline/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version: 3.3 Precision Medicine](https://img.shields.io/badge/Release-v3.3_Precision_Medicine-purple.svg)](https://github.com/goktugcici/computational-biology-pipeline/releases)
[![Tests Passing](https://img.shields.io/badge/Tests-25%2F25%20Passed-brightgreen.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions)

A state-of-the-art, production-grade computational biology and precision oncology platform integrating **Pharmacogenomics (PGx) Clinical Toxicity Modeling (CPIC Level A)**, **Synthetic Lethality & CRISPR DepMap Dependency Mapping**, **Early Cancer Detection & 100,000+ Patient Statistical Risk Engine**, **Unified Relational Multi-Omics Knowledgebase (SQLite 3 - 10 Tables)**, **3D Structural Crystallography**, **Genomic Pathogenicity Annotation**, **Pan-Cancer Targeted Oncology**, **Needleman-Wunsch Sequence Alignment**, **In Silico Drug Screening (Lipinski Ro5 & Binding Free Energy)**, **Protein-Protein Interaction (PPI) Interactome Networks**, **Drug Combination Synergy Modeling (Bliss / Loewe CI)**, **Immunotherapy Biomarker Response Prediction (TMB & MSI-H)**, **AlphaFold Structural Quality (pLDDT & Disordered Regions)**, and an **Interactive WebGL Dashboard**.

---

## 🏛️ Comprehensive Architecture (v3.3 Precision Medicine Edition)

```
                            +-----------------------------------------------+
                            |   Computational Biology Platform v3.3         |
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
                    |                               |                               |
                    v                               v                               v
    +-------------------------------+   +-----------------------+   +-------------------------------+
    |  Pharmacogenomics (PGx) Safety|   |  Synthetic Lethality  |   |  100k Cancer Epidemiology DB  |
    |  (CPIC DPYD, UGT1A1, CYP2D6)  |   |  (DepMap Chronos Map) |   |  (NICE NG12 Red Flags & Bayes)|
    +---------------+---------------+   +-----------+-----------+   +---------------+---------------+
                    |                               |                               |
                    +-------------------------------+-------------------------------+
                                                    |
                                                    v
                               +----------------------------------------+
                               |  Unified SQLite 3 Bio-Knowledgebase    |
                               |  (10 Relational Multi-Omics Tables)    |
                               +--------------------+-------------------+
                                                    |
                                                    v
    +-----------------------------------------------+       +---------------------------------------+
    |  Interactive WebGL Dashboard & PGx Console    |       |   Executive Publication Report v3.3   |
    |  (3Dmol.js, PGx Evaluator, 100k Risk Gauge)   |       |   (HTML / Printable PDF / JSON)       |
    +-----------------------------------------------+       +---------------------------------------+
```

---

## 💊 Pharmacogenomics (PGx) Safety Engine (`src/pharmacogenomics.py`)

Implements **CPIC Level A/B** and **PharmGKB** clinical guidelines to calculate patient-specific drug toxicity risks, metabolic phenotypes, and dosing modifications:

1. **`DPYD` (Fluoropyrimidines: 5-FU, Capecitabine)**:
   - Detects loss-of-function alleles: `*2A` ($c.1905+1G>A$), `*13` ($p.I560S$), $c.2846A>T$ ($p.D949V$), and `HapB3`.
   - Computes DPYD Activity Score ($0.0$ to $2.0$).
   - Flags **Intermediate Metabolizers** for a mandatory **50% initial dose reduction** with therapeutic drug monitoring (TDM).
   - Flags **Poor Metabolizers (Activity Score $\le 0.5$)** as **CONTRAINDICATED** due to lethal mucositis, pancytopenia, and septic shock risk.
2. **`UGT1A1` (Irinotecan / FOLFIRI / FOLFIRINOX)**:
   - Evaluates $(TA)_6$ vs $(TA)_7$ promoter repeats (`*28` allele, Gilbert syndrome variant).
   - Recommends at least **25–30% Irinotecan dose reduction** in homozygous `*28/*28` poor metabolizers to prevent life-threatening Grade 4 neutropenia and severe diarrhea.
3. **`CYP2D6` (Tamoxifen in ER+ Breast Cancer)**:
   - Evaluates bioactivation of prodrug Tamoxifen to potent active anti-estrogen **Endoxifen**.
   - Identifies poor metabolizers (`*4/*4`, `*4/*5`, `*5/*5`) with subtherapeutic endoxifen levels and recommends switching postmenopausal patients to **Aromatase Inhibitors (Anastrozole/Letrozole)**.

---

## 🎯 Synthetic Lethality & DepMap Dependency Engine (`src/synthetic_lethality.py`)

Models CRISPR knockout essentiality scores from the **Broad Institute Cancer Dependency Map (DepMap Chronos / Project Score)** to uncover exploitable genetic vulnerabilities:

- **`BRCA1/2 Deficiency` $\iff$ `PARP1` Dependency**: Homologous recombination deficiency (HRD+) causes hypersensitivity to PARP trapping (Chronos: $-1.24$, Vulnerability: $98.0\%$). Targeted by *Olaparib, Talazoparib*.
- **`MSI-H / dMMR` $\iff$ `WRN` RecQ Helicase Dependency**: Microsatellite-unstable tumors require WRN to resolve non-B DNA secondary hairpins (Chronos: $-1.45$, Vulnerability: $95.0\%$). Targeted by *HRO761, RO7589831*.
- **`MTAP Deletion (9p21)` $\iff$ `PRMT5 / MAT2A` Dependency**: MTA accumulation creates synthetic lethal vulnerability to PRMT5/MAT2A inhibition (Chronos: $-0.92$, Vulnerability: $88.0\%$). Targeted by *MRTX1719, AG-270*.
- **`ARID1A Mutation` $\iff$ `EZH2` Dependency**: SWI/SNF loss confers dependency on PRC2/EZH2 histone methyltransferase (Chronos: $-0.85$, Vulnerability: $84.0\%$). Targeted by *Tazemetostat*.
- **`KRAS G12C/D/V` $\iff$ `SHP2 / SOS1` Feedback Dependency**: Prevents RTK feedback reactivation during direct KRAS inhibition (Chronos: $-0.78$, Vulnerability: $82.0\%$). Targeted by *TNO155 + Sotorasib*.

---

## 🩺 Early Cancer Detection & 100,000+ Patient Risk Engine (`src/cancer_risk_predictor.py`)

- **100,000+ Patient Population Cohort (`data/cancer_epidemiology_100k.db`)**: Calibrated against **CDC WONDER**, **SEER Cancer Statistics Review**, and **WHO IARC**.
- **14 Clinical Alarm Symptoms ("Red Flags" - NICE NG12 & QCancer)**: Hemoptysis, persistent cough, unexplained weight loss, rectal bleeding, change in bowel habits, iron deficiency anemia, palpable breast lump, nipple abnormality, gross hematuria, urinary hesitancy, painless jaundice, progressive dysphagia, persistent dyspepsia, and suspicious mole changes (ABCDE).
- **NICE NG12 Urgent 2-Week-Wait (2WW) Referral Protocol**: Triggers immediate clinical pathways when symptoms exceed the clinical $\ge 3\%$ Positive Predictive Value (PPV) threshold.

---

## 🗄️ Unified Multi-Omics Relational Database (`src/database.py`)

An embedded, zero-configuration relational knowledgebase (`SQLite 3`) in `data/bio_knowledgebase.db` with **10 normalized tables**:
1. `genes` (10 target oncogenes & tumor suppressors)
2. `variants` (28 ClinVar & dbSNP pathogenic variants)
3. `structures` (7 crystallographic PDB coordinates & AlphaFold models)
4. `compounds` (11 targeted therapeutic small molecules)
5. `oncology_hotspots` (16 actionable driver mutations)
6. `pathways` (7 Reactome cellular signaling cascades)
7. `drug_synergy` (8 validated dual-agent combination regimens)
8. `literature_trials` (10 landmark NEJM Phase 2/3 clinical studies)
9. `pharmacogenomics_guidelines` (6 CPIC Level A dosing guidelines)
10. `synthetic_lethal_pairs` (6 DepMap validated genetic dependencies)

---

## 🌟 14 Integrated Multi-Omics Modules

1. **💊 Pharmacogenomics (PGx) Safety Engine (`src/pharmacogenomics.py`)**: CPIC Level A dosing rules for DPYD, UGT1A1, and CYP2D6.
2. **🎯 Synthetic Lethality & DepMap Engine (`src/synthetic_lethality.py`)**: CRISPR Chronos essentiality and synthetic lethal drug matching.
3. **🩺 Early Cancer Detection & 100k Risk Engine (`src/cancer_risk_predictor.py`)**: Population cohort matching and NICE 2WW criteria.
4. **🗄️ Relational Bio-Knowledgebase (`src/database.py`)**: SQLite 3 engine with 10 multi-omics tables and JSON exports.
5. **🔬 3D Structural Biology (`src/structural_biology.py`)**: Macromolecular coordinate analysis, catalytic histidines, and heme pockets (`4HHB`).
6. **🧬 Needleman-Wunsch Sequence Alignment (`src/sequence_aligner.py`)**: Dynamic programming with BLOSUM62 matrix for global amino acid alignments.
7. **💊 In Silico Drug Screening & Lipinski Scoring (`src/docking_scorer.py`)**: Lipinski Rule of 5, Veber criteria, QED proxy, and binding free energy estimation ($\Delta G$).
8. **💥 Combination Therapy & Drug Synergy Predictor (`src/synergy_predictor.py`)**: Loewe Combination Index (CI) and Bliss Independence models.
9. **🛡️ Immunotherapy Biomarker & TMB Response Predictor (`src/immunotherapy_response.py`)**: TMB ($\ge 10\ \text{mut/Mb}$), MSI-H/dMMR, and PD-L1 TPS evaluation.
10. **📐 AlphaFold Structural Quality & Disorder (`src/alphafold_quality.py`)**: Per-residue pLDDT confidence bands and Intrinsically Disordered Regions (IDRs).
11. **🕸️ Protein-Protein Interaction (PPI) Interactome (`src/interaction_network.py`)**: STRING-calibrated network topology, centrality hubs, and bottleneck proteins.
12. **🎯 Pan-Cancer Precision Oncology (`src/oncology_profiler.py`)**: Clinically actionable hotspot catalogs for *EGFR, KRAS, BRAF, TP53*.
13. **⚡ Dynamic Live Scientific API Search (`src/live_fetcher.py`)**: Real-time REST integrations for RCSB PDB, PubChem, and PubMed with smart caching.
14. **📑 Publication Report Generator (`src/report_generator.py`)**: Standalone publication-ready HTML reports with embedded charts and metrics.

---

## 🚀 Quick Start

### 1. Installation
```bash
git clone https://github.com/goktugcici/computational-biology-pipeline.git
cd computational-biology-pipeline
pip install -r requirements.txt
```

### 2. Precision Medicine & Pharmacogenomics CLI Usage

```bash
# 1. Screen DPYD diplotype for 5-FU / Capecitabine toxicity risk:
python src/main.py --dpyd_alleles *1 *2A

# 2. Screen UGT1A1 and CYP2D6 diplotypes:
python src/main.py --ugt1a1 "*28/*28" --cyp2d6 "*4/*4"

# 3. Query synthetic lethal dependency for BRCA1 or MSI-H:
python src/main.py --synthetic_lethal BRCA1

# 4. Evaluate early cancer risk with 100k cohort matching:
python src/main.py --predict_cancer --patient_age 62 --patient_gender M --patient_symptoms "hemoptysis,chronic_cough,weight_loss" --smoking_pack_years 30.0 --family_history

# 5. Full Pipeline Run with Pan-Cancer, Drug Synergy, TMB, and HTML report:
python src/main.py --pan_cancer EGFR --synergy Dabrafenib Trametinib --tmb 16.2 --msi MSI-H --db_query EGFR --generate_report
```

### 3. Interactive Web Dashboard
```bash
python -m http.server 8000
```
Open **[http://localhost:8000/dashboard/index.html](http://localhost:8000/dashboard/index.html)** in your browser:
- **💊 PGx & Sentetik Letalite**: Interactive DPYD, UGT1A1, and CYP2D6 diplotype toxicity evaluator, CPIC dosing calculator, and DepMap Chronos synthetic lethality visualizer.
- **🩺 Erken Teşhis (100k)**: Real-time risk gauge, 14 interactive NICE alarm symptom pills, 100,000-patient cohort matching card, and organ-specific progress bars.
- **🗄️ Veri Tabanı (Database Explorer)**: Instant search, category filtering, and JSON export over the 10 SQLite relational tables.
- **🧬 3B Molekül Görselleştirici**: Interactive 3Dmol.js WebGL crystallographic rendering.
- **🛡️ İmmünoterapi Simülatörü**: Interactive TMB, MSI, and PD-L1 response prediction sliders.

---

## 🧪 Automated Testing

Covering all 25 unit tests across pharmacogenomics CPIC rules, synthetic lethality DepMap scores, epidemiology cohort generation, NICE guidelines, multi-organ Bayesian risk models, SQLite queries, structural biology, and sequence alignments:
```bash
python -m unittest tests/test_pipeline.py
# or
pytest tests/
```
Output:
```
.........................
----------------------------------------------------------------------
Ran 25 tests in 0.787s

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
