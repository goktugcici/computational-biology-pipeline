# Computational Biology & Multi-Omics Science Pipeline 🧬🔬

[![CI Pipeline](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/ci.yml)
[![Pages Deployment](https://github.com/goktugcici/computational-biology-pipeline/actions/workflows/deploy-pages.yml/badge.svg)](https://goktugcici.github.io/computational-biology-pipeline/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version: 3.4 Clinical Oncology](https://img.shields.io/badge/Release-v3.4_Clinical_Oncology-purple.svg)](https://github.com/goktugcici/computational-biology-pipeline/releases)
[![Tests Passing](https://img.shields.io/badge/Tests-30%2F30%20Passed-brightgreen.svg)](https://github.com/goktugcici/computational-biology-pipeline/actions)

A state-of-the-art, production-grade computational biology, translational research, and precision oncology platform integrating **Kaplan-Meier Clinical Survival Analysis & Log-Rank Modeler**, **Liquid Biopsy & ctDNA Minimal Residual Disease (MRD) Tracking**, **Pharmacogenomics (PGx) Clinical Toxicity Modeling (CPIC Level A)**, **Synthetic Lethality & CRISPR DepMap Dependency Mapping**, **Early Cancer Detection & 100,000+ Patient Statistical Risk Engine**, **Unified Relational Multi-Omics Knowledgebase (SQLite 3 - 12 Tables)**, **3D Structural Crystallography**, **Genomic Pathogenicity Annotation**, **Pan-Cancer Targeted Oncology**, **Needleman-Wunsch Sequence Alignment**, **In Silico Drug Screening (Lipinski Ro5 & Binding Free Energy)**, **Protein-Protein Interaction (PPI) Interactome Networks**, **Drug Combination Synergy Modeling (Bliss / Loewe CI)**, **Immunotherapy Biomarker Response Prediction (TMB & MSI-H)**, **AlphaFold Structural Quality (pLDDT & Disordered Regions)**, and an **Interactive WebGL Dashboard**.

---

## 🏛️ Comprehensive Architecture (v3.4 Clinical Oncology Edition)

```
                            +-----------------------------------------------+
                            |   Computational Biology Platform v3.4         |
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
                    +-------------------------------+-------------------------------+
                    |                                                               |
                    v                                                               v
    +-------------------------------+                               +-------------------------------+
    |  Kaplan-Meier Survival Model  |                               |  Liquid Biopsy ctDNA MRD      |
    |  (Log-Rank, Greenwood SE, HR) |                               |  (VAF Dynamics & Resistance)  |
    +---------------+---------------+                               +---------------+---------------+
                    |                                                               |
                    +-------------------------------+-------------------------------+
                                                    |
                                                    v
                               +----------------------------------------+
                               |  Unified SQLite 3 Bio-Knowledgebase    |
                               |  (12 Relational Multi-Omics Tables)   |
                               +--------------------+-------------------+
                                                    |
                                                    v
    +-----------------------------------------------+       +---------------------------------------+
    |  Interactive WebGL Dashboard & Oncology Suite |       |   Executive Publication Report v3.4   |
    |  (3Dmol.js, KM Stepped Curve, ctDNA Monitor)  |       |   (HTML / Printable PDF / JSON)       |
    +-----------------------------------------------+       +---------------------------------------+
```

---

## 📈 Clinical Survival Analysis & Log-Rank Modeler (`src/survival_analysis.py`)

Implements clinical biostatistics and time-to-event survival modeling with zero external dependencies using pure standard mathematical algorithms:

- **Kaplan-Meier Non-Parametric Estimator**:
  $$S(t) = \prod_{t_i \le t} \left(1 - \frac{d_i}{n_i}\right)$$
  Calculates cumulative survival probabilities over time, handles censored observations, and reconstructs monotonic step functions.
- **Greenwood's Standard Error Formula**:
  $$\widehat{\text{Var}}(S(t)) = S(t)^2 \sum_{t_i \le t} \frac{d_i}{n_i(n_i - d_i)}$$
  Yields rigorous 95% confidence intervals $[S(t) \pm 1.96 \cdot \text{SE}(S(t))]$.
- **Mantel-Cox Log-Rank Hypothesis Test**:
  $$\chi^2 = \frac{(O_1 - E_1)^2}{V}$$
  Computes test statistic, asymptotic chi-square p-value, and verifies statistical superiority between experimental targeted regimens and standard control therapies.
- **Peto Hazard Ratio (HR) & 95% Confidence Intervals**:
  $$\text{HR} = \exp\left(\frac{O_1 - E_1}{V}\right), \quad 95\%\ \text{CI} = \exp\left(\ln(\text{HR}) \pm 1.96 \cdot \frac{1}{\sqrt{V}}\right)$$
- **Curated Phase 3 Benchmark Trial Cohorts**:
  1. **FLAURA (NEJM)**: Osimertinib vs Gefitinib in untreated EGFR-mutated NSCLC ($\text{HR} = 0.80$, $p = 0.046$).
  2. **KEYNOTE-042 (Lancet)**: Pembrolizumab monotherapy vs Chemotherapy in TMB-High / PD-L1 TPS $\ge 50\%$ NSCLC ($\text{HR} = 0.69$, $p = 0.0003$).
  3. **OlympiAD (NEJM)**: Olaparib vs Physician's choice chemotherapy in metastatic gBRCA-mutated breast cancer ($\text{HR} = 0.58$, $p = 0.0009$).
  4. **CodeBreaK 200 (Lancet)**: Sotorasib vs Docetaxel in pretreated KRAS G12C NSCLC ($\text{HR} = 0.66$, $p = 0.002$).

---

## 🩸 Liquid Biopsy & ctDNA Minimal Residual Disease (MRD) Tracking (`src/liquid_biopsy.py`)

Simulates and evaluates serial cell-free circulating tumor DNA (ctDNA) plasma monitoring for non-invasive longitudinal oncology:

- **Longitudinal Variant Allele Frequency (VAF %)**: Detects subclonal fluctuations with analytical sensitivity down to $0.01\%$.
- **Molecular Response Classification**:
  - **`mCR` (Molecular Complete Response)**: ctDNA clearance below detection limit ($<0.05\%$), indicating MRD negativity and favorable progression-free survival.
  - **`mPR` (Molecular Partial Response)**: Relative VAF reduction $\le -50\%$.
  - **`mSD` (Molecular Stable Disease)**: Persistent clonality within $-50\%$ to $+25\%$.
  - **`mPD` (Molecular Progressive Disease)**: Clonal expansion $>+25\%$, providing **3 to 6 months early warning lead time before radiologic RECIST recurrence**.
- **Acquired Resistance Clones & Actionable Bypass Therapies**:
  - **`EGFR C797S`** (Loss of covalent Osimertinib binding) $\implies$ Allosteric 4th-gen EGFR TKIs (*BLU-945, BDTX-1535*) or HER3-DXd (*Patritumab Deruxtecan*).
  - **`MET Amplification`** (ERBB3 / PI3K bypass activation) $\implies$ Combination TKI + MET inhibitor (*Osimertinib + Savolitinib*).
  - **`KRAS Y96D`** (Switch-II pocket steric resistance) $\implies$ Pan-KRAS inhibitors (*RMC-6236*) or SOS1/SHP2 inhibitors.
  - **`BRCA Reversion Mutations`** (Reading frame restoration) $\implies$ ATR inhibitors (*Ceralasertib*), WEE1 inhibitors (*Adavosertib*), or platinum re-challenge.

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

An embedded, zero-configuration relational knowledgebase (`SQLite 3`) in `data/bio_knowledgebase.db` with **12 normalized tables**:
1. `genes` (10 target oncogenes & tumor suppressors)
2. `variants` (42 ClinVar & dbSNP pathogenic variants)
3. `structures` (7 crystallographic PDB coordinates & AlphaFold models)
4. `compounds` (11 targeted therapeutic small molecules)
5. `oncology_hotspots` (24 actionable driver mutations)
6. `pathways` (7 Reactome cellular signaling cascades)
7. `drug_synergy` (12 validated dual-agent combination regimens)
8. `literature_trials` (15 landmark Phase 2/3 clinical studies)
9. `pharmacogenomics_guidelines` (12 CPIC Level A dosing guidelines)
10. `synthetic_lethal_pairs` (12 DepMap validated genetic dependencies)
11. `clinical_survival_benchmarks` (4 Phase 3 landmark survival cohorts)
12. `liquid_biopsy_clones` (4 ctDNA acquired resistance mechanisms)

---

## 🌟 16 Integrated Multi-Omics Modules

1. **📈 Kaplan-Meier Clinical Survival Analysis (`src/survival_analysis.py`)**: Non-parametric survival estimators, Greenwood SE, and Mantel-Cox Log-Rank test.
2. **🩸 Liquid Biopsy & ctDNA MRD Tracking (`src/liquid_biopsy.py`)**: Longitudinal VAF monitoring, molecular response criteria, and resistance bypass therapy.
3. **💊 Pharmacogenomics (PGx) Safety Engine (`src/pharmacogenomics.py`)**: CPIC Level A dosing rules for DPYD, UGT1A1, and CYP2D6.
4. **🎯 Synthetic Lethality & DepMap Engine (`src/synthetic_lethality.py`)**: CRISPR Chronos essentiality and synthetic lethal drug matching.
5. **🩺 Early Cancer Detection & 100k Risk Engine (`src/cancer_risk_predictor.py`)**: Population cohort matching and NICE 2WW criteria.
6. **🗄️ Relational Bio-Knowledgebase (`src/database.py`)**: SQLite 3 engine with 12 multi-omics tables and JSON exports.
7. **🔬 3D Structural Biology (`src/structural_biology.py`)**: Macromolecular coordinate analysis, catalytic histidines, and heme pockets (`4HHB`).
8. **🧬 Needleman-Wunsch Sequence Alignment (`src/sequence_aligner.py`)**: Dynamic programming with BLOSUM62 matrix for global amino acid alignments.
9. **💊 In Silico Drug Screening & Lipinski Scoring (`src/docking_scorer.py`)**: Lipinski Rule of 5, Veber criteria, QED proxy, and binding free energy estimation ($\Delta G$).
10. **💥 Combination Therapy & Drug Synergy Predictor (`src/synergy_predictor.py`)**: Loewe Combination Index (CI) and Bliss Independence models.
11. **🛡️ Immunotherapy Biomarker & TMB Response Predictor (`src/immunotherapy_response.py`)**: TMB ($\ge 10\ \text{mut/Mb}$), MSI-H/dMMR, and PD-L1 TPS evaluation.
12. **📐 AlphaFold Structural Quality & Disorder (`src/alphafold_quality.py`)**: Per-residue pLDDT confidence bands and Intrinsically Disordered Regions (IDRs).
13. **🕸️ Protein-Protein Interaction (PPI) Interactome (`src/interaction_network.py`)**: STRING-calibrated network topology, centrality hubs, and bottleneck proteins.
14. **🎯 Pan-Cancer Precision Oncology (`src/oncology_profiler.py`)**: Clinically actionable hotspot catalogs for *EGFR, KRAS, BRAF, TP53*.
15. **⚡ Dynamic Live Scientific API Search (`src/live_fetcher.py`)**: Real-time REST integrations for RCSB PDB, PubChem, and PubMed with smart caching.
16. **📑 Publication Report Generator (`src/report_generator.py`)**: Standalone publication-ready HTML reports with embedded charts and metrics.

---

## 🚀 Quick Start

### 1. Installation
```bash
git clone https://github.com/goktugcici/computational-biology-pipeline.git
cd computational-biology-pipeline
pip install -r requirements.txt
```

### 2. Clinical Oncology & Translational CLI Usage

```bash
# 1. Model Phase 3 Kaplan-Meier clinical trial survival (FLAURA Osimertinib):
python src/main.py --survival_trial FLAURA_EGFR

# 2. Track longitudinal ctDNA liquid biopsy response & resistance clones:
python src/main.py --ctdna_baseline 4.8 --ctdna_followup 0.02 --ctdna_mutations EGFR_C797S

# 3. Screen DPYD diplotype for 5-FU / Capecitabine toxicity risk:
python src/main.py --dpyd_alleles *1 *2A

# 4. Query synthetic lethal dependency for BRCA1 or MSI-H:
python src/main.py --synthetic_lethal BRCA1

# 5. Evaluate early cancer risk with 100k cohort matching:
python src/main.py --predict_cancer --patient_age 62 --patient_gender M --patient_symptoms "hemoptysis,chronic_cough,weight_loss" --smoking_pack_years 30.0 --family_history

# 6. Full Pipeline Run across all 16 modules:
python src/main.py --pan_cancer EGFR --synergy Dabrafenib Trametinib --tmb 16.2 --msi MSI-H --db_query EGFR --generate_report
```

### 3. Interactive Web Dashboard
```bash
python -m http.server 8000 --bind 127.0.0.1
```
Open **[http://localhost:8000/dashboard/index.html](http://localhost:8000/dashboard/index.html)** in your browser:
- **📈 Sağkalım (KM)**: Interactive Kaplan-Meier stepped survival curves for Phase 3 trials, median survival comparisons, and Log-Rank statistics.
- **🩸 Sıvı Biyopsi (ctDNA)**: Longitudinal VAF calculator, molecular response classification (mCR/mPR/mSD/mPD), relapse lead time, and resistance clone bypass therapies.
- **💊 PGx & Sentetik Letalite**: Interactive DPYD, UGT1A1, and CYP2D6 diplotype toxicity evaluator, CPIC dosing calculator, and DepMap Chronos synthetic lethality visualizer.
- **🩺 Erken Teşhis (100k)**: Real-time risk gauge, 14 interactive NICE alarm symptom pills, 100,000-patient cohort matching card, and organ-specific progress bars.
- **🗄️ Veri Tabanı (Database Explorer)**: Instant search, category filtering, and JSON export over the 12 SQLite relational tables.
- **🧬 3B Molekül Görselleştirici**: Interactive 3Dmol.js WebGL crystallographic rendering.

---

## 🧪 Automated Testing

Covering all 30 unit tests across Kaplan-Meier estimators, Log-Rank tests, ctDNA MRD response, resistance mutations, pharmacogenomics CPIC rules, synthetic lethality DepMap scores, epidemiology cohort generation, NICE guidelines, multi-organ Bayesian risk models, SQLite queries, structural biology, and sequence alignments:
```bash
python -m unittest tests/test_pipeline.py
# or
pytest tests/
```
Output:
```
..............................
----------------------------------------------------------------------
Ran 30 tests in 0.630s

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
