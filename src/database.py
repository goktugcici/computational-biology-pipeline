"""
BioKnowledgebase - Enterprise Unified Multi-Omics Relational Database.
Implements SQLite schema management, multi-omics dataset seeding,
and high-performance query interfaces for genes, variants, 3D structures,
pharmacology compounds, oncology hotspots, pathways, and clinical trials.
"""

import os
import json
import sqlite3
from typing import Dict, List, Any, Optional, Tuple

DEFAULT_DB_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "bio_knowledgebase.db"))

class BioKnowledgebase:
    """
    Unified relational knowledgebase engine for computational biology
    and precision oncology multi-omics analytics.
    """

    def __init__(self, db_path: str = DEFAULT_DB_PATH):
        self.db_path = db_path
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        self.conn: Optional[sqlite3.Connection] = None
        self._init_db()

    def _get_connection(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def _init_db(self):
        """Creates table schemas and indexes if they do not exist."""
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 1. Genes Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS genes (
                gene_symbol TEXT PRIMARY KEY,
                full_name TEXT NOT NULL,
                entrez_id INTEGER,
                chromosome TEXT,
                cytoband TEXT,
                molecular_function TEXT,
                uniprot_id TEXT,
                canonical_transcript TEXT
            );
            """)

            # 2. Variants Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_symbol TEXT NOT NULL,
                rs_id TEXT,
                hgvs_c TEXT,
                hgvs_p TEXT,
                chromosome TEXT,
                position INTEGER,
                ref_allele TEXT,
                alt_allele TEXT,
                clinical_significance TEXT NOT NULL,
                disease_name TEXT,
                clinvar_id TEXT,
                review_status TEXT,
                FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
            );
            """)
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene_symbol);")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_variants_sig ON variants(clinical_significance);")

            # 3. Structures Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS structures (
                pdb_id TEXT PRIMARY KEY,
                title TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                resolution_angstrom REAL,
                experimental_method TEXT,
                chains TEXT,
                bound_ligands TEXT,
                active_site_residues TEXT,
                mean_plddt REAL,
                deposition_year INTEGER,
                FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
            );
            """)
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_structures_gene ON structures(gene_symbol);")

            # 4. Pharmacology Compounds Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS compounds (
                compound_id TEXT PRIMARY KEY,
                common_name TEXT NOT NULL,
                pubchem_cid INTEGER,
                chembl_id TEXT,
                target_gene TEXT,
                primary_indication TEXT,
                molecular_formula TEXT,
                molecular_weight REAL,
                xlogp REAL,
                hbd_count INTEGER,
                hba_count INTEGER,
                tpsa REAL,
                rotatable_bonds INTEGER,
                lipinski_pass INTEGER,
                qed_score REAL,
                fda_status TEXT,
                FOREIGN KEY (target_gene) REFERENCES genes(gene_symbol)
            );
            """)
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_compounds_target ON compounds(target_gene);")

            # 5. Oncology Hotspots Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS oncology_hotspots (
                hotspot_id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_symbol TEXT NOT NULL,
                variant_nomenclature TEXT NOT NULL,
                exon TEXT,
                mutation_type TEXT,
                cancer_types TEXT,
                therapeutic_response TEXT,
                resistance_mechanism TEXT,
                FOREIGN KEY (gene_symbol) REFERENCES genes(gene_symbol)
            );
            """)
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_hotspots_gene ON oncology_hotspots(gene_symbol);")

            # 6. Biological Pathways Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS pathways (
                pathway_id TEXT PRIMARY KEY,
                pathway_name TEXT NOT NULL,
                species TEXT NOT NULL,
                category TEXT,
                key_genes TEXT,
                diagram_url TEXT
            );
            """)

            # 7. Drug Synergy Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS drug_synergy (
                synergy_id INTEGER PRIMARY KEY AUTOINCREMENT,
                drug_a TEXT NOT NULL,
                target_a TEXT NOT NULL,
                drug_b TEXT NOT NULL,
                target_b TEXT NOT NULL,
                indication TEXT,
                combination_index REAL,
                bliss_score REAL,
                mechanistic_class TEXT,
                clinical_trial TEXT
            );
            """)

            # 8. Literature & Clinical Trials Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS literature_trials (
                trial_id INTEGER PRIMARY KEY AUTOINCREMENT,
                identifier TEXT NOT NULL,
                title TEXT NOT NULL,
                journal_or_source TEXT,
                publication_year INTEGER,
                primary_gene TEXT,
                therapy_name TEXT,
                clinical_phase TEXT,
                key_finding TEXT
            );
            """)

            # 9. Pharmacogenomics Guidelines Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS pharmacogenomics_guidelines (
                pgx_id INTEGER PRIMARY KEY AUTOINCREMENT,
                gene_symbol TEXT NOT NULL,
                drug_name TEXT NOT NULL,
                diplotype_example TEXT NOT NULL,
                phenotype TEXT NOT NULL,
                toxicity_risk TEXT NOT NULL,
                dosing_recommendation TEXT NOT NULL,
                guideline_source TEXT NOT NULL
            );
            """)

            # 10. Synthetic Lethal Pairs Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS synthetic_lethal_pairs (
                sl_id INTEGER PRIMARY KEY AUTOINCREMENT,
                biomarker_gene TEXT NOT NULL,
                biomarker_state TEXT NOT NULL,
                synthetic_lethal_target TEXT NOT NULL,
                depmap_chronos REAL NOT NULL,
                vulnerability_score REAL NOT NULL,
                approved_therapy TEXT NOT NULL,
                clinical_status TEXT NOT NULL
            );
            """)

            # 11. Clinical Survival Benchmarks Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS clinical_survival_benchmarks (
                trial_key TEXT PRIMARY KEY,
                trial_name TEXT NOT NULL,
                target_gene TEXT NOT NULL,
                indication TEXT NOT NULL,
                endpoint TEXT NOT NULL,
                experimental_arm TEXT NOT NULL,
                experimental_n INTEGER NOT NULL,
                experimental_median_months REAL NOT NULL,
                control_arm TEXT NOT NULL,
                control_n INTEGER NOT NULL,
                control_median_months REAL NOT NULL,
                hazard_ratio REAL NOT NULL,
                hr_ci_lower REAL NOT NULL,
                hr_ci_upper REAL NOT NULL,
                p_value REAL NOT NULL,
                clinical_significance TEXT NOT NULL
            );
            """)

            # 12. Liquid Biopsy Clones & Acquired Resistance Table
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS liquid_biopsy_clones (
                clone_id INTEGER PRIMARY KEY AUTOINCREMENT,
                mutation_key TEXT NOT NULL UNIQUE,
                gene_symbol TEXT NOT NULL,
                target_therapy TEXT NOT NULL,
                mutation_name TEXT NOT NULL,
                resistance_mechanism TEXT NOT NULL,
                actionable_bypass TEXT NOT NULL
            );
            """)

            conn.commit()

        # Seed data if empty
        self.seed_if_empty()

    def seed_if_empty(self):
        """Seeds the database with extensive multi-omics datasets if empty."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM genes;")
            g_count = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(*) FROM pharmacogenomics_guidelines;")
            p_count = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(*) FROM clinical_survival_benchmarks;")
            c_count = cursor.fetchone()[0]
            if g_count > 0 and p_count > 0 and c_count > 0:
                return # Already fully seeded

        self.seed_comprehensive_data()

    def seed_comprehensive_data(self):
        """Populates comprehensive curated datasets across all tables."""
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 1. Genes
            genes_data = [
                ("EGFR", "Epidermal Growth Factor Receptor", 1956, "7", "7p11.2", "Receptor Tyrosine Kinase (Cell Growth)", "P00533", "NM_005228.5"),
                ("KRAS", "KRAS Proto-Oncogene, GTPase", 3845, "12", "12p12.1", "GTPase Signal Transduction Switch", "P01116", "NM_033360.4"),
                ("BRAF", "B-Raf Proto-Oncogene, Serine/Threonine Kinase", 673, "7", "7q34", "MAPK/ERK Signaling Cascade Kinase", "P15056", "NM_004333.6"),
                ("TP53", "Tumor Protein P53", 7157, "17", "17p13.1", "DNA Damage Response Transcription Factor", "P04637", "NM_000546.6"),
                ("HBB", "Hemoglobin Subunit Beta", 3043, "11", "11p15.4", "Oxygen Transport & Heme Binding", "P68871", "NM_000518.5"),
                ("BRCA1", "BRCA1 DNA Repair Associated", 672, "17", "17q21.31", "Homologous Recombination DNA Repair", "P38398", "NM_007294.4"),
                ("BRCA2", "BRCA2 DNA Repair Associated", 675, "13", "13q13.1", "Rad51 Recombinase Recruitment", "P51587", "NM_000059.4"),
                ("MET", "MET Proto-Oncogene, Receptor Tyrosine Kinase", 4233, "7", "7q31.2", "Hepatocyte Growth Factor Receptor", "P08581", "NM_000245.4"),
                ("ALK", "ALK Receptor Tyrosine Kinase", 238, "2", "2p23.2", "Neurogenesis & Oncogenic Fusion Kinase", "Q9UM73", "NM_004304.5"),
                ("PIK3CA", "Phosphatidylinositol-4,5-Bisphosphate 3-Kinase Catalytic Subunit Alpha", 5290, "3", "3q26.32", "PI3K/AKT/mTOR Survival Signaling", "P42336", "NM_006218.4")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO genes (gene_symbol, full_name, entrez_id, chromosome, cytoband, molecular_function, uniprot_id, canonical_transcript)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?);
            """, genes_data)

            # 2. Variants
            variants_data = [
                ("HBB", "rs334", "c.20A>T", "p.Glu7Val", "11", 5227002, "A", "T", "Pathogenic", "Sickle Cell Anemia / HbS", "15333", "criteria provided, single submitter"),
                ("HBB", "rs33930165", "c.20A>G", "p.Glu7Gly", "11", 5227002, "A", "G", "Pathogenic", "Hemoglobin G-Makassar", "15334", "criteria provided, single submitter"),
                ("HBB", "rs33950507", "c.19G>A", "p.Glu7Lys", "11", 5227001, "G", "A", "Pathogenic", "Hemoglobin C Disease", "15335", "reviewed by expert panel"),
                ("HBB", "rs35004220", "c.92+1G>A", "Splice Donor", "11", 5226929, "G", "A", "Pathogenic", "Beta-Thalassemia Major", "15340", "criteria provided, multiple submitters"),
                ("EGFR", "rs121434568", "c.2573T>G", "p.Leu858Arg", "7", 55191822, "T", "G", "Pathogenic", "Non-Small Cell Lung Cancer (TKI Sensitive)", "16609", "reviewed by expert panel"),
                ("EGFR", "rs121434569", "c.2369C>T", "p.Thr790M", "7", 55181378, "C", "T", "Pathogenic", "Acquired 1st Gen EGFR TKI Resistance", "16610", "reviewed by expert panel"),
                ("KRAS", "rs121913529", "c.34G>T", "p.Gly12Cys", "12", 25245350, "G", "T", "Pathogenic", "Lung & Colorectal Carcinoma", "12582", "reviewed by expert panel"),
                ("KRAS", "rs121913530", "c.35G>A", "p.Gly12Asp", "12", 25245347, "G", "A", "Pathogenic", "Pancreatic Ductal Adenocarcinoma", "12583", "reviewed by expert panel"),
                ("BRAF", "rs113488022", "c.1799T>A", "p.Val600Glu", "7", 140753336, "T", "A", "Pathogenic", "Malignant Melanoma & Papillary Thyroid", "13961", "reviewed by expert panel"),
                ("TP53", "rs28934578", "c.524G>A", "p.Arg175His", "17", 7675088, "G", "A", "Pathogenic", "Li-Fraumeni Syndrome / Diverse Carcinomas", "12374", "reviewed by expert panel"),
                ("TP53", "rs11540652", "c.743G>A", "p.Arg248Gln", "17", 7674220, "G", "A", "Pathogenic", "Li-Fraumeni Syndrome", "12375", "reviewed by expert panel"),
                ("BRCA1", "rs80357906", "c.68_69delAG", "p.Glu23fs", "17", 43124031, "AG", "-", "Pathogenic", "Hereditary Breast & Ovarian Cancer", "17660", "reviewed by expert panel"),
                ("BRCA2", "rs80359550", "c.5946delT", "p.Ser1982fs", "13", 32340300, "T", "-", "Pathogenic", "Hereditary Breast & Ovarian Cancer", "9325", "reviewed by expert panel"),
                ("BRCA2", "rs28897700", "c.198A>G", "p.Lys66Lys", "13", 32316499, "A", "G", "Benign", "Synonymous Control Variant", "51234", "reviewed by expert panel")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO variants (gene_symbol, rs_id, hgvs_c, hgvs_p, chromosome, position, ref_allele, alt_allele, clinical_significance, disease_name, clinvar_id, review_status)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """, variants_data)

            # 3. Structures
            structures_data = [
                ("4HHB", "Human Deoxyhemoglobin at 1.74 A resolution", "HBB", 1.74, "X-RAY DIFFRACTION", "A,B,C,D", "HEM, PHO, OXY", "His87(A), His58(A), Glu6(B)", 94.8, 1984),
                ("1TSR", "Crystal Structure of p53 Tumor Suppressor-DNA Complex", "TP53", 2.20, "X-RAY DIFFRACTION", "A,B,C", "ZN", "Cys176, His179, Cys238, Cys242", 92.5, 1994),
                ("2J6M", "Structure of EGFR Kinase Domain in Complex with Gefitinib", "EGFR", 2.80, "X-RAY DIFFRACTION", "A,B", "IRE (Gefitinib)", "Met793, Cys797, Leu858", 93.8, 2007),
                ("6OIM", "KRAS G12C bound to covalent inhibitor ARS-1620", "KRAS", 1.39, "X-RAY DIFFRACTION", "A", "GDP, ARS", "Cys12, Asp54, Gly60", 98.2, 2019),
                ("6LU7", "The crystal structure of COVID-19 main protease in complex with an inhibitor N3", "TP53", 2.16, "X-RAY DIFFRACTION", "A", "N3", "His41, Cys145", 89.4, 2020),
                ("4ZAU", "Crystal Structure of BRAF V600E Kinase Domain Complexed with Dabrafenib", "BRAF", 2.50, "X-RAY DIFFRACTION", "A", "13D (Dabrafenib)", "Cys532, Asp594, Glu600", 94.2, 2015),
                ("6N4G", "Structure of human MET kinase domain with Savolitinib", "MET", 2.15, "X-RAY DIFFRACTION", "A", "SAV", "Met1160, Tyr1230", 91.0, 2019)
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO structures (pdb_id, title, gene_symbol, resolution_angstrom, experimental_method, chains, bound_ligands, active_site_residues, mean_plddt, deposition_year)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """, structures_data)

            # 4. Pharmacology Compounds
            compounds_data = [
                ("OSIMERTINIB", "Osimertinib (Tagrisso)", 71496458, "CHEMBL3353410", "EGFR", "EGFR T790M/L858R NSCLC", "C28H33N7O2", 499.6, 3.7, 1, 7, 75.3, 9, 1, 82.4, "Approved"),
                ("GEFITINIB", "Gefitinib (Iressa)", 123631, "CHEMBL939", "EGFR", "Metastatic EGFR mutant NSCLC", "C22H24ClFN4O3", 446.9, 4.1, 1, 6, 68.7, 7, 1, 86.1, "Approved"),
                ("SOTORASIB", "Sotorasib (Lumakras)", 137278711, "CHEMBL4475865", "KRAS", "KRAS G12C NSCLC & CRC", "C30H30F2N6O3", 560.6, 2.8, 0, 7, 91.5, 5, 1, 70.8, "Approved"),
                ("ADAGRASIB", "Adagrasib (Krazati)", 137450376, "CHEMBL4650047", "KRAS", "KRAS G12C Advanced NSCLC", "C32H35Cl2FN6O2", 604.1, 3.9, 0, 6, 65.2, 5, 1, 68.4, "Approved"),
                ("DABRAFENIB", "Dabrafenib (Tafinlar)", 44462760, "CHEMBL1963831", "BRAF", "BRAF V600E Melanoma / NSCLC", "C23H20F3N5O2S2", 519.6, 4.3, 2, 6, 105.0, 6, 1, 72.5, "Approved"),
                ("VOXELOTOR", "Voxelotor (Oxbryta)", 71602803, "CHEMBL3707348", "HBB", "Sickle Cell Disease HbS Polymerization", "C19H19N3O3", 337.4, 2.4, 1, 5, 62.1, 5, 1, 89.2, "Approved"),
                ("METFORMIN", "Metformin", 4091, "CHEMBL1431", "PIK3CA", "Type 2 Diabetes & AMPK Activation", "C4H11N5", 129.16, -1.4, 3, 2, 88.0, 0, 1, 80.0, "Approved"),
                ("IMATINIB", "Imatinib (Gleevec)", 5291, "CHEMBL941", "EGFR", "Chronic Myeloid Leukemia & GIST", "C29H31N7O", 493.6, 3.5, 2, 7, 86.3, 7, 1, 81.9, "Approved"),
                ("SAVOLITINIB", "Savolitinib (Orpathys)", 53342379, "CHEMBL3622434", "MET", "MET-Amplified EGFR Resistance", "C17H15N9", 345.36, 1.8, 1, 8, 88.3, 3, 1, 78.4, "Approved"),
                ("TRAMETINIB", "Trametinib (Mekinist)", 11707110, "CHEMBL2048451", "BRAF", "MEK1/MEK2 Kinase Inhibition", "C26H23FIN5O4", 615.4, 3.2, 1, 7, 98.7, 4, 1, 67.2, "Approved"),
                ("ASPIRIN", "Acetylsalicylic Acid (Aspirin)", 2244, "CHEMBL25", "HBB", "Antiplatelet & Anti-inflammatory", "C9H8O4", 180.16, 1.2, 1, 3, 63.6, 3, 1, 85.4, "Approved")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO compounds (compound_id, common_name, pubchem_cid, chembl_id, target_gene, primary_indication, molecular_formula, molecular_weight, xlogp, hbd_count, hba_count, tpsa, rotatable_bonds, lipinski_pass, qed_score, fda_status)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """, compounds_data)

            # 5. Oncology Hotspots
            hotspots_data = [
                ("EGFR", "p.L858R", "Exon 21", "Missense Substitution", "Non-Small Cell Lung Cancer", "1st/3rd gen TKI Duyarlı (Osimertinib, Gefitinib)", "T790M veya MET amplifikasyonu"),
                ("EGFR", "p.T790M", "Exon 20", "Gatekeeper Missense", "NSCLC (Secondary Resistance)", "3. Nesil Osimertinib Duyarlı", "C797S mutasyonu"),
                ("EGFR", "Exon 19 del", "Exon 19", "In-frame Deletion", "NSCLC", "EGFR TKI Dramatik Klinik Yanıt", "T790M ikincil mutasyonu"),
                ("KRAS", "p.G12C", "Exon 2", "Codon 12 Missense", "NSCLC, Kolorektal, Pankreas", "Sotorasib ve Adagrasib Doğrudan Hedefleme", "Geri bildirim EGFR veya RTK aktivasyonu"),
                ("KRAS", "p.G12D", "Exon 2", "Codon 12 Missense", "Pankreas Duktal Adenokarsinomu", "Geliştirilmekte Olan Pan-KRAS / Degrader İlaçlar", "GTPaz hidroliz kilitlenmesi"),
                ("BRAF", "p.V600E", "Exon 15", "Activation Loop Missense", "Malign Melanom, Kolorektal, Tiroid", "Dabrafenib + Trametinib Kombinasyonu", "MEK baypas ve CRAF dimerizasyonu"),
                ("TP53", "p.R175H", "Exon 5", "Zinc Binding Structural Defect", "Kolorektal, Akciğer, Meme", "Eprenetapopt (APR-246) Yapı Düzelticisi", "Dominant-negatif tetramer disfonksiyonu"),
                ("TP53", "p.R248Q", "Exon 7", "DNA Contact Mutation", "Glioblastoma, Meme, Over", "WEE1/ATR Sentetik Ölümcüllük Adayı", "DNA temas kaybı")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO oncology_hotspots (gene_symbol, variant_nomenclature, exon, mutation_type, cancer_types, therapeutic_response, resistance_mechanism)
                VALUES (?, ?, ?, ?, ?, ?, ?);
            """, hotspots_data)

            # 6. Biological Pathways
            pathways_data = [
                ("R-HSA-1247673", "Erythrocytes take up oxygen and release carbon dioxide", "Homo sapiens", "Hemostasis & Transport", "HBB, HBA1, HBA2, BPGM", "https://reactome.org/PathwayBrowser/#/R-HSA-1247673"),
                ("R-HSA-1480926", "O2/CO2 exchange in erythrocytes", "Homo sapiens", "Cellular Transport", "HBB, HBA1, CA1, CA2", "https://reactome.org/PathwayBrowser/#/R-HSA-1480926"),
                ("R-HSA-9927020", "Heme assimilation", "Homo sapiens", "Metabolism", "HBB, FECH, ALAS2, HMBS", "https://reactome.org/PathwayBrowser/#/R-HSA-9927020"),
                ("R-HSA-177929", "Signaling by EGFR", "Homo sapiens", "Signal Transduction", "EGFR, ERBB2, GRB2, SOS1, SHC1", "https://reactome.org/PathwayBrowser/#/R-HSA-177929"),
                ("R-HSA-5683057", "MAPK family signaling cascades", "Homo sapiens", "Signal Transduction", "KRAS, BRAF, RAF1, MAP2K1, MAPK1", "https://reactome.org/PathwayBrowser/#/R-HSA-5683057"),
                ("R-HSA-69620", "Cell Cycle Checkpoints", "Homo sapiens", "Cell Cycle", "TP53, CDKN1A, ATM, CHEK2, MDM2", "https://reactome.org/PathwayBrowser/#/R-HSA-69620"),
                ("R-HSA-5693532", "DNA Double-Strand Break Repair", "Homo sapiens", "DNA Repair", "BRCA1, BRCA2, RAD51, PALB2, ATM", "https://reactome.org/PathwayBrowser/#/R-HSA-5693532")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO pathways (pathway_id, pathway_name, species, category, key_genes, diagram_url)
                VALUES (?, ?, ?, ?, ?, ?);
            """, pathways_data)

            # 7. Drug Synergy
            synergy_data = [
                ("Dabrafenib", "BRAF V600E", "Trametinib", "MEK1/2", "BRAF V600E Melanom & NSCLC", 0.42, 28.5, "Vertical MAPK Pathway Blockade", "COMBI-d / COMBI-v Phase 3"),
                ("Osimertinib", "EGFR T790M", "Savolitinib", "MET Receptor", "MET-Amplified EGFR Resistance", 0.51, 24.0, "Parallel Receptor Bypass Inhibition", "TATTON / SAVANNAH Phase 2"),
                ("Sotorasib", "KRAS G12C", "Panitumumab", "EGFR mAb", "KRAS G12C Colorectal Carcinoma", 0.48, 26.2, "Upstream RTK Feedback Loop Suppression", "CodeBreaK 300 Phase 3"),
                ("Pembrolizumab", "PD-1", "Lenvatinib", "VEGFR/FGFR", "Renal Cell Carcinoma & Endometrial", 0.58, 19.8, "Immuno-Angiogenic Microenvironment", "CLEAR Phase 3")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO drug_synergy (drug_a, target_a, drug_b, target_b, indication, combination_index, bliss_score, mechanistic_class, clinical_trial)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);
            """, synergy_data)

            # 8. Literature & Clinical Trials
            trials_data = [
                ("PMID: 38661449", "Exagamglogene Autotemcel for Severe Sickle Cell Disease", "New England Journal of Medicine", 2024, "HBB", "Exa-cel (Casgevy)", "Phase 3", "CRISPR-Cas9 editing of BCL11A enhancer produces durable elimination of vaso-occlusive crises in 97% of patients."),
                ("PMID: 29151359", "Osimertinib in Untreated EGFR-Mutated Advanced Non-Small-Cell Lung Cancer", "New England Journal of Medicine", 2018, "EGFR", "Osimertinib (Tagrisso)", "Phase 3 (FLAURA)", "Significantly superior progression-free survival (18.9 vs 10.2 months) and central nervous system efficacy compared to standard EGFR TKIs."),
                ("PMID: 34096690", "Overall Survival with Sotorasib in Advanced KRAS G12C-Mutated Non-Small-Cell Lung Cancer", "New England Journal of Medicine", 2021, "KRAS", "Sotorasib (Lumakras)", "Phase 2 (CodeBreaK 100)", "First targeted therapeutic agent demonstrating confirmed objective response in KRAS G12C mutant solid tumors."),
                ("PMID: 25399551", "Combined BRAF and MEK Inhibition in Melanoma with BRAF V600 Mutations", "New England Journal of Medicine", 2014, "BRAF", "Dabrafenib + Trametinib", "Phase 3 (COMBI-v)", "Combination therapy significantly improved overall survival compared to vemurafenib monotherapy alone."),
                ("PMID: 31922372", "Voxelotor for the Treatment of Sickle Cell Disease", "New England Journal of Medicine", 2019, "HBB", "Voxelotor (Oxbryta)", "Phase 3 (HOPE)", "Dose-dependent increase in hemoglobin and reduction in hemolytic markers without causing tissue hypoxia.")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO literature_trials (identifier, title, journal_or_source, publication_year, primary_gene, therapy_name, clinical_phase, key_finding)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?);
            """, trials_data)

            # 9. Pharmacogenomics (PGx) Guidelines
            pgx_data = [
                ("DPYD", "Fluorouracil (5-FU) / Capecitabine", "*1/*2A (c.1905+1G>A)", "Intermediate Metabolizer (AS: 1.0)", "Moderate/High Toxicity", "%50 Başlangıç Dozu Azaltımı Önerilir. TDM takibi yapınız.", "CPIC Level A / PharmGKB"),
                ("DPYD", "Fluorouracil (5-FU) / Capecitabine", "*2A/*2A (c.1905+1G>A)", "Poor Metabolizer (AS: 0.0)", "Fatal / Lethal Toxicity Risk", "KONTRENDİKE! 5-FU ve Kapesitabin kullanımından tamamen kaçınınız.", "CPIC Level A / PharmGKB"),
                ("UGT1A1", "Irinotecan (Camptosar)", "*28/*28 ((TA)7 repeat)", "Poor Metabolizer", "Severe Grade 4 Neutropenia / Diarrhea", "İrinotekan başlangıç dozunda %25-30 azaltım önerilir.", "CPIC Level A / FDA"),
                ("UGT1A1", "Irinotecan (Camptosar)", "*1/*28 ((TA)6/(TA)7)", "Intermediate Metabolizer", "Moderate Toxicity Risk", "Standart başlangıç dozu uygulanabilir; haftalık tam kan sayımı izlenmelidir.", "CPIC Level A"),
                ("CYP2D6", "Tamoxifen (Nolvadex)", "*4/*4 (Splice Defect)", "Poor Metabolizer", "Therapeutic Resistance / Subtherapeutic Endoxifen", "Tamoksifen aktif metabolit üretemez. Postmenopozal kadınlarda Aromataz İnhibitörüne geçiniz.", "CPIC Level A"),
                ("TPMT", "6-Mercaptopurine / Azathioprine", "*1/*3A (c.460G>A, c.719A>G)", "Intermediate Metabolizer", "High Myelosuppression", "%30-50 doz azaltımı ve haftalık lökosit takibi önerilir.", "CPIC Level A")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO pharmacogenomics_guidelines (gene_symbol, drug_name, diplotype_example, phenotype, toxicity_risk, dosing_recommendation, guideline_source)
                VALUES (?, ?, ?, ?, ?, ?, ?);
            """, pgx_data)

            # 10. Synthetic Lethal Pairs (DepMap Dependency)
            sl_data = [
                ("BRCA1", "HRD+ / Inactivating Mutation", "PARP1", -1.24, 98.0, "Olaparib / Talazoparib", "FDA Onaylı"),
                ("BRCA2", "HRD+ / Truncating Mutation", "PARP1", -1.18, 96.0, "Olaparib / Rucaparib", "FDA Onaylı"),
                ("MLH1 / MSH2", "dMMR / MSI-High", "WRN", -1.45, 95.0, "HRO761 / RO7589831", "Faz 1/2 Klinik"),
                ("MTAP", "Homozygous 9p21 Deletion", "PRMT5 / MAT2A", -0.92, 88.0, "MRTX1719 / AG-270", "Faz 1/2 Klinik"),
                ("ARID1A", "Truncating SWI/SNF Loss", "EZH2", -0.85, 84.0, "Tazemetostat (Tazverik)", "FDA Onaylı"),
                ("KRAS", "G12C / G12D / G12V Mutation", "SHP2 / SOS1", -0.78, 82.0, "TNO155 + Sotorasib / BI-1701963", "Faz 1/2 Klinik")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO synthetic_lethal_pairs (biomarker_gene, biomarker_state, synthetic_lethal_target, depmap_chronos, vulnerability_score, approved_therapy, clinical_status)
                VALUES (?, ?, ?, ?, ?, ?, ?);
            """, sl_data)

            # 11. Clinical Survival Benchmarks (Phase 3 Landmark Trials)
            survival_data = [
                ("FLAURA_EGFR", "FLAURA (NEJM 2018 / 2020)", "EGFR (Exon 19 del / L858R)", "Untreated Advanced EGFR-Mutated NSCLC", "Overall Survival (OS)", "Osimertinib (Tagrisso 80mg)", 279, 38.6, "Standard EGFR-TKI (Gefitinib / Erlotinib)", 277, 31.8, 0.80, 0.64, 1.00, 0.046, "Osimertinib, standart 1. jenerasyon TKI lere kıyasla genel sağkalımda istatistiksel ve klinik olarak anlamlı 6.8 aylık uzama sağlar."),
                ("KEYNOTE_042_TMB_HIGH", "KEYNOTE-042 (Lancet 2019 / NEJM)", "TMB-High (TMB >= 10 mut/Mb) & PD-L1 TPS >= 50%", "First-Line Advanced Non-Small-Cell Lung Cancer", "Overall Survival (OS)", "Pembrolizumab Monoterapisi (Keytruda)", 120, 20.0, "Platin Bazlı Kemoterapi (Karboplatin + Pemetreksed)", 120, 12.2, 0.69, 0.56, 0.85, 0.0003, "TMB-High ve PD-L1 pozitif NSCLC hastalarında immünoterapi, kemoterapiye kıyasla ölüm riskini %31 oranında azaltır."),
                ("OLYMPIAD_BRCA", "OlympiAD (NEJM 2017)", "Germline BRCA1 / BRCA2 Mutation (HRD+)", "gBRCAm HER2-Negative Metastatic Breast Cancer", "Progression-Free Survival (PFS)", "Olaparib Tablet (Lynparza 300mg BID)", 205, 7.0, "Hekim Seçimi Standart Kemoterapi (Kapesitabin/Vinorelbin)", 97, 4.2, 0.58, 0.43, 0.80, 0.0009, "PARP inhibitörü Olaparib, gBRCA mutant metastatik meme kanserinde progresyonsuz sağkalımı belirgin şekilde uzatarak progresyon riskini %42 düşürür."),
                ("CODEBREAK_200_KRAS", "CodeBreaK 200 (Lancet 2023)", "KRAS p.G12C", "Pretreated KRAS G12C Advanced NSCLC", "Progression-Free Survival (PFS)", "Sotorasib (Lumakras 960mg QD)", 171, 5.6, "Docetaxel Kemoterapisi", 174, 4.5, 0.66, 0.51, 0.86, 0.002, "Doğrudan KRAS G12C inhibitörü Sotorasib, önceden tedavi görmüş NSCLC hastalarında dosetaksele kıyasla progresyonsuz sağkalımı anlamlı şekilde iyileştirir.")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO clinical_survival_benchmarks (trial_key, trial_name, target_gene, indication, endpoint, experimental_arm, experimental_n, experimental_median_months, control_arm, control_n, control_median_months, hazard_ratio, hr_ci_lower, hr_ci_upper, p_value, clinical_significance)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """, survival_data)

            # 12. Liquid Biopsy Clones & Acquired Resistance Mechanisms
            ctdna_clones_data = [
                ("EGFR_C797S", "EGFR", "Osimertinib (3rd-gen TKI)", "p.Cys797Ser (c.2389T>A)", "Kovalent bağlanma cebindeki sistein 797 kalıntısının serin ile değişmesi; Osimertinibin kovalent bağ kurmasını engeller.", "Allosterik 4. jenerasyon EGFR TKI (örn. BLU-945, BDTX-1535) veya antikor-ilaç konjugatı (Patritumab Deruxtecan - HER3-DXd)."),
                ("MET_AMPLIFICATION", "MET", "Osimertinib / Gefitinib", "MET Gene Copy Number Gain (Copy Number >= 5)", "EGFR yolağını baypas ederek ERBB3 fosforilasyonunu ve downstream PI3K/AKT sinyalini devam ettirir.", "Kombinasyon: Osimertinib + MET inhibitörü Savolitinib (SAVANNAH Faz 2 çalışması) veya Tepotinib."),
                ("KRAS_Y96D", "KRAS", "Sotorasib / Adagrasib (KRAS G12C İnhibitörleri)", "p.Tyr96Asp (Switch-II cebi direnci)", "Switch-II allosterik cebindeki Y96 mutasyonu, kovalent inhibitörlerin bağlanmasını sterik olarak bloke eder.", "Pan-KRAS inhibitörleri (RMC-6236, Darovasertib) veya SOS1/SHP2 inhibitör kombinasyonları."),
                ("BRCA_REVERSION", "BRCA1 / BRCA2", "Olaparib / Talazoparib (PARP İnhibitörleri)", "Secondary in-frame restore / Reversion mutation", "İkincil mutasyon okuma çerçevesini yeniden düzelterek homolog rekombinasyon (HR) kapasitesini geri kazandırır.", "ATR inhibitörleri (Ceralasertib), WEE1 inhibitörleri (Adavosertib) veya Platin bazlı ikili kemoterapi.")
            ]
            cursor.executemany("""
                INSERT OR REPLACE INTO liquid_biopsy_clones (mutation_key, gene_symbol, target_therapy, mutation_name, resistance_mechanism, actionable_bypass)
                VALUES (?, ?, ?, ?, ?, ?);
            """, ctdna_clones_data)

            conn.commit()

    # Query APIs
    def query_genes(self, chromosome: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if chromosome:
                cursor.execute("SELECT * FROM genes WHERE chromosome = ? ORDER BY gene_symbol ASC;", (str(chromosome),))
            else:
                cursor.execute("SELECT * FROM genes ORDER BY gene_symbol ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_variants(self, gene: Optional[str] = None, significance: Optional[str] = None, limit: int = 50) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            query = "SELECT * FROM variants WHERE 1=1"
            params = []
            if gene:
                query += " AND gene_symbol = ?"
                params.append(gene.upper())
            if significance:
                query += " AND clinical_significance LIKE ?"
                params.append(f"%{significance}%")
            query += " ORDER BY position ASC LIMIT ?;"
            params.append(limit)
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]

    def query_structures(self, gene: Optional[str] = None, max_resolution: Optional[float] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            query = "SELECT * FROM structures WHERE 1=1"
            params = []
            if gene:
                query += " AND gene_symbol = ?"
                params.append(gene.upper())
            if max_resolution is not None:
                query += " AND resolution_angstrom <= ?"
                params.append(max_resolution)
            query += " ORDER BY resolution_angstrom ASC;"
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]

    def query_compounds(self, target: Optional[str] = None, max_mw: Optional[float] = None, lipinski_only: bool = False) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            query = "SELECT * FROM compounds WHERE 1=1"
            params = []
            if target:
                query += " AND target_gene = ?"
                params.append(target.upper())
            if max_mw is not None:
                query += " AND molecular_weight <= ?"
                params.append(max_mw)
            if lipinski_only:
                query += " AND lipinski_pass = 1"
            query += " ORDER BY qed_score DESC;"
            cursor.execute(query, params)
            return [dict(row) for row in cursor.fetchall()]

    def query_hotspots(self, gene: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if gene:
                cursor.execute("SELECT * FROM oncology_hotspots WHERE gene_symbol = ? ORDER BY hotspot_id ASC;", (gene.upper(),))
            else:
                cursor.execute("SELECT * FROM oncology_hotspots ORDER BY gene_symbol ASC, hotspot_id ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_synergies(self, drug: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if drug:
                cursor.execute("""
                    SELECT * FROM drug_synergy 
                    WHERE drug_a LIKE ? OR drug_b LIKE ? 
                    ORDER BY combination_index ASC;
                """, (f"%{drug}%", f"%{drug}%"))
            else:
                cursor.execute("SELECT * FROM drug_synergy ORDER BY combination_index ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_pathways(self, gene: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if gene:
                cursor.execute("SELECT * FROM pathways WHERE key_genes LIKE ? ORDER BY pathway_name ASC;", (f"%{gene.upper()}%",))
            else:
                cursor.execute("SELECT * FROM pathways ORDER BY pathway_name ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_literature(self, gene: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if gene:
                cursor.execute("SELECT * FROM literature_trials WHERE primary_gene = ? ORDER BY publication_year DESC;", (gene.upper(),))
            else:
                cursor.execute("SELECT * FROM literature_trials ORDER BY publication_year DESC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_pgx(self, gene: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if gene:
                cursor.execute("SELECT * FROM pharmacogenomics_guidelines WHERE gene_symbol = ?;", (gene.upper(),))
            else:
                cursor.execute("SELECT * FROM pharmacogenomics_guidelines ORDER BY gene_symbol ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_synthetic_lethality(self, biomarker: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if biomarker:
                cursor.execute("SELECT * FROM synthetic_lethal_pairs WHERE biomarker_gene LIKE ? OR synthetic_lethal_target LIKE ?;", (f"%{biomarker.upper()}%", f"%{biomarker.upper()}%"))
            else:
                cursor.execute("SELECT * FROM synthetic_lethal_pairs ORDER BY vulnerability_score DESC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_survival_benchmarks(self, trial_key: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if trial_key:
                cursor.execute("SELECT * FROM clinical_survival_benchmarks WHERE trial_key = ?;", (trial_key,))
            else:
                cursor.execute("SELECT * FROM clinical_survival_benchmarks ORDER BY trial_name ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def query_liquid_biopsy_clones(self, gene: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if gene:
                cursor.execute("SELECT * FROM liquid_biopsy_clones WHERE gene_symbol = ?;", (gene.upper(),))
            else:
                cursor.execute("SELECT * FROM liquid_biopsy_clones ORDER BY gene_symbol ASC, clone_id ASC;")
            return [dict(row) for row in cursor.fetchall()]

    def get_database_stats(self) -> Dict[str, Any]:
        """Calculates global statistical summaries across the entire knowledgebase."""
        with self._get_connection() as conn:
            cursor = conn.cursor()

            stats = {}
            for table in [
                'genes', 'variants', 'structures', 'compounds', 'oncology_hotspots',
                'pathways', 'drug_synergy', 'literature_trials', 'pharmacogenomics_guidelines',
                'synthetic_lethal_pairs', 'clinical_survival_benchmarks', 'liquid_biopsy_clones'
            ]:
                cursor.execute(f"SELECT COUNT(*) FROM {table};")
                stats[f"{table}_count"] = cursor.fetchone()[0]

            # Variant breakdown
            cursor.execute("SELECT clinical_significance, COUNT(*) FROM variants GROUP BY clinical_significance;")
            stats['variants_by_significance'] = {row[0]: row[1] for row in cursor.fetchall()}

            # Average structure resolution
            cursor.execute("SELECT AVG(resolution_angstrom), AVG(mean_plddt) FROM structures;")
            row = cursor.fetchone()
            stats['avg_structure_resolution_angstrom'] = round(row[0], 2) if row[0] else 0.0
            stats['avg_alphafold_plddt'] = round(row[1], 1) if row[1] else 0.0

            # Pharmacology summary
            cursor.execute("SELECT AVG(molecular_weight), AVG(qed_score) FROM compounds;")
            row = cursor.fetchone()
            stats['avg_compound_mw'] = round(row[0], 1) if row[0] else 0.0
            stats['avg_qed_score'] = round(row[1], 1) if row[1] else 0.0

            return stats

    def get_epidemiology_stats(self) -> Dict[str, Any]:
        """Returns statistics on the 100,000+ cancer epidemiology cohort if available."""
        from cancer_risk_predictor import DEFAULT_COHORT_DB_PATH
        if not os.path.exists(DEFAULT_COHORT_DB_PATH):
            return {'total_records': 0, 'status': 'not_generated'}
        try:
            conn = sqlite3.connect(DEFAULT_COHORT_DB_PATH)
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*), SUM(cancer_diagnosed) FROM cancer_cohort;")
            row = cursor.fetchone()
            total_records = row[0] if row else 0
            cancer_cases = row[1] if row else 0
            cursor.execute("SELECT primary_cancer_type, COUNT(*) FROM cancer_cohort WHERE cancer_diagnosed = 1 GROUP BY primary_cancer_type ORDER BY COUNT(*) DESC;")
            breakdown = {r[0]: r[1] for r in cursor.fetchall()}
            conn.close()
            return {
                'total_records': total_records,
                'total_cancer_cases': cancer_cases,
                'overall_cancer_prevalence_percent': round((cancer_cases / total_records) * 100, 2) if total_records > 0 else 0.0,
                'cancer_distribution': breakdown,
                'status': 'active'
            }
        except Exception as e:
            return {'total_records': 0, 'error': str(e)}

    def export_full_json(self, target_file: str) -> str:
        """Exports the complete contents of all relational tables into a JSON document."""
        full_kb = {
            'metadata': {
                'engine': 'BioKnowledgebase SQLite 3',
                'database_file': self.db_path,
                'stats': self.get_database_stats()
            },
            'genes': self.query_genes(),
            'variants': self.query_variants(limit=500),
            'structures': self.query_structures(),
            'compounds': self.query_compounds(),
            'oncology_hotspots': self.query_hotspots(),
            'pathways': self.query_pathways(),
            'drug_synergy': self.query_synergies(),
            'literature_trials': self.query_literature(),
            'pharmacogenomics_guidelines': self.query_pgx(),
            'synthetic_lethal_pairs': self.query_synthetic_lethality(),
            'clinical_survival_benchmarks': self.query_survival_benchmarks(),
            'liquid_biopsy_clones': self.query_liquid_biopsy_clones()
        }

        os.makedirs(os.path.dirname(target_file), exist_ok=True)
        with open(target_file, 'w', encoding='utf-8') as f:
            json.dump(full_kb, f, indent=2, ensure_ascii=False)
        return target_file


if __name__ == '__main__':
    kb = BioKnowledgebase()
    stats = kb.get_database_stats()
    print("=" * 60)
    print("BIO-KNOWLEDGEBASE STATISTICAL OVERVIEW")
    print("=" * 60)
    print(f"Database File: {kb.db_path}")
    for k, v in stats.items():
        print(f"  * {k}: {v}")
    
    print("\nSample Gene Query (EGFR):")
    vars_egfr = kb.query_variants(gene="EGFR")
    for v in vars_egfr:
        print(f"  -> {v['gene_symbol']} {v['hgvs_p']} ({v['rs_id']}) - {v['clinical_significance']}")
