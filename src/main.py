"""
Computational Biology & Science Pipeline - Main Runner (v3.1 Database Edition)
Executes an end-to-end multi-omics analysis spanning structural biology,
genomics, pharmacology, biological pathways, literature intelligence,
pairwise sequence alignments, in silico drug screening, PPI interactomes,
drug combination synergy prediction, immunotherapy response (TMB/MSI),
AlphaFold structural confidence metrics (pLDDT), and unified SQLite database management.
"""

import sys
import json
import argparse
from structural_biology import StructureAnalyzer
from variant_analyzer import VariantAnalyzer
from pharmacology import PharmacologyAnalyzer
from pathway_mapper import PathwayMapper
from literature_miner import LiteratureMiner
from live_fetcher import LiveFetcher
from oncology_profiler import OncologyProfiler
from sequence_aligner import SequenceAligner, KNOWN_VARIANTS
from docking_scorer import DockingScorer
from interaction_network import InteractionNetwork
from synergy_predictor import SynergyPredictor
from immunotherapy_response import ImmunotherapyPredictor
from alphafold_quality import AlphaFoldQualityAnalyzer
from database import BioKnowledgebase
from cancer_risk_predictor import CancerRiskPredictor
from pharmacogenomics import PharmacogenomicsEngine
from synthetic_lethality import SyntheticLethalityEngine
from survival_analysis import SurvivalAnalysisModeler, BENCHMARK_CLINICAL_COHORTS
from liquid_biopsy import LiquidBiopsyEngine, ACQUIRED_RESISTANCE_MUTATIONS
from report_generator import ReportGenerator

def run_pipeline(data_dir: str = "data", live_pdb: str = None, live_drug: str = None, 
                 live_pubmed: str = None, pan_cancer: str = None, align_variant: str = None,
                 screen_drugs_target: str = None, ppi_gene: str = None,
                 synergy_pair: list = None, tmb_val: float = None, msi_status: str = "MSS",
                 pdl1_val: float = 0.0, alphafold_gene: str = None,
                 db_query_gene: str = None, show_db_stats: bool = False,
                 predict_cancer_risk: bool = False, patient_age: int = 50,
                 patient_gender: str = "M", patient_symptoms: str = None,
                 smoking_pack_years: float = 0.0, family_history: bool = False,
                 build_cohort_size: int = None,
                 pgx_check: bool = False, dpyd_alleles: list = None,
                 ugt1a1_diplotype: str = None, cyp2d6_diplotype: str = None,
                 synthetic_lethal_target: str = None,
                 survival_trial: str = "FLAURA_EGFR",
                 ctdna_baseline_vaf: float = 4.8,
                 ctdna_followup_vaf: float = 0.02,
                 ctdna_mutations: list = None,
                 export_db_json_path: str = None, export_json_path: str = None,
                 generate_report: bool = False):
    print("=" * 70)
    print("COMPUTATIONAL BIOLOGY & SCIENCE PIPELINE v3.4 CLINICAL ONCOLOGY")
    print("Multi-Omics Structural, Genomic, Kaplan-Meier Survival & ctDNA MRD Tracking")
    print("=" * 70)

    pipeline_summary = {}

    # Initialize / connect to BioKnowledgebase
    kb = BioKnowledgebase(f"{data_dir}/bio_knowledgebase.db")

    # 1. Structural Biology
    print("\n[1/16] 3D STRUCTURAL BIOLOGY (Hemoglobin 4HHB)")
    struct = StructureAnalyzer(f"{data_dir}/pdb_00004hhb.cif")
    print(f"  Loaded file: {struct.file_path}")
    print(f"  Detected ligands/co-factors: {struct.ligands[:10]}...")
    pipeline_summary['structural_biology'] = {'file': struct.file_path, 'ligand_count': len(struct.ligands)}

    # 2. Genomic Variant Analysis
    print("\n[2/16] GENOMIC VARIANT ANALYSIS (HBB rs334 & ClinVar Pathogenicity)")
    var_analyzer = VariantAnalyzer(data_dir)
    res = var_analyzer.analyze_pathogenicity()
    print(f"  Total variants screened: {res['total_variants']}")
    print(f"  Pathogenic: {res['pathogenic_count']} | Benign: {res['benign_count']} | VUS/Other: {res['vus_other_count']}")
    if res['pathogenic_samples']:
        top_p = res['pathogenic_samples'][0]
        print(f"  Key Pathogenic Mutation: {top_p['title']} [{top_p['clinical_significance']}]")
    pipeline_summary['genomic_variants'] = res

    # 3. Pharmacology & Cheminformatics
    print("\n[3/16] PHARMACOLOGY & DRUG DISCOVERY (Voxelotor & Metformin)")
    pharm = PharmacologyAnalyzer(data_dir)
    vox = pharm.get_compound_profile("voxelotor")
    props = vox.get("properties", {})
    print(f"  Drug: {vox['name']} (Formula: {props.get('MolecularFormula')}, MW: {props.get('MolecularWeight')} g/mol)")
    print(f"  Target Mechanism: {vox.get('mechanism')[:110]}...")
    met = pharm.get_chembl_molecule_summary()
    if met:
        print(f"  ChEMBL Molecule: {met.get('name')} ({met.get('chembl_id')}) - Type: {met.get('type')}")
    pipeline_summary['pharmacology'] = {'voxelotor': props, 'chembl': met}

    # 4. Biological Pathways
    print("\n[4/16] BIOLOGICAL PATHWAYS (Reactome HBB Mapping)")
    pm = PathwayMapper(data_dir)
    pathways = pm.get_gene_pathways("hbb")
    print(f"  Mapped Pathways: {len(pathways)} cellular events identified")
    for p in pathways[:3]:
        print(f"   * [{p['stId']}] {p['name']}")
    pipeline_summary['pathways_count'] = len(pathways)

    # 5. Scientific Literature Intelligence
    print("\n[5/16] SCIENTIFIC LITERATURE & CRISPR BREAKTHROUGH (PubMed)")
    lm = LiteratureMiner(data_dir)
    lit = lm.get_abstract_summary()
    print(f"  Landmark Study (PMID: {lit.get('pmid')}): {lit.get('title')}")
    print(f"  Journal: {lit.get('journal')}")
    pipeline_summary['literature'] = {'pmid': lit.get('pmid'), 'title': lit.get('title')}

    # 6. Pairwise Sequence Alignment & Mutation Consequence
    print("\n[6/16] SEQUENCE ALIGNMENT & BLOSUM62 CONSERVATION")
    aligner = SequenceAligner()
    chosen_variant = align_variant if align_variant else 'HBB_SICKLE'
    if chosen_variant in KNOWN_VARIANTS:
        al_res = aligner.analyze_known_variant(chosen_variant)
        print(f"  Variant Target: {chosen_variant} ({al_res['metadata']['condition']})")
        print(f"  Identity: {al_res['identity_percent']}% | Score: {al_res['alignment_score']}")
        if al_res['mutations_detected']:
            mut = al_res['mutations_detected'][0]
            print(f"  Detected Mutation: {mut['description']} ({mut['type']}, BLOSUM62: {mut['blosum62_score']})")
        pipeline_summary['sequence_alignment'] = al_res

    # 7. In Silico Drug Screening & Lipinski Scoring
    print("\n[7/16] IN SILICO DRUG SCREENING & LIPINSKI PROFILING")
    docking = DockingScorer()
    screen_target = screen_drugs_target if screen_drugs_target else (pan_cancer if pan_cancer else 'EGFR')
    ranked_drugs = docking.screen_library_against_target(screen_target)
    print(f"  Target Pocket: {ranked_drugs[0]['target_pocket']}")
    print(f"  Top In Silico Hit: #{ranked_drugs[0]['screen_rank']} {ranked_drugs[0]['drug_name']}")
    print(f"  Estimated dG: {ranked_drugs[0]['estimated_delta_g_kcal_mol']} kcal/mol (Fitness: {ranked_drugs[0]['contact_fitness_score']}%)")
    pipeline_summary['drug_screening'] = ranked_drugs[:5]

    # 8. Drug Combination Synergy Prediction
    print("\n[8/16] COMBINATION THERAPY & DRUG SYNERGY PREDICTION")
    synergy = SynergyPredictor()
    d_a = synergy_pair[0] if synergy_pair and len(synergy_pair) >= 1 else "Dabrafenib"
    d_b = synergy_pair[1] if synergy_pair and len(synergy_pair) >= 2 else "Trametinib"
    syn_res = synergy.evaluate_combination(d_a, d_b)
    print(f"  Evaluated Regimen: {syn_res['drug_a']} + {syn_res['drug_b']}")
    print(f"  Combination Index (CI): {syn_res['combination_index']} -> {syn_res['synergy_level']}")
    print(f"  Bliss Excess Score: +{syn_res['bliss_score']} | Class: {syn_res['mechanistic_class']}")
    pipeline_summary['drug_synergy'] = syn_res

    # 9. Immunotherapy Efficacy & TMB/MSI
    print("\n[9/16] IMMUNOTHERAPY RESPONSE PREDICTION (TMB & MSI-H)")
    immuno = ImmunotherapyPredictor()
    cur_tmb = tmb_val if tmb_val is not None else 14.5
    imm_res = immuno.evaluate_response(cur_tmb, msi_status, pdl1_val)
    print(f"  Biomarkers: TMB={imm_res['tmb_mut_per_mb']} mut/Mb | MSI={imm_res['msi_status']} | PD-L1={imm_res['pdl1_tps_percent']}%")
    print(f"  Prediction: {imm_res['response_probability_percent']}% Efficacy -> {imm_res['classification']}")
    print(f"  Actionable Regimen: {imm_res['recommended_regimen']}")
    pipeline_summary['immunotherapy'] = imm_res

    # 10. AlphaFold Structure Confidence & IDR Prediction
    print("\n[10/16] ALPHAFOLD STRUCTURAL CONFIDENCE (pLDDT Analysis)")
    af_analyzer = AlphaFoldQualityAnalyzer()
    target_gene = alphafold_gene if alphafold_gene else "TP53"
    af_res = af_analyzer.analyze_protein(target_gene)
    print(f"  Target: {af_res['gene']} ({af_res['protein_name']}) - UniProt: {af_res['uniprot_id']}")
    print(f"  Mean pLDDT: {af_res['mean_plddt']} ({af_res['overall_band']}) | Disordered Fraction: {af_res['disorder_fraction_percent']}%")
    pipeline_summary['alphafold_quality'] = af_res

    # 11. Relational Database Engine
    print("\n[11/16] UNIFIED RELATIONAL BIO-KNOWLEDGEBASE (SQLite 3)")
    db_stats = kb.get_database_stats()
    print(f"  Database File: {kb.db_path}")
    print(f"  Entities: {db_stats['genes_count']} Genes, {db_stats['variants_count']} Variants, {db_stats['structures_count']} Structures, {db_stats['compounds_count']} Compounds, {db_stats['pathways_count']} Pathways")
    print(f"  Clinical Benchmarks: {db_stats['clinical_survival_benchmarks_count']} Trials | ctDNA Clones: {db_stats['liquid_biopsy_clones_count']}")
    pipeline_summary['database_stats'] = db_stats

    # 12. Early Cancer Detection & 100k Population Epidemiology Engine
    print("\n[12/16] EARLY CANCER DETECTION & POPULATION EPIDEMIOLOGY (100k Cohort)")
    risk_pred = CancerRiskPredictor()
    cohort_count = risk_pred.get_cohort_count()
    print(f"  Cohort Database: {risk_pred.db_path} ({cohort_count:,} Records)")
    pipeline_summary['cohort_database'] = {'path': risk_pred.db_path, 'records': cohort_count}

    # 13. Pharmacogenomics (PGx) Safety Engine
    print("\n[13/16] PHARMACOGENOMICS (PGx) CLINICAL TOXICITY (CPIC Level A)")
    dpyd_a1 = dpyd_alleles[0] if dpyd_alleles else "*1"
    dpyd_a2 = dpyd_alleles[1] if dpyd_alleles else "*1"
    pgx_profile = PharmacogenomicsEngine.profile_all_oncology_pgx(
        dpyd_diplotype=(dpyd_a1, dpyd_a2),
        ugt1a1_diplotype=ugt1a1_diplotype if ugt1a1_diplotype else "*1/*1",
        cyp2d6_diplotype=cyp2d6_diplotype if cyp2d6_diplotype else "*1/*1"
    )
    dpyd_eval = pgx_profile['dpyd_fluoropyrimidine']
    print(f"  DPYD Diplotype: {dpyd_eval['diplotype']} (Activity Score: {dpyd_eval['activity_score']}) -> {dpyd_eval['phenotype']}")
    print(f"  Fluoropyrimidine Toxicity: {dpyd_eval['risk_level']} | Guideline: {dpyd_eval['clinical_action']}")
    pipeline_summary['pharmacogenomics'] = pgx_profile

    # 14. Synthetic Lethality & DepMap Essentiality
    print("\n[14/16] SYNTHETIC LETHALITY & CRISPR DEPMAP DEPENDENCY")
    sl_query = synthetic_lethal_target if synthetic_lethal_target else "BRCA1"
    sl_res = SyntheticLethalityEngine.query_vulnerabilities(sl_query)
    print(f"  Biomarker: {sl_res.get('biomarker_gene', sl_query)} -> Target: {sl_res.get('synthetic_lethal_target', 'N/A')}")
    print(f"  DepMap Chronos Score: {sl_res.get('depmap_chronos_score', 0.0)} | Vulnerability: %{sl_res.get('vulnerability_score', 0.0)}")
    if sl_res.get('clinical_therapies'):
        print(f"  Primary Targeted Drug: {sl_res['clinical_therapies'][0]['drug']} ({sl_res['clinical_therapies'][0]['fda_status']})")
    pipeline_summary['synthetic_lethality'] = sl_res

    # 15. Clinical Survival Analysis & Log-Rank Modeler
    print("\n[15/16] CLINICAL SURVIVAL ANALYSIS & LOG-RANK HYPOTHESIS MODELER (Kaplan-Meier)")
    surv_modeler = SurvivalAnalysisModeler()
    chosen_trial = survival_trial if survival_trial else "FLAURA_EGFR"
    surv_res = surv_modeler.evaluate_clinical_trial_cohort(chosen_trial)
    if "error" in surv_res:
        print(f"  Trial: {chosen_trial} (Custom / Not in Preset Benchmarks)")
    else:
        print(f"  Benchmark Trial: {surv_res['trial_name']} ({surv_res['target_gene']})")
        print(f"  Indication: {surv_res['indication']} | Endpoint: {surv_res['endpoint']}")
        print(f"  Experimental Arm: {surv_res['experimental_arm']} (Median: {surv_res['experimental_median_months']} mo)")
        print(f"  Control Arm: {surv_res['control_arm']} (Median: {surv_res['control_median_months']} mo)")
        print(f"  Hazard Ratio (HR): {surv_res['hazard_ratio']} (95% CI: [{surv_res['hr_ci_95'][0]}, {surv_res['hr_ci_95'][1]}])")
        print(f"  Log-Rank p-value: {surv_res['p_value']} | Significant: {surv_res['statistically_significant']}")
        print(f"  Clinical Conclusion: {surv_res['clinical_conclusion']}")
    pipeline_summary['survival_analysis'] = surv_res

    # 16. Liquid Biopsy & ctDNA Minimal Residual Disease
    print("\n[16/16] LIQUID BIOPSY & CIRCULATING TUMOR DNA (ctDNA) MRD DYNAMICS")
    lb_engine = LiquidBiopsyEngine()
    ctdna_muts = ctdna_mutations if ctdna_mutations else ["EGFR_C797S"]
    lb_res = lb_engine.evaluate_longitudinal_ctdna(
        baseline_vaf_percent=ctdna_baseline_vaf if ctdna_baseline_vaf is not None else 4.8,
        followup_vaf_percent=ctdna_followup_vaf if ctdna_followup_vaf is not None else 0.02,
        primary_gene="EGFR",
        detected_new_mutations=ctdna_muts
    )
    print(f"  Longitudinal VAF: {lb_res['baseline_vaf_percent']}% -> {lb_res['followup_vaf_percent']}% (Delta: {lb_res['delta_vaf_percent']}%)")
    print(f"  Molecular Response: {lb_res['response_classification']} [{lb_res['response_code']}]")
    print(f"  MRD Status: {lb_res['mrd_status']} | Relapse Risk: {lb_res['relapse_risk']}")
    print(f"  Clinical Recommendation: {lb_res['clinical_action']}")
    if lb_res['resistance_findings']:
        print(f"  Acquired Resistance Clone: {lb_res['resistance_findings'][0]['gene']} {lb_res['resistance_findings'][0]['mutation']}")
        print(f"  Actionable Bypass Therapy: {lb_res['resistance_findings'][0]['actionable_bypass']}")
    pipeline_summary['liquid_biopsy'] = lb_res

    # OPTIONAL: Build / Rebuild Cohort
    if build_cohort_size:
        print("\n" + "-" * 70)
        print(f"GENERATING EPIDEMIOLOGICAL COHORT ({build_cohort_size:,} Records)")
        print("-" * 70)
        risk_pred.generate_and_seed_cohort(total_records=build_cohort_size)
        print(f"  Cohort successfully generated with {build_cohort_size:,} population records!")

    # OPTIONAL: Cancer Risk Evaluation
    if predict_cancer_risk or patient_symptoms:
        print("\n" + "-" * 70)
        print("INDIVIDUAL EARLY CANCER RISK EVALUATION (NICE NG12 Guidelines)")
        print("-" * 70)
        symptoms_list = [s.strip() for s in patient_symptoms.split(",")] if patient_symptoms else []
        eval_result = risk_pred.evaluate_patient(
            age=patient_age,
            gender=patient_gender,
            symptoms=symptoms_list,
            smoking_pack_years=smoking_pack_years,
            family_history=family_history
        )
        p_info = eval_result['patient_profile']
        r_eval = eval_result['risk_evaluation']
        c_evid = eval_result['cohort_empirical_evidence']

        print(f"  Patient Profile: {p_info['age']} y/o {p_info['gender']} | Smoking: {p_info['smoking_pack_years']} pack-years | Family Hx: {p_info['family_history']}")
        print(f"  Active Symptoms: {', '.join(p_info['active_symptoms']) if p_info['active_symptoms'] else 'None (Asymptomatic)'}")
        print(f"\n  [PREDICTION RESULT]")
        print(f"  Composite Cancer Risk: {r_eval['composite_risk_percent']}% -> {r_eval['risk_tier']}")
        print(f"  NICE 2-Week-Wait (2WW) Urgent Referral: {'REQUIRED (YES)' if r_eval['meets_nice_2ww_referral'] else 'Not required'}")
        if r_eval['urgent_reasons']:
            print(f"    * Referral Triggers: {', '.join(r_eval['urgent_reasons'])}")
        print(f"  Clinical Action: {r_eval['clinical_action']}")
        print(f"  Primary Recommended Test: {r_eval['primary_diagnostic_recommendation']}")
        print(f"\n  [100,000-PATIENT COHORT EVIDENCE]")
        print(f"  Matched Cohort Subgroup: {c_evid['matched_cohort_patients']:,} patients (Age {c_evid['age_bracket']})")
        print(f"  Observed Cohort Cancer Rate: %{c_evid['cohort_incidence_percent']} ({c_evid['cohort_diagnosed_cancers']} cases)")
        print(f"\n  Top Suspected Cancer Types:")
        for org in eval_result['organ_probabilities']:
            print(f"    -> {org['cancer_name']}: %{org['probability_percent']} | Recommended: {org['recommended_test']}")

        pipeline_summary['cancer_risk_evaluation'] = eval_result

    # OPTIONAL: Database Query
    if db_query_gene or show_db_stats:
        print("\n" + "-" * 70)
        print("RELATIONAL DATABASE QUERY CONSOLE")
        print("-" * 70)
        if show_db_stats:
            print(f"  Knowledgebase Stats: {json.dumps(db_stats, indent=2)}")
        if db_query_gene:
            print(f"\n  [QUERY GENE: {db_query_gene.upper()}]")
            g_vars = kb.query_variants(gene=db_query_gene)
            print(f"  Found {len(g_vars)} Clinical Variants:")
            for v in g_vars:
                print(f"    * {v['hgvs_p']} ({v['rs_id']}) - {v['clinical_significance']} [{v['disease_name']}]")
            g_structs = kb.query_structures(gene=db_query_gene)
            print(f"  Found {len(g_structs)} 3D Structures:")
            for s in g_structs:
                print(f"    * PDB: {s['pdb_id']} ({s['resolution_angstrom']} A) - {s['title'][:55]}...")
            g_drugs = kb.query_compounds(target=db_query_gene)
            print(f"  Found {len(g_drugs)} Targeted Compounds:")
            for d in g_drugs:
                print(f"    * {d['common_name']} ({d['molecular_formula']}) - MW: {d['molecular_weight']} | FDA: {d['fda_status']}")

    # OPTIONAL: Export Database to JSON
    if export_db_json_path:
        out_f = kb.export_full_json(export_db_json_path)
        print(f"\n  Knowledgebase exported to JSON: {out_f}")

    # OPTIONAL: PPI Network Analysis
    if ppi_gene:
        print("\n" + "-" * 70)
        print("PROTEIN-PROTEIN INTERACTION (PPI) NETWORK TOPOLOGY")
        print("-" * 70)
        ppi = InteractionNetwork(ppi_gene)
        ppi_data = ppi.compute_network_metrics()
        print(f"  Network: {ppi_data['network_name']}")
        print(f"  Nodes: {ppi_data['total_nodes']} | Edges: {ppi_data['total_edges']} | Density: {ppi_data['network_density']}")
        print(f"  Primary Bottleneck Hub: {ppi_data['primary_hub_protein']} | Centrality Hubs: {', '.join(ppi_data['top_ranked_hubs'])}")
        pipeline_summary['ppi_network'] = ppi_data

    # OPTIONAL: Dynamic Live Fetcher Extension
    if live_pdb or live_drug or live_pubmed:
        print("\n" + "-" * 70)
        print("LIVE SCIENTIFIC API QUERIES (Dynamic Search Engine)")
        print("-" * 70)
        fetcher = LiveFetcher(f"{data_dir}/cache")

        if live_pdb:
            print(f"\n[LIVE PDB] Querying: {live_pdb}")
            pdb_info = fetcher.fetch_pdb_entry(live_pdb)
            from_cache = " (Cached)" if pdb_info.get("_from_cache") else " (Live API)"
            print(f"  Title: {pdb_info.get('title')}{from_cache}")
            print(f"  Resolution: {pdb_info.get('resolution')} A | Method: {pdb_info.get('method')}")

        if live_drug:
            print(f"\n[LIVE PUBCHEM] Querying: {live_drug}")
            drug_info = fetcher.fetch_pubchem_drug(live_drug)
            from_cache = " (Cached)" if drug_info.get("_from_cache") else " (Live API)"
            print(f"  CID: {drug_info.get('cid')} | Formula: {drug_info.get('formula')}{from_cache}")
            print(f"  Molecular Weight: {drug_info.get('molecular_weight')} g/mol | XLogP: {drug_info.get('xlogp')}")

        if live_pubmed:
            print(f"\n[LIVE PUBMED] Searching: '{live_pubmed}'")
            pm_info = fetcher.fetch_pubmed_articles(live_pubmed, max_results=2)
            articles = pm_info.get("articles", [])
            for a in articles:
                print(f"   * [PMID: {a.get('pmid')}] {a.get('title')} ({a.get('source')})")

    # OPTIONAL: Pan-Cancer Oncology Profiler
    if pan_cancer:
        print("\n" + "-" * 70)
        print("PAN-CANCER ONCOLOGY & TARGETED THERAPEUTICS")
        print("-" * 70)
        profiler = OncologyProfiler()
        profile = profiler.profile_gene(pan_cancer)
        if "error" in profile:
            print(f"  Error: {profile['error']}")
        else:
            print(f"  Target Gene: {profile['gene']} ({profile['full_name']})")
            print(f"  Associated Malignancies: {', '.join(profile['cancer_types'])}")
            print(f"  Cataloged Hotspot Mutations ({profile['hotspot_count']}):")
            for h in profile['hotspots']:
                print(f"    * {h['variant']} [{h.get('type')}]: {h.get('response')}")
            print("  FDA-Approved Targeted Therapies:")
            for d in profile['targeted_therapies']:
                print(f"    -> {d['drug']} ({d.get('line', d.get('status'))})")
        pipeline_summary['pan_cancer_profile'] = profile

    # OPTIONAL: Export Pipeline Run JSON
    if export_json_path:
        with open(export_json_path, 'w', encoding='utf-8') as f:
            json.dump(pipeline_summary, f, indent=2, ensure_ascii=False)
        print(f"\n  Pipeline Run JSON saved to: {export_json_path}")

    # OPTIONAL: Publication HTML Report Generation
    if generate_report:
        print("\n" + "-" * 70)
        print("PUBLICATION REPORT GENERATION")
        print("-" * 70)
        rg = ReportGenerator(data_dir=data_dir, output_dir="reports")
        rep_file = rg.generate_html_report(target_gene=screen_target)
        print(f"  Executive Report Generated: {rep_file}")

    print("\n" + "=" * 70)
    print("PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
    print("=" * 70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Computational Biology Pipeline v3.1 Database Edition")
    parser.add_argument("--data_dir", default="data", help="Path to data directory")
    parser.add_argument("--live_pdb", default=None, help="Fetch live PDB metadata (e.g. 1TSR, 6LU7)")
    parser.add_argument("--live_drug", default=None, help="Fetch live PubChem drug properties (e.g. Imatinib, Aspirin)")
    parser.add_argument("--live_pubmed", default=None, help="Search live PubMed literature (e.g. 'CAR-T lymphoma')")
    parser.add_argument("--pan_cancer", default=None, help="Profile oncogene (e.g. EGFR, KRAS, BRAF, TP53)")
    parser.add_argument("--align", default=None, help="Pairwise sequence alignment (e.g. HBB_SICKLE, EGFR_L858R, KRAS_G12C)")
    parser.add_argument("--screen_drugs", default=None, help="Screen candidate drugs against target pocket (e.g. EGFR, KRAS, BRAF, HBB)")
    parser.add_argument("--ppi", default=None, help="Analyze PPI interactome network (e.g. EGFR, KRAS, TP53, HBB)")
    parser.add_argument("--synergy", nargs=2, default=None, help="Evaluate combination therapy synergy (e.g. --synergy Dabrafenib Trametinib)")
    parser.add_argument("--tmb", type=float, default=None, help="Tumor Mutational Burden (mut/Mb, e.g. 15.5)")
    parser.add_argument("--msi", default="MSS", help="Microsatellite Instability status (MSS or MSI-H)")
    parser.add_argument("--pdl1", type=float, default=0.0, help="PD-L1 Tumor Proportion Score (TPS %)")
    parser.add_argument("--alphafold_conf", default=None, help="Analyze AlphaFold pLDDT confidence (e.g. TP53, EGFR, HBB)")
    parser.add_argument("--db_query", default=None, help="Query relational knowledgebase for a gene (e.g. EGFR, KRAS, HBB)")
    parser.add_argument("--db_stats", action="store_true", help="Display global relational database statistics")
    parser.add_argument("--export_db", default=None, help="Export entire SQLite knowledgebase to JSON file")
    parser.add_argument("--predict_cancer", action="store_true", help="Run early cancer detection statistical risk engine")
    parser.add_argument("--patient_age", type=int, default=50, help="Patient age in years (e.g. 58)")
    parser.add_argument("--patient_gender", default="M", help="Patient biological sex (M or F)")
    parser.add_argument("--patient_symptoms", default=None, help="Comma-separated clinical symptoms (e.g. 'hemoptysis,weight_loss,chronic_cough')")
    parser.add_argument("--smoking_pack_years", type=float, default=0.0, help="Smoking history in pack-years (e.g. 25.0)")
    parser.add_argument("--family_history", action="store_true", help="First-degree family history of malignant neoplasm")
    parser.add_argument("--build_cohort", type=int, default=None, help="Generate / seed epidemiological population cohort size (e.g. 100000)")
    parser.add_argument("--pgx_check", action="store_true", help="Run pharmacogenomics toxicity and CPIC dosing analysis")
    parser.add_argument("--dpyd_alleles", nargs=2, default=None, help="DPYD alleles, e.g. --dpyd_alleles *1 *2A")
    parser.add_argument("--ugt1a1", default=None, help="UGT1A1 diplotype, e.g. --ugt1a1 '*28/*28'")
    parser.add_argument("--cyp2d6", default=None, help="CYP2D6 diplotype, e.g. --cyp2d6 '*4/*4'")
    parser.add_argument("--synthetic_lethal", default=None, help="Query synthetic lethality target dependency (e.g. BRCA1, MSI-H, MTAP, ARID1A, KRAS)")
    parser.add_argument("--survival_trial", default="FLAURA_EGFR", help="Benchmark trial for survival analysis (FLAURA_EGFR, KEYNOTE_042_TMB_HIGH, OLYMPIAD_BRCA, CODEBREAK_200_KRAS)")
    parser.add_argument("--ctdna_baseline", type=float, default=4.8, help="Baseline ctDNA Variant Allele Frequency (VAF %)")
    parser.add_argument("--ctdna_followup", type=float, default=0.02, help="Follow-up ctDNA Variant Allele Frequency (VAF %)")
    parser.add_argument("--ctdna_mutations", nargs="+", default=None, help="Emerging resistance mutations in liquid biopsy, e.g. EGFR_C797S MET_AMPLIFICATION")
    parser.add_argument("--export_json", default=None, help="Export complete pipeline results to JSON file")
    parser.add_argument("--generate_report", action="store_true", help="Generate publication-grade HTML report")
    args = parser.parse_args()
    
    run_pipeline(
        data_dir=args.data_dir,
        live_pdb=args.live_pdb,
        live_drug=args.live_drug,
        live_pubmed=args.live_pubmed,
        pan_cancer=args.pan_cancer,
        align_variant=args.align,
        screen_drugs_target=args.screen_drugs,
        ppi_gene=args.ppi,
        synergy_pair=args.synergy,
        tmb_val=args.tmb,
        msi_status=args.msi,
        pdl1_val=args.pdl1,
        alphafold_gene=args.alphafold_conf,
        db_query_gene=args.db_query,
        show_db_stats=args.db_stats,
        predict_cancer_risk=args.predict_cancer,
        patient_age=args.patient_age,
        patient_gender=args.patient_gender,
        patient_symptoms=args.patient_symptoms,
        smoking_pack_years=args.smoking_pack_years,
        family_history=args.family_history,
        build_cohort_size=args.build_cohort,
        pgx_check=args.pgx_check,
        dpyd_alleles=args.dpyd_alleles,
        ugt1a1_diplotype=args.ugt1a1,
        cyp2d6_diplotype=args.cyp2d6,
        synthetic_lethal_target=args.synthetic_lethal,
        survival_trial=args.survival_trial,
        ctdna_baseline_vaf=args.ctdna_baseline,
        ctdna_followup_vaf=args.ctdna_followup,
        ctdna_mutations=args.ctdna_mutations,
        export_db_json_path=args.export_db,
        export_json_path=args.export_json,
        generate_report=args.generate_report
    )
