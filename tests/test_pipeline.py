"""
Unit Tests for Computational Biology Pipeline (v3.1 Database Edition)
Tests structural, variant, pharmacology, pathway, literature,
sequence alignment, docking scorer, interaction network,
drug combination synergy, immunotherapy response, AlphaFold quality,
and relational database query capabilities.
"""

import os
import unittest
import sys

# src dizinini path'e ekle
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from structural_biology import StructureAnalyzer
from variant_analyzer import VariantAnalyzer
from pharmacology import PharmacologyAnalyzer
from pathway_mapper import PathwayMapper
from literature_miner import LiteratureMiner
from live_fetcher import LiveFetcher
from oncology_profiler import OncologyProfiler
from sequence_aligner import SequenceAligner
from docking_scorer import DockingScorer
from interaction_network import InteractionNetwork
from synergy_predictor import SynergyPredictor
from immunotherapy_response import ImmunotherapyPredictor
from alphafold_quality import AlphaFoldQualityAnalyzer
from database import BioKnowledgebase
from report_generator import ReportGenerator

class TestComputationalBiologyPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))
        cls.output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "reports"))
        cls.db_path = os.path.join(cls.data_dir, "bio_knowledgebase.db")

    @classmethod
    def tearDownClass(cls):
        test_file = os.path.join(cls.output_dir, "test_summary.html")
        if os.path.exists(test_file):
            try:
                os.remove(test_file)
            except OSError:
                pass

    def test_structural_analyzer(self):
        cif_path = os.path.join(self.data_dir, "pdb_00004hhb.cif")
        analyzer = StructureAnalyzer(cif_path)
        summary = analyzer.get_summary()
        self.assertIn("file", summary)
        self.assertEqual(summary["file"], "pdb_00004hhb.cif")
        self.assertGreater(summary["size_bytes"], 0)

    def test_variant_analyzer(self):
        analyzer = VariantAnalyzer(self.data_dir)
        res = analyzer.analyze_pathogenicity()
        self.assertGreaterEqual(res["total_variants"], 1)
        self.assertGreaterEqual(res["pathogenic_count"], 1)
        self.assertIsInstance(res["pathogenic_samples"], list)

    def test_pharmacology_analyzer(self):
        analyzer = PharmacologyAnalyzer(self.data_dir)
        vox = analyzer.get_compound_profile("voxelotor")
        self.assertEqual(vox["name"], "Voxelotor")
        props = vox.get("properties", {})
        self.assertIn("MolecularFormula", props)
        self.assertEqual(props.get("MolecularFormula"), "C19H19N3O3")

        # ChEMBL Metformin check
        chembl = analyzer.get_chembl_molecule_summary()
        self.assertEqual(chembl.get("chembl_id"), "CHEMBL1431")
        self.assertEqual(chembl.get("name"), "METFORMIN")

    def test_pathway_mapper(self):
        mapper = PathwayMapper(self.data_dir)
        pathways = mapper.get_gene_pathways("hbb")
        self.assertIsInstance(pathways, list)
        self.assertGreater(len(pathways), 0)
        st_ids = [p["stId"] for p in pathways]
        self.assertTrue(any("R-HSA" in sid for sid in st_ids))

    def test_literature_miner(self):
        miner = LiteratureMiner(self.data_dir)
        lit = miner.get_abstract_summary()
        self.assertEqual(lit.get("pmid"), "38661449")
        self.assertIn("Exagamglogene", lit.get("title", ""))

    def test_live_fetcher_cached(self):
        fetcher = LiveFetcher(os.path.join(self.data_dir, "cache"))
        drug = fetcher.fetch_pubchem_drug("Imatinib")
        self.assertIn("name", drug)
        self.assertEqual(drug["name"], "Imatinib")
        self.assertEqual(drug.get("cid"), 5291)

    def test_oncology_profiler(self):
        profiler = OncologyProfiler()
        egfr = profiler.profile_gene("EGFR")
        self.assertEqual(egfr["gene"], "EGFR")
        self.assertGreater(egfr["hotspot_count"], 0)
        self.assertIn("Osimertinib", egfr["targeted_therapies"][0]["drug"])

        pan = profiler.pan_cancer_summary()
        self.assertIn("EGFR", pan["genes_covered"])
        self.assertIn("KRAS", pan["genes_covered"])

    def test_sequence_aligner(self):
        aligner = SequenceAligner()
        res = aligner.analyze_known_variant("HBB_SICKLE")
        self.assertGreater(res["alignment_score"], 0)
        self.assertGreater(res["identity_percent"], 95.0)
        self.assertEqual(len(res["mutations_detected"]), 1)
        self.assertEqual(res["mutations_detected"][0]["ref"], "E")
        self.assertEqual(res["mutations_detected"][0]["alt"], "V")

    def test_docking_scorer(self):
        docking = DockingScorer()
        ranked = docking.screen_library_against_target("EGFR")
        self.assertGreaterEqual(len(ranked), 5)
        self.assertEqual(ranked[0]["drug_id"], "OSIMERTINIB")
        self.assertLess(ranked[0]["estimated_delta_g_kcal_mol"], -10.0)
        self.assertTrue(ranked[0]["lipinski_pass"])

    def test_interaction_network(self):
        ppi = InteractionNetwork("EGFR")
        metrics = ppi.compute_network_metrics()
        self.assertEqual(metrics["gene"], "EGFR")
        self.assertGreaterEqual(metrics["total_nodes"], 5)
        self.assertGreaterEqual(metrics["total_edges"], 5)
        self.assertEqual(metrics["primary_hub_protein"], "EGFR")
        self.assertIn("EGFR", metrics["node_centrality"])

    def test_synergy_predictor(self):
        synergy = SynergyPredictor()
        combo = synergy.evaluate_combination("Dabrafenib", "Trametinib")
        self.assertLess(combo["combination_index"], 0.8) # Synergistic
        self.assertGreater(combo["bliss_score"], 10.0)
        self.assertTrue(combo["is_curated_trial"])

    def test_immunotherapy_response(self):
        immuno = ImmunotherapyPredictor()
        res_high = immuno.evaluate_response(tmb_mut_per_mb=22.0, msi_status="MSI-H", pdl1_tps_percent=75.0)
        self.assertTrue(res_high["is_tmb_high"])
        self.assertTrue(res_high["is_msi_high"])
        self.assertGreater(res_high["response_probability_percent"], 80.0)

        res_low = immuno.evaluate_response(tmb_mut_per_mb=2.0, msi_status="MSS", pdl1_tps_percent=0.0)
        self.assertFalse(res_low["is_tmb_high"])
        self.assertLess(res_low["response_probability_percent"], 35.0)

    def test_alphafold_quality(self):
        af = AlphaFoldQualityAnalyzer()
        res = af.analyze_protein("TP53")
        self.assertEqual(res["gene"], "TP53")
        self.assertGreater(res["mean_plddt"], 50.0)
        self.assertGreater(len(res["structured_domains"]), 0)
        self.assertGreater(res["disorder_fraction_percent"], 20.0)

    def test_database_init_and_stats(self):
        kb = BioKnowledgebase(self.db_path)
        stats = kb.get_database_stats()
        self.assertGreaterEqual(stats["genes_count"], 5)
        self.assertGreaterEqual(stats["variants_count"], 10)
        self.assertGreaterEqual(stats["structures_count"], 5)
        self.assertGreaterEqual(stats["compounds_count"], 5)
        self.assertIn("Pathogenic", stats["variants_by_significance"])

    def test_database_queries(self):
        kb = BioKnowledgebase(self.db_path)
        # Query gene
        egfr_vars = kb.query_variants(gene="EGFR")
        self.assertGreaterEqual(len(egfr_vars), 1)
        self.assertEqual(egfr_vars[0]["gene_symbol"], "EGFR")

        # Query drug
        drugs = kb.query_compounds(target="KRAS")
        self.assertGreaterEqual(len(drugs), 1)
        self.assertIn("Sotorasib", drugs[0]["common_name"])

        # Query structure
        structs = kb.query_structures(gene="HBB")
        self.assertGreaterEqual(len(structs), 1)
        self.assertEqual(structs[0]["pdb_id"], "4HHB")

    def test_report_generator(self):
        rg = ReportGenerator(data_dir=self.data_dir, output_dir=self.output_dir)
        report_path = rg.generate_html_report("test_summary.html")
        self.assertTrue(os.path.exists(report_path))
        self.assertGreater(os.path.getsize(report_path), 500)

    def test_cancer_cohort_initialization(self):
        from cancer_risk_predictor import CancerRiskPredictor
        predictor = CancerRiskPredictor()
        self.assertGreaterEqual(predictor.get_cohort_count(), 100000)

    def test_cancer_risk_high_risk_smoker(self):
        from cancer_risk_predictor import CancerRiskPredictor
        predictor = CancerRiskPredictor()
        eval_res = predictor.evaluate_patient(
            age=65,
            gender='M',
            symptoms=['hemoptysis', 'chronic_cough', 'unexplained_weight_loss'],
            smoking_pack_years=35.0,
            family_history=True
        )
        self.assertGreaterEqual(eval_res['risk_evaluation']['composite_risk_percent'], 75.0)
        self.assertTrue(eval_res['risk_evaluation']['meets_nice_2ww_referral'])
        self.assertEqual(eval_res['risk_evaluation']['tier_code'], 'CRITICAL')
        self.assertIn("Akciğer", eval_res['organ_probabilities'][0]['cancer_name'])

    def test_cancer_risk_asymptomatic_low_risk(self):
        from cancer_risk_predictor import CancerRiskPredictor
        predictor = CancerRiskPredictor()
        eval_res = predictor.evaluate_patient(
            age=25,
            gender='F',
            symptoms=[],
            smoking_pack_years=0.0
        )
        self.assertLess(eval_res['risk_evaluation']['composite_risk_percent'], 3.0)
        self.assertFalse(eval_res['risk_evaluation']['meets_nice_2ww_referral'])
        self.assertEqual(eval_res['risk_evaluation']['tier_code'], 'LOW')

    def test_cancer_risk_organ_specificity(self):
        from cancer_risk_predictor import CancerRiskPredictor
        predictor = CancerRiskPredictor()

        # Colorectal test
        col_res = predictor.evaluate_patient(
            age=60,
            gender='M',
            symptoms=['rectal_bleeding', 'iron_deficiency_anemia']
        )
        self.assertIn("Kolorektal", col_res['organ_probabilities'][0]['cancer_name'])

        # Breast test
        br_res = predictor.evaluate_patient(
            age=45,
            gender='F',
            symptoms=['breast_lump', 'nipple_abnormality']
        )
        self.assertIn("Meme", br_res['organ_probabilities'][0]['cancer_name'])

    def test_pharmacogenomics_dpyd_toxicity(self):
        from pharmacogenomics import PharmacogenomicsEngine
        im_res = PharmacogenomicsEngine.evaluate_dpyd('*1', '*2A')
        self.assertEqual(im_res['activity_score'], 1.0)
        self.assertEqual(im_res['phenotype_code'], 'IM')
        self.assertIn("50", im_res['clinical_action'])

        pm_res = PharmacogenomicsEngine.evaluate_dpyd('*2A', '*2A')
        self.assertEqual(pm_res['activity_score'], 0.0)
        self.assertEqual(pm_res['phenotype_code'], 'PM')
        self.assertEqual(pm_res['action_code'], 'CONTRAINDICATED')

    def test_pharmacogenomics_ugt1a1_and_cyp2d6(self):
        from pharmacogenomics import PharmacogenomicsEngine
        ugt_res = PharmacogenomicsEngine.evaluate_ugt1a1('*28/*28')
        self.assertEqual(ugt_res['phenotype_code'], 'PM')
        self.assertIn("azalt", ugt_res['clinical_action'])

        cyp_res = PharmacogenomicsEngine.evaluate_cyp2d6_tamoxifen('*4/*4')
        self.assertEqual(cyp_res['phenotype_code'], 'PM')
        self.assertIn("Aromataz", cyp_res['clinical_action'])

    def test_synthetic_lethality_brca_parp(self):
        from synthetic_lethality import SyntheticLethalityEngine
        sl_res = SyntheticLethalityEngine.query_vulnerabilities('BRCA1')
        self.assertEqual(sl_res['status'], 'synthetic_lethal_identified')
        self.assertEqual(sl_res['synthetic_lethal_target'], 'PARP1')
        self.assertLess(sl_res['depmap_chronos_score'], -1.0)
        self.assertGreaterEqual(sl_res['vulnerability_score'], 90.0)
        self.assertIn("Olaparib", sl_res['clinical_therapies'][0]['drug'])

    def test_synthetic_lethality_msi_wrn(self):
        from synthetic_lethality import SyntheticLethalityEngine
        wrn_res = SyntheticLethalityEngine.query_vulnerabilities('MSI-H')
        self.assertEqual(wrn_res['status'], 'synthetic_lethal_identified')
        self.assertEqual(wrn_res['synthetic_lethal_target'], 'WRN')
        self.assertLess(wrn_res['depmap_chronos_score'], -1.2)

    def test_database_pgx_and_synthetic_lethality(self):
        kb = BioKnowledgebase(self.db_path)
        dpyd_guidelines = kb.query_pgx("DPYD")
        self.assertGreaterEqual(len(dpyd_guidelines), 1)
        self.assertIn("Fluorouracil", dpyd_guidelines[0]['drug_name'])

        sl_pairs = kb.query_synthetic_lethality("BRCA1")
        self.assertGreaterEqual(len(sl_pairs), 1)
        self.assertEqual(sl_pairs[0]['synthetic_lethal_target'], 'PARP1')

    def test_survival_kaplan_meier_product_limit(self):
        from survival_analysis import SurvivalAnalyzer
        durations = [3.0, 6.0, 9.0, 12.0, 15.0, 18.0]
        events = [1, 1, 0, 1, 0, 1]
        km_res = SurvivalAnalyzer.calculate_kaplan_meier(durations, events)
        self.assertEqual(km_res['total_patients'], 6)
        self.assertEqual(km_res['total_events'], 4)
        self.assertIn('curve_points', km_res)
        pts = km_res['curve_points']
        self.assertGreater(len(pts), 0)
        for i in range(len(pts) - 1):
            self.assertGreaterEqual(pts[i]['survival_probability'], pts[i + 1]['survival_probability'])

    def test_survival_log_rank_test_statistics(self):
        from survival_analysis import SurvivalAnalyzer
        exp_d = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
        exp_e = [0, 0, 1, 0, 1, 0]
        ctrl_d = [2.0, 4.0, 6.0, 8.0, 12.0, 16.0]
        ctrl_e = [1, 1, 1, 1, 1, 1]
        logrank = SurvivalAnalyzer.compute_log_rank_test(exp_d, exp_e, ctrl_d, ctrl_e)
        self.assertIn('log_rank_chisquare', logrank)
        self.assertIn('hazard_ratio', logrank)
        self.assertLess(logrank['hazard_ratio'], 1.0)
        self.assertTrue(logrank['statistically_significant'])

    def test_survival_benchmark_cohort_retrieval(self):
        from survival_analysis import SurvivalAnalyzer
        flaura = SurvivalAnalyzer.get_benchmark_cohort('FLAURA_EGFR')
        self.assertEqual(flaura['hazard_ratio'], 0.80)
        self.assertIn("Osimertinib", flaura['experimental_arm']['name'])

        keynote = SurvivalAnalyzer.get_benchmark_cohort('KEYNOTE_042_TMB_HIGH')
        self.assertEqual(keynote['hazard_ratio'], 0.69)
        self.assertIn("Pembrolizumab", keynote['experimental_arm']['name'])

    def test_liquid_biopsy_longitudinal_evaluation(self):
        from liquid_biopsy import LiquidBiopsyEngine
        # Case 1: ctDNA clearance
        res_mcr = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(
            baseline_vaf_percent=5.2,
            followup_vaf_percent=0.01,
            primary_gene="EGFR"
        )
        self.assertEqual(res_mcr['response_code'], 'MCR')
        self.assertIn("MRD Negatif", res_mcr['mrd_status'])

        # Case 2: Resistance mutation C797S
        res_prog = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(
            baseline_vaf_percent=0.5,
            followup_vaf_percent=3.8,
            primary_gene="EGFR",
            detected_new_mutations=["EGFR_C797S"]
        )
        self.assertEqual(res_prog['response_code'], 'MPD')
        self.assertTrue(res_prog['acquired_resistance_detected'])
        self.assertEqual(res_prog['resistance_clones'][0]['mutation_key'], 'EGFR_C797S')

    def test_database_survival_and_liquid_biopsy_tables(self):
        kb = BioKnowledgebase(self.db_path)
        benchmarks = kb.query_survival_benchmarks("FLAURA_EGFR")
        self.assertEqual(len(benchmarks), 1)
        self.assertEqual(benchmarks[0]['target_gene'], "EGFR (Exon 19 del / L858R)")

        clones = kb.query_liquid_biopsy_clones("EGFR")
        self.assertGreaterEqual(len(clones), 1)
        self.assertEqual(clones[0]['mutation_key'], "EGFR_C797S")

if __name__ == "__main__":
    unittest.main()
