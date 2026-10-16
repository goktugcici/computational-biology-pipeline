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

if __name__ == "__main__":
    unittest.main()
