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

        # Gynecologic test (Postmenopausal bleeding in female)
        gyn_res = predictor.evaluate_patient(
            age=56,
            gender='F',
            symptoms=['postmenopausal_bleeding']
        )
        self.assertEqual(gyn_res['organ_probabilities'][0]['cancer_code'], 'GYNECOLOGIC')
        self.assertTrue(gyn_res['risk_evaluation']['meets_nice_2ww_referral'])

        # Head & Neck test (Hoarseness + Oral ulceration)
        hn_res = predictor.evaluate_patient(
            age=62,
            gender='M',
            symptoms=['persistent_hoarseness', 'oral_ulceration'],
            smoking_pack_years=35.0
        )
        self.assertEqual(hn_res['organ_probabilities'][0]['cancer_code'], 'HEAD_NECK')

        # Hematologic test (Night sweats + Lymphadenopathy)
        hem_res = predictor.evaluate_patient(
            age=40,
            gender='M',
            symptoms=['night_sweats_fever', 'cervical_lymphadenopathy']
        )
        self.assertEqual(hem_res['organ_probabilities'][0]['cancer_code'], 'HEMATOLOGIC')

        # CNS Brain test
        cns_res = predictor.evaluate_patient(
            age=50,
            gender='M',
            symptoms=['neurological_deficits_headache']
        )
        self.assertEqual(cns_res['organ_probabilities'][0]['cancer_code'], 'CNS_BRAIN')

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

    def test_diagnostic_chatbot_nlp_extraction(self):
        from diagnostic_chatbot import DiagnosticChatbot
        bot = DiagnosticChatbot()
        
        # Test Turkish clinical phrase extraction
        text1 = "Yaklaşık 3 aydır inatçı bir öksürük var ve sabahları balgamımda kan tükürüyorum, ayrıca çok zayıfladım."
        extracted1 = bot.extract_symptoms_from_text(text1)
        self.assertIn("hemoptysis", extracted1)
        self.assertIn("chronic_cough", extracted1)
        self.assertIn("unexplained_weight_loss", extracted1)

        # Test gastrointestinal and breast symptoms
        text2 = "Büyük abdestimde kan gördüm ve katı gıdaları yutmakta zorlanıyorum."
        extracted2 = bot.extract_symptoms_from_text(text2)
        self.assertIn("rectal_bleeding", extracted2)
        self.assertIn("dysphagia", extracted2)

        # Test breast lump and hematuria
        text3 = "Sol mememde elime sert bir kitle geliyor, ayrıca idrar yaparken kırmızı kan fark ettim."
        extracted3 = bot.extract_symptoms_from_text(text3)
        self.assertIn("breast_lump", extracted3)
        self.assertIn("hematuria", extracted3)

        # Test expanded red-flag alarm symptoms (NICE NG12)
        text4 = "Menopoz sonrası vajinal kanama ve lekelenme başladı, karnımda da geçmeyen şişlik ve asit var."
        extracted4 = bot.extract_symptoms_from_text(text4)
        self.assertIn("postmenopausal_bleeding", extracted4)
        self.assertIn("persistent_bloating_ascites", extracted4)

        text5 = "3 haftadır ses kısıklığı var, boynumda sert lenf bezi şişliği fark ettim ve geceleri su gibi terliyorum."
        extracted5 = bot.extract_symptoms_from_text(text5)
        self.assertIn("persistent_hoarseness", extracted5)
        self.assertIn("cervical_lymphadenopathy", extracted5)
        self.assertIn("night_sweats_fever", extracted5)

        text6 = "Ağzımda 4 haftadır geçmeyen yara var, ayrıca sabahları şiddetli baş ağrısı ve kolumda uyuşma felç hissediyorum."
        extracted6 = bot.extract_symptoms_from_text(text6)
        self.assertIn("oral_ulceration", extracted6)
        self.assertIn("neurological_deficits_headache", extracted6)

        text7 = "Bacağımda ve kalçamda geceleri uykudan uyandıran derin kemik ağrısı ve ele gelen sert kemik kitle var."
        extracted7 = bot.extract_symptoms_from_text(text7)
        self.assertIn("bone_pain_mass", extracted7)

    def test_diagnostic_chatbot_multiturn_flow(self):
        from diagnostic_chatbot import DiagnosticChatbot
        bot = DiagnosticChatbot()
        sess_init = bot.start_session("patient_test_01")
        self.assertEqual(sess_init['state'], bot.STATE_SYMPTOMS)
        self.assertIn("Erken Kanser Teşhisi", sess_init['message'])

        # Turn 1: Presenting symptoms
        res1 = bot.process_message("patient_test_01", "2 aydır öksürüyorum ve kanlı balgam var, ayrıca kilo verdim")
        self.assertEqual(res1['state'], bot.STATE_DEMOGRAPHICS)
        self.assertIn("hemoptysis", res1['session_data']['symptoms'])
        self.assertIn("unexplained_weight_loss", res1['session_data']['symptoms'])

        # Turn 2: Demographics
        res2 = bot.process_message("patient_test_01", "58 yaşındayım, erkek")
        self.assertEqual(res2['state'], bot.STATE_SMOKING)
        self.assertEqual(res2['session_data']['age'], 58)
        self.assertEqual(res2['session_data']['gender'], 'M')

        # Turn 3: Smoking pack-years
        res3 = bot.process_message("patient_test_01", "Günde 1 paket, yaklaşık 25 yıl içtim")
        self.assertEqual(res3['state'], bot.STATE_FAMILY)
        self.assertAlmostEqual(res3['session_data']['smoking_pack_years'], 25.0)

        # Turn 4: Family history
        res4 = bot.process_message("patient_test_01", "Evet, babamda akciğer kanseri vardı")
        self.assertEqual(res4['state'], bot.STATE_DURATION)
        self.assertTrue(res4['session_data']['family_history'])

        # Turn 5: Symptom duration -> Final Evaluation
        res5 = bot.process_message("patient_test_01", "Yaklaşık 2 aydır devam ediyor")
        self.assertEqual(res5['state'], bot.STATE_COMPLETED)
        self.assertIn("evaluation", res5)
        eval_data = res5['evaluation']
        risk_eval = eval_data['risk_evaluation']
        self.assertGreater(risk_eval['composite_risk_percent'], 50.0)
        self.assertTrue(risk_eval['meets_nice_2ww_referral'])
        self.assertEqual(eval_data['organ_probabilities'][0]['cancer_code'], 'LUNG')
        self.assertIn("Toraks BT", risk_eval['primary_diagnostic_recommendation'])

    def test_diagnostic_chatbot_risk_generation(self):
        from diagnostic_chatbot import DiagnosticChatbot
        bot = DiagnosticChatbot()
        bot.start_session("patient_eval_test")
        
        # Test low-risk healthy checkup scenario
        bot.process_message("patient_eval_test", "Hiçbir şikayetim yok, sadece rutin kontrol yaptırmak istiyorum")
        bot.process_message("patient_eval_test", "30 yaşındayım, kadın")
        bot.process_message("patient_eval_test", "Hiç sigara içmedim")
        bot.process_message("patient_eval_test", "Hayır, ailemde kanser öyküsü yok")
        final_res = bot.process_message("patient_eval_test", "Şikayetim bulunmamaktadır")
        
        eval_data = final_res['evaluation']
        risk_eval = eval_data['risk_evaluation']
        self.assertLess(risk_eval['composite_risk_percent'], 10.0)
        self.assertFalse(risk_eval['meets_nice_2ww_referral'])
        self.assertEqual(risk_eval['tier_code'], 'LOW')

    def test_diagnostic_chatbot_offline_qa(self):
        """Verifies 100% offline question answering for various oncological inquiries without database dependency."""
        from diagnostic_chatbot import DiagnosticChatbot
        bot = DiagnosticChatbot() # Completely standalone, no db_path

        # 1. Lung cancer questions
        res_lung = bot.process_message("Akciğer kanseri belirtileri nelerdir?")
        self.assertIn("Düşük Doz", res_lung['message'])
        self.assertIn("Toraks BT", res_lung['message'])

        # 2. Colonoscopy inquiry
        res_colo = bot.process_message("Kolonoskopi ne zaman gerekir ve acıtır mı?")
        self.assertIn("Kolonoskopi", res_colo['message'])
        self.assertIn("polip", res_colo['message'].lower())

        # 3. Biopsy myth
        res_biopsy = bot.process_message("Biyopsi yaptırmak kanseri yayar mı?")
        self.assertIn("Yaymaz", res_biopsy['message'])
        self.assertIn("biyopsi", res_biopsy['message'].lower())

        # 4. NICE 2WW referral guidelines
        res_nice = bot.process_message("NICE 2-Week-Wait kuralı nedir?")
        self.assertIn("14", res_nice['message'])
        self.assertIn("%3", res_nice['message'])

        # 5. Hemoptysis alarm
        res_hemo = bot.process_message("Balgamda kan görülmesi tehlikeli mi?")
        self.assertIn("Hemoptizi", res_hemo['message'])
        self.assertIn("Toraks BT", res_hemo['message'])

        # 6. Head & Neck / Hoarseness
        res_hn = bot.process_message("Geçmeyen ses kısıklığı kanser belirtisi olabilir mi?")
        self.assertIn("Larenks", res_hn['message'])

        # 7. Lymphoma / B-symptoms
        res_lymph = bot.process_message("Gece terlemesi kanser işareti midir?")
        self.assertIn("Lenfoma", res_lymph['message'])

        # 8. Gynecologic / Postmenopausal bleeding
        res_gyn = bot.process_message("Menopoz sonrası vajinal kanama tehlikeli midir?")
        self.assertIn("Endometrium", res_gyn['message'])

        # 9. Brain Tumor / Morning headache
        res_cns = bot.process_message("Sabah baş ağrısı ve kusma beyin tümörü belirtisi midir?")
        self.assertIn("Beyin", res_cns['message'])

        # 10. Sarcoma / Bone pain
        res_sarc = bot.process_message("Geceleri derin kemik ağrısı kemik kanseri midir?")
        self.assertIn("Sarkom", res_sarc['message'])

        # 11. Pharmacogenomics DPYD
        res_dpyd = bot.process_message("Kemoterapide DPYD testi neden yapılır?")
        self.assertIn("Fluorouracil", res_dpyd['message'])

        # 12. Synthetic Lethality PARP
        res_sl = bot.process_message("Sentetik letalite ve PARP inhibitörleri nasıl çalışır?")
        self.assertIn("BRCA", res_sl['message'])

        # 13. Immunotherapy & Keytruda
        res_immuno = bot.process_message("İmmünoterapi ve Keytruda nasıl çalışır?")
        self.assertIn("Keytruda", res_immuno['message'])
        self.assertIn("PD-L1", res_immuno['message'])

        # 14. Targeted Therapy & EGFR Osimertinib
        res_targ = bot.process_message("Osimertinib Tagrisso ilacı ne işe yarar?")
        self.assertIn("EGFR", res_targ['message'])

        # 15. PET-CT & SUVmax
        res_pet = bot.process_message("PET-BT çekiminde SUVmax değeri ne anlama gelir?")
        self.assertIn("SUVmax", res_pet['message'])

        # 16. Cancer Staging & Metastasis
        res_stage = bot.process_message("Kanserde Evre 4 ve metastaz ne demektir?")
        self.assertIn("Metastatik", res_stage['message'])

        # 17. HPV Vaccine
        res_hpv = bot.process_message("HPV aşısı Gardasil rahim ağzı kanserini önler mi?")
        self.assertIn("Serviks", res_hpv['message'])

        # 18. Greetings & Gratitude
        res_greet = bot.process_message("Merhaba, günaydın!")
        self.assertIn("Merhaba", res_greet['message'])
        res_thanks = bot.process_message("Çok teşekkür ederim, eline sağlık!")
        self.assertIn("Rica ederim", res_thanks['message'])

        # 19. Adaptive Medical Synthesizer (Unmatched question with anatomical organ & finding)
        res_liver = bot.process_message("Karaciğerimde 3 cm lezyon saptandı ne yapmalıyım?")
        self.assertIn("Karaciğer", res_liver['message'])
        self.assertIn("Gastroenteroloji", res_liver['message'])
        self.assertIn("MRG", res_liver['message'])

    def test_diagnostic_chatbot_dynamic_evaluation_standalone(self):
        """Verifies pure-Python dynamic risk evaluation without external database files."""
        from diagnostic_chatbot import DiagnosticChatbot, clean_for_console
        bot = DiagnosticChatbot()

        # High risk patient (65yo male, smoker, hemoptysis, weight loss)
        eval_high = bot.evaluate_patient_dynamically(
            age=65,
            gender='M',
            symptoms=['hemoptysis', 'unexplained_weight_loss', 'chronic_cough'],
            smoking_pack_years=35.0,
            family_history=True
        )
        self.assertGreater(eval_high['risk_evaluation']['composite_risk_percent'], 60.0)
        self.assertTrue(eval_high['risk_evaluation']['meets_nice_2ww_referral'])
        self.assertEqual(eval_high['organ_probabilities'][0]['cancer_code'], 'LUNG')

        # Clean console test for Windows
        clean_text = clean_for_console("🏥 Test [KLINIK] 🩺")
        self.assertNotIn("🏥", clean_text)
        self.assertIn("[KLINIK]", clean_text)

if __name__ == "__main__":
    unittest.main()

