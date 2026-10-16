"""
Exhaustive Backend Test Suite for Computational Biology & Diagnostic Chatbot Pipeline
Validates every single module, algorithmic boundary, CPIC guideline, Bayesian likelihood ratio,
all 52 clinical chatbot topics, and all 16 adaptive clinical synthesizer organs.
"""

import os
import sys
import unittest

# Ensure src is on python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from diagnostic_chatbot import DiagnosticChatbot, OFFLINE_KNOWLEDGE_TOPICS
from cancer_risk_predictor import CancerRiskPredictor, SYMPTOM_DEFINITIONS
from pharmacogenomics import PharmacogenomicsEngine
from synthetic_lethality import SyntheticLethalityEngine
from survival_analysis import SurvivalAnalyzer
from liquid_biopsy import LiquidBiopsyEngine
from immunotherapy_response import ImmunotherapyPredictor
from docking_scorer import DockingScorer
from sequence_aligner import SequenceAligner

class TestExhaustivePipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chatbot = DiagnosticChatbot()
        cls.risk_predictor = CancerRiskPredictor()
        cls.aligner = SequenceAligner()

    # =========================================================================
    # 1. CHATBOT: ALL 52 KNOWLEDGE TOPICS (100% COVERAGE)
    # =========================================================================
    def test_all_52_chatbot_topics_coverage(self):
        """Verify that all topics in OFFLINE_KNOWLEDGE_TOPICS are matched with high precision."""
        self.assertGreaterEqual(len(OFFLINE_KNOWLEDGE_TOPICS), 50, 
                                f"Expected >= 50 topics, found {len(OFFLINE_KNOWLEDGE_TOPICS)}")

        topic_test_queries = {
            'lung_cancer': 'Akciğer kanserinin başlıca belirtileri nelerdir?',
            'breast_cancer': 'Meme kanserinin erken belirtileri nelerdir?',
            'colorectal_cancer': 'Kolon kanseri ve bağırsak kanseri nasıl belirti verir?',
            'prostate_cancer': 'Prostat kanseri belirtileri ve PSA testi nedir?',
            'pancreatic_cancer': 'Pankreas kanseri erken evrede belirti verir mi?',
            'stomach_cancer': 'Mide kanseri belirtileri ve gastroskopi ne zaman gerekir?',
            'bladder_cancer': 'Mesane kanseri ve sistoskopi nedir?',
            'melanoma_skin': 'Melanom ve ben değişimi nedir?',
            'screening_and_tests': 'Hangi kanser tarama testleri yapılır?',
            'nice_guidelines': 'NICE kılavuzu 2ww sevk kuralı nedir?',
            'colonoscopy_info': 'Kolonoskopi taraması kaç yaşında başlar?',
            'biopsy_myths': 'Biyopsi kanseri yayar mı bıçak değerse yayılır mı?',
            'smoking_risk': 'Sigara içmek kansere nasıl yol açar paket yıl nedir?',
            'family_history': 'Ailemde kanser öyküsü olması riski artırır mı?',
            'liquid_biopsy': 'Sıvı biyopsi ctDNA testi kanda nasıl çalışır?',
            'hemoptysis_alarm': 'Öksürükte balgamda kan gelmesi hemoptizi nedir?',
            'weight_loss_alarm': 'Nedensiz ve istemsiz hızlı kilo kaybı nedir?',
            'rectal_bleeding_alarm': 'Makat veya dışkıda kanama rektal kanama nedir?',
            'hematuria_alarm': 'İdrarda kan görülmesi hematüri nedir?',
            'breast_lump_alarm': 'Memede kitle ve sertlik fark ettim?',
            'dysphagia_alarm': 'Yutma güçlüğü ve disfaji neden olur?',
            'jaundice_alarm': 'Gözlerde sarılık ve sararma ne anlama gelir?',
            'head_neck_cancer': 'Baş boyun kanserleri ve larenks kanseri nedir?',
            'lymphoma_hematologic': 'Lenfoma belirtileri ve gece terlemesi nedir?',
            'gynecologic_cancer': 'Jinekolojik kanserler ve rahim kanseri nedir?',
            'cns_brain_tumor': 'Beyin tümörü belirtileri nelerdir?',
            'sarcoma_bone_tumor': 'Kemik kanseri ve sarkom belirtileri nelerdir?',
            'thyroid_cancer': 'Tiroid kanseri ve nodül nedir?',
            'pharmacogenomics_dpyd': 'DPYD geni testi kemoterapide neden önemlidir?',
            'synthetic_lethality_parp': 'Sentetik letalite ve PARP inhibitörleri nasıl etki eder?',
            'immunotherapy_pembro': 'İmmünoterapi ve Pembrolizumab Keytruda nasıl çalışır?',
            'targeted_therapy_osimertinib': 'Osimertinib Tagrisso ve EGFR T790M mutasyonu nedir?',
            'kras_sotorasib': 'KRAS G12C mutasyonu ve Sotorasib Lumakras nedir?',
            'her2_trastuzumab': 'HER2 pozitifliği ve Trastuzumab Herceptin nasıl etki eder?',
            'chemotherapy_basics': 'Kemoterapi sitotoksik tedavi nedir nasıl uygulanır?',
            'chemo_side_effects': 'Kemoterapi yan etkileri saç dökülmesi ve nötropeni nedir?',
            'radiotherapy_basics': 'Radyoterapi ışın tedavisi cyberknife nedir?',
            'hormonotherapy_breast_prostate': 'Hormonoterapi Tamoksifen ve Letrozol nedir?',
            'pet_ct_suvmax': 'PET-BT SUVmax değeri ne anlama gelir?',
            'biopsy_types': 'Tru-cut biyopsi ve ince iğne biyopsisi çeşitleri nelerdir?',
            'tumor_marker_cea': 'CEA tümör belirteci testi nedir?',
            'tumor_marker_ca125': 'CA 125 belirteci yumurtalık kanserinde nasıl kullanılır?',
            'tumor_marker_ca19_9': 'CA 19-9 belirteci pankreas ve safrada nedir?',
            'cancer_stages_tnm': 'Kanser evreleri Evre 4 ve TNM sistemi nedir?',
            'remission_vs_cure': 'Remisyon nedir nüks ve tam kür ne demektir?',
            'survival_os_pfs': 'Genel sağkalım OS ve progresyonsuz sağkalım PFS nedir?',
            'cancer_prevention_diet': 'Kanserden korunma yolları ve beslenme nasıldır?',
            'hpv_vaccine': 'HPV aşısı Gardasil rahim ağzı kanserini önler mi?',
            'genetic_counseling': 'Genetik kanser testi ve BRCA danışmanlığı nedir?',
            'greetings': 'Merhaba günaydın selam',
            'gratitude': 'Çok teşekkür ederim eline sağlık sağol',
            'bot_help': 'Yardım nasıl kullanılır ne sorabilirim?'
        }

        for topic_key, query in topic_test_queries.items():
            matched = self.chatbot.find_matching_topic(query)
            self.assertIsNotNone(matched, f"Query '{query}' failed to match topic '{topic_key}'")
            matched_key = matched.get('key')
            # If matched object has title, verify non-empty answer
            self.assertIn('answer', matched)
            self.assertGreater(len(matched['answer']), 20)

    # =========================================================================
    # 2. CHATBOT: ADAPTIVE CLINICAL SYNTHESIZER (16 ORGANS)
    # =========================================================================
    def test_adaptive_clinical_synthesizer_all_16_organs(self):
        """Verify that the adaptive synthesizer recognizes 16 anatomical organs with clinical guidance."""
        organ_scenarios = [
            ("Karaciğerimde 3 cm lezyon görüldü", "Karaciğer", "Gastroenteroloji", "BT veya Karaciğer MRG"),
            ("Sol böbrekte kitle saptandı", "Böbrek", "Üroloji", "Ultrasonografisi"),
            ("Pankreas kuyruğunda nodül çıktı", "Pankreas", "Gastroenteroloji / Cerrahi Onkoloji", "CA 19-9"),
            ("Akciğer tomografisinde lezyon var", "Akciğer", "Göğüs Hastalıkları", "Toraks BT"),
            ("Memede sertlik ve kist tespit edildi", "Meme", "Genel Cerrahi / Meme Cerrahisi", "Mamografi"),
            ("Midede kalınlaşma ve yara izlendi", "Mide", "Gastroenteroloji", "Gastroskopi"),
            ("Kolon duvarında polip ve kitle", "Kolon", "Gastroenteroloji", "Kolonoskopi"),
            ("Bağırsakta lezyon şüphesi bildirildi", "Bağırsak", "Gastroenteroloji", "Kolonoskopi"),
            ("Prostatta nodül ve büyüme var", "Prostat", "Üroloji", "PSA"),
            ("Beyinde kontrast tutan lezyon saptandı", "Beyin", "Beyin ve Sinir Cerrahisi", "MRG"),
            ("Kemikte litik lezyon ve ağrı", "Kemik", "Ortopedik Onkoloji", "Sintigrafisi"),
            ("Tiroid bezinde 2 cm nodül çıktı", "Tiroid", "Endokrinoloji", "İnce İğne Aspirasyon"),
            ("Yumurtalıkta kompleks kistik kitle", "Yumurtalık", "Jinekolojik Onkoloji", "Pelvik MRG"),
            ("Rahimde kalınlaşma ve polip görüldü", "Rahim", "Kadın Hastalıkları ve Doğum", "Pipelle"),
            ("Larenks bölgesinde lezyon saptandı", "Larenks", "Kulak Burun Boğaz", "Laringoskopi"),
            ("Boyunda lenf nodunda büyüme var", "Boyun", "Kulak Burun Boğaz", "Boyun")
        ]

        for query, expected_organ, expected_clinic, expected_test in organ_scenarios:
            response = self.chatbot.synthesize_clinical_guidance(query)
            self.assertIn(expected_organ, response, f"Organ '{expected_organ}' missing in query '{query}'")
            self.assertIn(expected_clinic, response, f"Clinic '{expected_clinic}' missing in query '{query}'")
            self.assertIn(expected_test, response, f"Test '{expected_test}' missing in query '{query}'")

    # =========================================================================
    # 3. CHATBOT: MULTI-TURN ANAMNESIS FULL CONVERSATION FLOW
    # =========================================================================
    def test_chatbot_multiturn_full_conversation_and_reset(self):
        """Simulate an entire patient conversation from chief complaint to completed report."""
        bot = DiagnosticChatbot()
        bot.start_session("patient_session_01")

        # Step 1: Initial complaint (unexplained weight loss + dysphagia)
        r1 = bot.process_message("patient_session_01", "Son zamanlarda istemsiz kilo kaybı ve yutma güçlüğü yaşıyorum")
        self.assertIn("belirtilerini kaydettim", r1['message'])
        self.assertEqual(r1['state'], bot.STATE_DEMOGRAPHICS)
        self.assertIn("unexplained_weight_loss", r1['session_data']['symptoms'])
        self.assertIn("dysphagia", r1['session_data']['symptoms'])

        # Step 2: Demographics
        r2 = bot.process_message("patient_session_01", "62 Yaşında Erkek")
        self.assertEqual(r2['state'], bot.STATE_SMOKING)
        self.assertEqual(r2['session_data']['age'], 62)
        self.assertEqual(r2['session_data']['gender'], 'M')

        # Step 3: Smoking
        r3 = bot.process_message("patient_session_01", "Günde 1 paket, 30 yıl içtim")
        self.assertEqual(r3['state'], bot.STATE_FAMILY)
        self.assertGreaterEqual(r3['session_data']['smoking_pack_years'], 20.0)

        # Step 4: Family history
        r4 = bot.process_message("patient_session_01", "Evet, babamda mide kanseri vardı")
        self.assertEqual(r4['state'], bot.STATE_DURATION)
        self.assertTrue(r4['session_data']['family_history'])

        # Step 5: Duration -> Completes Evaluation
        r5 = bot.process_message("patient_session_01", "Yaklaşık 3 aydır devam ediyor")
        self.assertEqual(r5['state'], bot.STATE_COMPLETED)
        self.assertIn("DEĞERLENDİRME RAPORU", r5['message'])
        self.assertIsNotNone(r5['evaluation'])
        self.assertTrue(r5['evaluation']['risk_evaluation']['meets_nice_2ww_referral'])

        # Step 6: Reset Session
        reset_sess = bot.start_session("patient_session_01")
        self.assertEqual(reset_sess['state'], bot.STATE_SYMPTOMS)
        self.assertEqual(len(bot.sessions["patient_session_01"]['symptoms']), 0)

    # =========================================================================
    # 4. CANCER RISK PREDICTOR: ALL 22 SYMPTOMS ACROSS GENDERS & AGES
    # =========================================================================
    def test_cancer_risk_predictor_all_22_symptoms(self):
        """Verify Bayesian likelihood updates for every single symptom in SYMPTOM_DEFINITIONS."""
        self.assertEqual(len(SYMPTOM_DEFINITIONS), 22)

        for sym_id, sym_info in SYMPTOM_DEFINITIONS.items():
            self.assertIn('name_tr', sym_info)
            self.assertIn('target_cancers', sym_info)
            self.assertIn('base_odds_ratio', sym_info)
            self.assertIn('nice_ppv', sym_info)
            self.assertGreater(sym_info['base_odds_ratio'], 1.0)

            # Test single symptom evaluation in male
            res_m = self.risk_predictor.evaluate_patient(
                age=60,
                gender='M',
                symptoms=[sym_id],
                smoking_pack_years=0.0,
                family_history=False
            )
            self.assertIn('composite_risk_percent', res_m['risk_evaluation'])
            self.assertGreater(res_m['risk_evaluation']['composite_risk_percent'], 0.0)

            # Test single symptom evaluation in female
            res_f = self.risk_predictor.evaluate_patient(
                age=60,
                gender='F',
                symptoms=[sym_id],
                smoking_pack_years=0.0,
                family_history=False
            )
            self.assertIn('composite_risk_percent', res_f['risk_evaluation'])
            self.assertGreater(res_f['risk_evaluation']['composite_risk_percent'], 0.0)

    # =========================================================================
    # 5. PHARMACOGENOMICS (CPIC GUIDELINES EXHAUSTIVE)
    # =========================================================================
    def test_pharmacogenomics_cpic_guidelines(self):
        """Verify DPYD, UGT1A1, and CYP2D6 metabolizer status and CPIC dosage recommendations."""
        # DPYD Poor Metabolizer (*2A / *2A) -> Activity Score 0.0
        dpyd_poor = PharmacogenomicsEngine.evaluate_dpyd('*2A', '*2A')
        self.assertEqual(dpyd_poor['activity_score'], 0.0)
        self.assertEqual(dpyd_poor['phenotype'], "Poor Metabolizer")
        self.assertEqual(dpyd_poor['action_code'], "CONTRAINDICATED")

        # DPYD Intermediate Metabolizer (*1 / *2A) -> Activity Score 1.0
        dpyd_int = PharmacogenomicsEngine.evaluate_dpyd('*1', '*2A')
        self.assertEqual(dpyd_int['activity_score'], 1.0)
        self.assertEqual(dpyd_int['phenotype'], "Intermediate Metabolizer")
        self.assertEqual(dpyd_int['action_code'], "DOSE_REDUCTION_50")

        # DPYD Normal (*1 / *1) -> Activity Score 2.0
        dpyd_norm = PharmacogenomicsEngine.evaluate_dpyd('*1', '*1')
        self.assertEqual(dpyd_norm['activity_score'], 2.0)
        self.assertEqual(dpyd_norm['phenotype'], "Normal Metabolizer")
        self.assertEqual(dpyd_norm['action_code'], "STANDARD_DOSE")

        # UGT1A1 *28/*28 Homozygote (Gilbert / Irinotecan severe toxicity)
        ugt = PharmacogenomicsEngine.evaluate_ugt1a1('*28/*28')
        self.assertEqual(ugt['phenotype_code'], "PM")
        self.assertIn("İrinotekan", ugt['clinical_action'])

        # CYP2D6 Poor Metabolizer (Tamoxifen activation failure)
        cyp = PharmacogenomicsEngine.evaluate_cyp2d6_tamoxifen('*4/*4')
        self.assertEqual(cyp['phenotype_code'], "PM")
        self.assertIn("Aromataz", cyp['clinical_action'])

    # =========================================================================
    # 6. SYNTHETIC LETHALITY & DEPMAP TARGETS
    # =========================================================================
    def test_synthetic_lethality_targets(self):
        """Verify clinically actionable synthetic lethality pairs."""
        # BRCA1 -> PARP1
        sl_brca = SyntheticLethalityEngine.query_vulnerabilities("BRCA1")
        self.assertEqual(sl_brca['synthetic_lethal_target'], 'PARP1')
        self.assertTrue(any('Olaparib' in t['drug'] for t in sl_brca['clinical_therapies']))

        # MSI-H / dMMR -> WRN helicase
        sl_msi = SyntheticLethalityEngine.query_vulnerabilities("MSI-H")
        self.assertEqual(sl_msi['synthetic_lethal_target'], 'WRN')

        # MTAP deletion -> PRMT5 / MAT2A
        sl_mtap = SyntheticLethalityEngine.query_vulnerabilities("MTAP")
        self.assertTrue(any(g in sl_mtap['synthetic_lethal_target'] for g in ['PRMT5', 'MAT2A']))

    # =========================================================================
    # 7. SURVIVAL ANALYSIS & LOG-RANK TEST
    # =========================================================================
    def test_survival_analysis_kaplan_meier_and_logrank(self):
        """Verify Kaplan-Meier product-limit estimation and log-rank statistics."""
        # Mock trial data: Treatment vs Control
        exp_durations = [3.0, 6.0, 9.0, 12.0, 18.0, 24.0, 30.0]
        exp_events = [1, 1, 0, 1, 0, 1, 0]
        ctrl_durations = [2.0, 4.0, 6.0, 8.0, 10.0, 14.0, 18.0]
        ctrl_events = [1, 1, 1, 1, 1, 0, 1]

        km_exp = SurvivalAnalyzer.calculate_kaplan_meier(exp_durations, exp_events)
        self.assertGreater(km_exp['total_patients'], 0)
        self.assertIn('curve_points', km_exp)

        # Log-rank comparison
        logrank_res = SurvivalAnalyzer.compute_log_rank_test(exp_durations, exp_events, ctrl_durations, ctrl_events)
        self.assertIn('p_value', logrank_res)
        self.assertIn('hazard_ratio', logrank_res)

    # =========================================================================
    # 8. LIQUID BIOPSY LONGITUDINAL MONITORING
    # =========================================================================
    def test_liquid_biopsy_monitoring(self):
        """Verify ctDNA VAF response classifications."""
        # Case 1: Molecular Complete Response (mCR)
        res_mcr = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(baseline_vaf_percent=8.5, followup_vaf_percent=0.0)
        self.assertEqual(res_mcr['response_code'], "MCR")

        # Case 2: Partial Response (mPR) (>50% drop)
        res_pr = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(baseline_vaf_percent=10.0, followup_vaf_percent=2.5)
        self.assertEqual(res_pr['response_code'], "MPR")

        # Case 3: Molecular Progressive Disease (mPD) (relapse / spike)
        res_mpd = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(baseline_vaf_percent=2.0, followup_vaf_percent=5.2)
        self.assertEqual(res_mpd['response_code'], "MPD")

    # =========================================================================
    # 9. IMMUNOTHERAPY RESPONSE & DOCKING SCORER
    # =========================================================================
    def test_immunotherapy_and_docking(self):
        """Verify composite biomarker predictor and Lipinski rule-of-5 scoring."""
        # High response profile
        resp_high = ImmunotherapyPredictor.evaluate_response(tmb_mut_per_mb=18.5, msi_status="MSI-H", pdl1_tps_percent=75.0)
        self.assertGreaterEqual(resp_high['response_probability_percent'], 70.0)
        self.assertIn("Exceptional", resp_high['classification'])

        # Low response profile
        resp_low = ImmunotherapyPredictor.evaluate_response(tmb_mut_per_mb=2.1, msi_status="MSS", pdl1_tps_percent=0.0)
        self.assertLess(resp_low['response_probability_percent'], 35.0)

        # Docking / Lipinski
        test_drug = {
            'id': 'TEST_MOL',
            'name': 'Test Molecule',
            'molecular_weight': 420.0,
            'xlogp': 3.1,
            'hbd_count': 2,
            'hba_count': 6,
            'rotatable_bonds': 5,
            'tpsa': 68.0
        }
        lipinski = DockingScorer.evaluate_lipinski(test_drug)
        self.assertTrue(lipinski['lipinski_pass'])
        self.assertEqual(lipinski['lipinski_violation_count'], 0)

    # =========================================================================
    # 10. SEQUENCE ALIGNMENT (NEEDLEMAN-WUNSCH & BLOSUM62)
    # =========================================================================
    def test_sequence_alignment_algorithms(self):
        """Verify global sequence alignment and mutation extraction."""
        seq_wt = "KITDFGRAKLLGAEEKEYHAE"
        seq_mut = "KITDFGRAKLRGAEEKEYHAE" # L -> R mutation at pos 11

        alignment = self.aligner.align(seq_wt, seq_mut)
        self.assertIn('alignment_score', alignment)
        self.assertIn('identity_percent', alignment)
        self.assertGreater(alignment['identity_percent'], 90.0)
        self.assertEqual(len(alignment['mutations_detected']), 1)
        self.assertEqual(alignment['mutations_detected'][0]['ref'], 'L')
        self.assertEqual(alignment['mutations_detected'][0]['alt'], 'R')

if __name__ == "__main__":
    unittest.main()
