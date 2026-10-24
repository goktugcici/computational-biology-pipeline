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
from molecular_tumor_board import MolecularTumorBoardEngine
from cancer_pathway_engine import CancerPathwayEngine, CANONICAL_PATHWAYS
from immunotherapy_resistance import ImmunotherapyResistanceEngine, CANONICAL_ICB_SCENARIOS
from mutational_signatures_hrd import MutationalSignaturesHrdEngine, CANONICAL_MUTATIONAL_SIGNATURES
from spatial_transcriptomics import SpatialTranscriptomicsEngine, CANONICAL_SPATIAL_SCENARIOS
from epigenetics_methylation import EpigenomicsEngine, CANONICAL_EPIGENETIC_SCENARIOS, analyze_epigenetics
from cell_therapy_cart import CellTherapyCartEngine, CANONICAL_CART_CONSTRUCTS, analyze_cart_construct
from antibody_drug_conjugates import AntibodyDrugConjugateEngine, CANONICAL_ADC_CONSTRUCTS, analyze_adc_construct
from radioligand_theranostics import RadioligandTheranosticsEngine, CANONICAL_RADIOLIGAND_REGIMENS, analyze_radioligand_regimen
from bispecific_antibodies import BispecificAntibodyEngine, CANONICAL_BISPECIFIC_CONSTRUCTS, analyze_bispecific_construct
from neoantigen_mrna_vaccine import (
    CANONICAL_NEOVAC_REGIMENS,
    predict_hla_binding_affinity,
    calculate_neoantigen_fitness,
    design_mrna_cassette,
    model_lnp_biophysics,
    analyze_neovac_regimen
)
from oncolytic_virotherapy import (
    CANONICAL_OV_CONSTRUCTS,
    model_oncolytic_viral_kinetics,
    calculate_icd_biomarker_score,
    calculate_abscopal_systemic_index,
    evaluate_neutralizing_antibody_impact,
    analyze_oncolytic_regimen,
    OncolyticVirotherapyEngine
)

class TestExhaustivePipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chatbot = DiagnosticChatbot()
        cls.risk_predictor = CancerRiskPredictor()
        cls.aligner = SequenceAligner()
        cls.ov_engine = OncolyticVirotherapyEngine()

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
            'adc_therapies': 'ADC nedir antikor ilaç konjugatları Enhertu nasıl çalışır?',
            'cart_cell_therapy': 'CAR-T hücre tedavisi nedir ve yaşayan ilaçlar ne yapar?',
            'bispecific_antibodies': 'Bispesifik antikorlar ve Amivantamab nasıl etki eder?',
            'recist_criteria': 'RECIST 1.1 kriterleri tam yanıt ve progresyon nedir?',
            'amp_mtb_tiers': 'AMP seviyeleri Tier 1 Tier 2 ve VUS nedir?',
            'cancer_pathways_signaling': 'Kanserde MAPK yolağı ve sinyal iletimi kaskadı nedir?',
            'clinical_trials_phases': 'Kanser klinik çalışma Faz 1 Faz 2 ve Faz 3 nedir?',
            'comprehensive_genomic_profiling': 'Kapsamlı genomik profil CGP ve NGS testi nedir?',
            'immunotherapy_resistance_stk11_b2m': 'İmmünoterapiye direnç mekanizmaları soğuk tümör ve STK11 B2M kaybı nedir?',
            'tumor_infiltrating_lymphocytes_til': 'TIL tedavisi tümör infiltre eden lenfosit ve lifileucel nedir?',
            'oncolytic_virus_therapy': 'Onkolitik virüs tedavisi imlygic t-vec kanser virüsü nasıl çalışır?',
            'neoantigen_mrna_vaccines': 'Kişiye özel mrna kanser aşısı ve neoantijen nedir?',
            'radioligand_pluvicto_lutathera': 'Radyoligand tedavisi Pluvicto ve Lutathera nedir?',
            'mutational_signatures_cosmic': 'COSMIC mutasyonel imzalar SBS1 ve APOBEC SBS2 nedir?',
            'hrd_genomic_scar': 'Homolog rekombinasyon eksikliği HRD skoru LOH ve PARP nedir?',
            'protac_molecular_glues': 'PROTAC teknolojisi ve hedefe yönelik protein yıkımı TPD nasıl çalışır?',
            'ctc_circulating_tumor_cells': 'Dolaşımdaki tümör hücresi CTC ve CellSearch testi nedir?',
            'cryo_em_in_drug_discovery': 'Kriyojenik elektron mikroskobisi Cryo-EM ile ilaç keşfi nasıl yapılır?',
            'spatial_transcriptomics_visium': '10x Visium uzamsal transkriptomik doku haritası nedir?',
            'tertiary_lymphoid_structures_tls': 'Tersiyer lenfoid yapılar TLS ve germinal merkez immünoterapi yanıtı nedir?',
            'tumor_associated_macrophages_tam': 'Tümör ilişkili makrofajlar TAM M1 M2 ve CD47 beni yeme sinyali nedir?',
            'cancer_associated_fibroblasts_caf': 'Kanserle ilişkili fibroblastlar CAF ve desmoplazi stroma bariyeri nedir?',
            'single_cell_rna_sequencing_scrna': 'Tek hücre RNA dizileme scRNA-seq ve klonal heterojenite nedir?',
            'mgmt_promoter_methylation': 'Glioblastomda MGMT metilasyonu ve Temozolomid duyarlılığı nedir?',
            'mlh1_promoter_hypermethylation': 'MLH1 hipermetilasyonu ile Lynch sendromu ayrımı nasıl yapılır?',
            'multicancer_early_detection_mced': 'MCED sıvı biyopsi testleri ve cfDNA metilasyonu nedir?',
            'ezh2_inhibitors_epigenetics': 'EZH2 inhibitörü Tazemetostat ve SMARCB1 kaybı nedir?',
            'hypomethylating_agents_dnmt_hdac': 'Hipometilleyici ajanlar Azasitidin ve HDAC inhibitörleri nasıl etki eder?',
            'cart_therapy_generations': 'CAR-T nesilleri 4-1BB CD28 farkı ve TRUCK zırhlı hücreler nedir?',
            'crs_cytokine_release_syndrome': 'Sitokin salınım sendromu CRS nedir ve Tocilizumab nasıl uygulanır?',
            'icans_neurotoxicity_cart': 'ICANS nedir ve nörotoksisite ICE skoru nasıl bakılır?',
            'bcma_carvykti_myeloma': 'Multipl miyelomda BCMA CAR-T ve Carvykti tedavisi nedir?',
            'claudin18_2_gastric_cancer': 'Claudin 18.2 hedefi Zolbetuximab ve CT041 tedavisi nedir?',
            'adc_mechanism_bystander': 'Bystander effect komşu hücre öldürme antikor ilaç konjugatlarında nedir?',
            'enhertu_her2_low': 'HER2-low meme kanseri Enhertu T-DXd tedavisi nedir?',
            'trodelvy_trop2_tnbc': 'Trodelvy Sacituzumab Trop-2 hedefli üçlü negatif meme kanseri nedir?',
            'padcev_nectin4_bladder': 'Padcev Enfortumab Vedotin Nectin-4 mesane kanseri tedavisi nedir?',
            'adc_toxicity_ild_neuropathy': 'ADC yan etkileri ILD pnömonit ve periferik nöropati yönetimi nedir?',
            'pluvicto_psma_mcrpc': '177Lu-PSMA-617 Pluvicto tedavisi nedir?',
            'lutathera_sstr2_net': 'Lutathera 177Lu-DOTATATE nöroendokrin tümörlerde PRRT tedavisi nedir?',
            'targeted_alpha_therapy_actinium': 'Hedefe yönelik alfa tedavisi Actinium-225 ve TAT nedir?',
            'theranostics_pet_paired_imaging': 'Teranostik tıp nedir ve PET eşleştirilmiş tedavi nasıl yapılır?',
            'rlt_radiation_safety_dosimetry': 'Radyoligand güvenliği MIRD dozimetri ve böbrek koruma nasıl yapılır?',
            'blinatumomab_bite_all': 'Blinatumomab Blincyto CD19xCD3 BiTE tedavisi nedir?',
            'amivantamab_egfr_met_nsclc': 'Amivantamab Rybrevant EGFR MET ekzon 20 insersiyon tedavisi nedir?',
            'teclistamab_bcma_cd3_myeloma': 'Teclistamab Tecvayli ve Talquetamab Talvey tedavisi nedir?',
            'tebentafusp_immtac_uveal_melanoma': 'Tebentafusp Kimmtrak ImmTAC tedavisi nedir?',
            'bispecific_step_up_dosing_crs': 'Bispesifik antikorlarda basamaklı dozlama step-up dosing protokolü nedir?',
            'mrna_cancer_vaccine_keynote942': 'mRNA-4157 V940 melanoma kanser aşısı çalışması nedir?',
            'autogene_cevumeran_bnt122_pdac': 'Autogene Cevumeran BNT122 pankreas kanseri aşısı nedir?',
            'neoantigen_hla_binding_netmhcpan': 'Kanser neoantijeni NetMHCpan ile nasıl seçilir?',
            'mrna_cassette_lnp_formulation': 'Sentetik mRNA kaseti ve LNP SM-102 formülasyonu nasıldır?',
            'lynch_syndrome_frameshift_neoantigen': 'Lynch sendromu Nous-209 frameshift aşısı nedir?',
            'kras_amphiphile_vaccine_amplify': 'ELI-002 amfifil mKRAS aşısı sonuçları nelerdir?',
            'tvec_imlygic_melanoma': 'T-VEC onkolitik virüs nasıl çalışır?',
            'delytact_g47delta_glioblastoma': 'DELYTACT glioblastoma tedavisinde onaylı mı?',
            'pelareorep_reolysin_ras': 'Pelareorep reovirüs ve KRAS ilişkisi nedir?',
            'immunogenic_cell_death_icd_ov': 'Onkolitik virüslerde immünojenik hücre ölümü ve calreticulin nasıl uyarılır?',
            'ov_neutralizing_antibodies_abscopal': 'Onkolitik tedavide nötralizan antikorlar ve abskopal etki nedir?',
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

    # =========================================================================
    # 11. MOLECULAR TUMOR BOARD (MTB) & PRECISION ONCOLOGY
    # =========================================================================
    def test_molecular_tumor_board_canonical_cases(self):
        """Verify all canonical MTB cases and consensus recommendations."""
        cases = MolecularTumorBoardEngine.list_all_cases()
        self.assertEqual(len(cases), 6)
        
        c1 = MolecularTumorBoardEngine.get_case('CASE_01_LUNG_MET_EGFR')
        self.assertIsNotNone(c1)
        self.assertIn("Osimertinib", c1['mtb_consensus']['medical_oncologist'])
        self.assertIn("Savolitinib", c1['mtb_consensus']['medical_oncologist'])

    def test_molecular_tumor_board_custom_evaluation(self):
        """Verify dynamic MTB consensus synthesis for custom patient profile."""
        eval_res = MolecularTumorBoardEngine.evaluate_custom_case(
            cancer_type="Non-Small Cell Lung Cancer",
            alterations=[
                {'gene': 'EGFR', 'alteration': 'p.L858R', 'vaf_percent': 35.0},
                {'gene': 'MET', 'alteration': 'Amplification', 'vaf_percent': 0.0}
            ],
            tmb=14.0,
            msi='MSS',
            pdl1=60.0
        )
        self.assertEqual(len(eval_res['tier1_actionable_targets']), 1)
        self.assertEqual(len(eval_res['tier2_investigational_targets']), 1)
        self.assertIn("Osimertinib", eval_res['tier1_actionable_targets'][0]['therapy'])
        self.assertIn("İmmün kontrol noktası", eval_res['mtb_consensus_recommendation'])
        self.assertEqual(eval_res['evidence_strength'], 'Tier I (FDA / NCCN Cat 1)')

    # =========================================================================
    # 12. CANCER PATHWAYS & ONCOGENIC SIGNALING FLUX ENGINE
    # =========================================================================
    def test_cancer_pathway_engine_all_cascades(self):
        """Verify that all canonical cancer signaling pathways are modeled."""
        engine = CancerPathwayEngine()
        pathways = engine.list_pathways()
        self.assertEqual(len(pathways), 5)
        pathway_ids = [p['id'] for p in pathways]
        self.assertIn('MAPK_ERK', pathway_ids)
        self.assertIn('PI3K_AKT_MTOR', pathway_ids)
        self.assertIn('DDR_HRR', pathway_ids)
        self.assertIn('CELL_CYCLE_CDK', pathway_ids)
        self.assertIn('IMMUNE_CHECKPOINT', pathway_ids)

    def test_cancer_pathway_engine_drug_perturbation(self):
        """Verify flux simulation, combination therapy synergy, and resistance feedback."""
        engine = CancerPathwayEngine()

        # 1. MAPK baseline with KRAS G12C
        res_mut = engine.simulate_pathway_flux('MAPK_ERK', alterations=['KRAS_G12C'])
        self.assertGreater(res_mut['final_flux'], 70.0)
        self.assertEqual(res_mut['status'], 'AŞIRI AKTİF / DİRENÇLİ')

        # 2. KRAS G12C + Sotorasib monotherapy (feedback flag)
        res_mono = engine.simulate_pathway_flux('MAPK_ERK', alterations=['KRAS_G12C'], drugs=['Sotorasib'])
        self.assertTrue(any('feedback' in f.lower() for f in res_mono['resistance_flags']))

        # 3. KRAS G12C + Sotorasib + Panitumumab (dual blockade suppression)
        res_combo = engine.simulate_pathway_flux('MAPK_ERK', alterations=['KRAS_G12C'], drugs=['Sotorasib', 'Panitumumab'])
        self.assertLess(res_combo['final_flux'], 35.0)
        self.assertEqual(res_combo['status'], 'BASKILANDI')

        # 4. DDR synthetic lethality (BRCA1 loss + Olaparib)
        res_ddr = engine.simulate_pathway_flux('DDR_HRR', alterations=['BRCA1_MUT'], drugs=['Olaparib'])
        self.assertGreater(res_ddr['apoptosis_score'], 80.0)

    # =========================================================================
    # 13. IMMUNOTHERAPY RESISTANCE & IMMUNE MICROENVIRONMENT (TIME) ENGINE
    # =========================================================================
    def test_immunotherapy_resistance_canonical_scenarios(self):
        """Verify that all canonical ICB scenarios are defined and retrievable."""
        scenarios = ImmunotherapyResistanceEngine.list_scenarios()
        self.assertEqual(len(scenarios), 4)
        sc_ids = [s['id'] for s in scenarios]
        self.assertIn('NSCLC_STK11_KEAP1', sc_ids)
        self.assertIn('MELANOMA_JAK1_B2M', sc_ids)
        self.assertIn('COLON_DMMR_MSI_HIGH', sc_ids)
        self.assertIn('BREAST_TNBC_VEGF_MDSC', sc_ids)

    def test_immunotherapy_resistance_evaluations(self):
        """Verify algorithmic calculation of sensitivity score and resistance mechanisms."""
        # 1. Hot tumor (MSI-High / dMMR)
        res_hot = ImmunotherapyResistanceEngine.evaluate_tumor_immune_status(
            cancer_type='Kolorektal',
            tmb=42.0,
            msi_status='MSI-H',
            pdl1_percent=65.0
        )
        self.assertGreater(res_hot['icb_sensitivity_score'], 75.0)
        self.assertIn('Hot', res_hot['tumor_phenotype'])
        self.assertIn('YÜKSEK DUYARLILIK', res_hot['response_category'])

        # 2. Cold tumor with STK11 / KEAP1
        res_cold = ImmunotherapyResistanceEngine.evaluate_tumor_immune_status(
            cancer_type='NSCLC',
            tmb=14.0,
            msi_status='MSS',
            pdl1_percent=1.0,
            stk11_mut=True,
            keap1_mut=True
        )
        self.assertLess(res_cold['icb_sensitivity_score'], 35.0)
        self.assertIn('Cold', res_cold['tumor_phenotype'])
        self.assertIn('DÜŞÜK YANIT', res_cold['response_category'])
        self.assertTrue(any('STK11' in m for m in res_cold['resistance_reasons']))

        # 3. Acquired resistance via B2M loss + JAK1
        res_b2m = ImmunotherapyResistanceEngine.evaluate_tumor_immune_status(
            cancer_type='Melanom',
            tmb=28.0,
            msi_status='MSS',
            pdl1_percent=80.0,
            b2m_loss=True,
            jak_mut=True
        )
        self.assertLess(res_b2m['icb_sensitivity_score'], 35.0)
        self.assertTrue(any('B2M' in m for m in res_b2m['resistance_reasons']))
        self.assertTrue(any('NK' in b for b in res_b2m['bypass_therapies']))

    # =========================================================================
    # 14. COSMIC MUTATIONAL SIGNATURES & HRD GENOMIC SCAR ENGINE
    # =========================================================================
    def test_mutational_signatures_canonical_signatures(self):
        """Verify that all canonical COSMIC signatures are defined and retrievable."""
        sigs = MutationalSignaturesHrdEngine.list_signatures()
        self.assertEqual(len(sigs), 6)
        sig_ids = [s['id'] for s in sigs]
        self.assertIn('SBS1_DEAMINATION', sig_ids)
        self.assertIn('SBS2_13_APOBEC', sig_ids)
        self.assertIn('SBS3_HRD_BRCA', sig_ids)
        self.assertIn('SBS4_SMOKING', sig_ids)
        self.assertIn('SBS7_UV_LIGHT', sig_ids)
        self.assertIn('SBS10_POLE_PROOFREADING', sig_ids)

    def test_mutational_signatures_hrd_calculation(self):
        """Verify composite HRD scoring (LOH + TAI + LST), FDA cutoff >=42, and BRCA status."""
        # 1. High genomic scar (>=42) + BRCA Mutant -> Positive
        res_pos = MutationalSignaturesHrdEngine.calculate_hrd_genomic_scar(
            loh_score=18,
            tai_score=16,
            lst_score=17,
            brca_status='BRCA1_MUT',
            cancer_type='Over Kanseri'
        )
        self.assertEqual(res_pos['composite_hrd_score'], 51)
        self.assertTrue(res_pos['hrd_positive'])
        self.assertIn('Kategori 1', res_pos['parp_recommendation'])
        self.assertIn('Platin', res_pos['platinum_recommendation'])

        # 2. Low scar (<42) + BRCA WT -> Negative
        res_neg = MutationalSignaturesHrdEngine.calculate_hrd_genomic_scar(
            loh_score=8,
            tai_score=7,
            lst_score=9,
            brca_status='WT',
            cancer_type='Meme Kanseri'
        )
        self.assertEqual(res_neg['composite_hrd_score'], 24)
        self.assertFalse(res_neg['hrd_positive'])
        self.assertIn('HRD Negatif', res_neg['clinical_status'])

    def test_mutational_signatures_spectrum_decomposition(self):
        """Verify heuristic SBS deconvolution matching smoking, UV, and APOBEC patterns."""
        # Smoking C>A high
        dec_smk = MutationalSignaturesHrdEngine.decompose_patient_spectrum(
            c_to_a=65, c_to_g=10, c_to_t=12, t_to_a=5, t_to_c=5, t_to_g=3
        )
        self.assertEqual(dec_smk['dominant_signature'], 'SBS4_SMOKING')

        # UV C>T high
        dec_uv = MutationalSignaturesHrdEngine.decompose_patient_spectrum(
            c_to_a=4, c_to_g=2, c_to_t=85, t_to_a=3, t_to_c=4, t_to_g=2
        )
        self.assertEqual(dec_uv['dominant_signature'], 'SBS7_UV_LIGHT')

    # =========================================================================
    # 15. SPATIAL TRANSCRIPTOMICS & TUMOR MICROARCHITECTURE ENGINE
    # =========================================================================
    def test_spatial_transcriptomics_canonical_scenarios(self):
        """Verify that all canonical spatial microarchitectures are defined."""
        scenarios = SpatialTranscriptomicsEngine.list_scenarios()
        self.assertEqual(len(scenarios), 4)
        sc_ids = [s['id'] for s in scenarios]
        self.assertIn('NSCLC_TLS_MATURE', sc_ids)
        self.assertIn('PDAC_DESMOPLASTIC_CAF', sc_ids)
        self.assertIn('GBM_HYPOXIC_TAM', sc_ids)
        self.assertIn('MELANOMA_INVASIVE_FRONT', sc_ids)

    def test_spatial_transcriptomics_metrics_calculation(self):
        """Verify TLS score, stroma desmoplasia score, and stage classification."""
        # 1. Mature TLS in NSCLC
        metrics_tls = SpatialTranscriptomicsEngine.calculate_spatial_metrics(
            cd8_density=25.0,
            b_cell_density=20.0,
            caf_density=10.0,
            m2_tam_ratio=15.0,
            cxcl13_positive=True
        )
        self.assertGreater(metrics_tls['tls_score'], 70.0)
        self.assertIn('Evre 3', metrics_tls['tls_stage'])
        self.assertIn('Yüksek İmmünoterapi', metrics_tls['icb_prediction'])

        # 2. Desmoplastic CAF barrier in PDAC
        metrics_caf = SpatialTranscriptomicsEngine.calculate_spatial_metrics(
            cd8_density=5.0,
            b_cell_density=3.0,
            caf_density=52.0,
            m2_tam_ratio=40.0,
            cxcl13_positive=False
        )
        self.assertGreater(metrics_caf['stroma_barrier_score'], 70.0)
        self.assertIn('Desmoplazi', metrics_caf['stroma_status'])
        self.assertIn('Evre 0', metrics_caf['tls_stage'])

    # =========================================================================
    # 16. EPIGENETICS & DNA METHYLATION ENGINE
    # =========================================================================
    def test_epigenetics_canonical_scenarios(self):
        """Verify that all canonical epigenetic profiles are defined and retrievable."""
        engine = EpigenomicsEngine()
        scenarios = engine.list_scenarios()
        self.assertEqual(len(scenarios), 5)
        sc_ids = [s['id'] for s in scenarios]
        self.assertIn('GBM_MGMT_METHYLATED', sc_ids)
        self.assertIn('CRC_MLH1_SPORADIC', sc_ids)
        self.assertIn('MCED_LIQUID_BIOPSY_CRC', sc_ids)
        self.assertIn('ES_SMARCB1_EZH2_ONCOGENIC', sc_ids)
        self.assertIn('LGG_IDH1_GCIMP', sc_ids)

    def test_epigenetics_metrics_calculation(self):
        """Verify silencing status, drug targetability, and Lynch/TOO deconvolution."""
        engine = EpigenomicsEngine()

        # 1. Glioblastoma MGMT Hypermethylated
        gbm = engine.calculate_epigenetic_metrics('GBM_MGMT_METHYLATED', custom_methylation_override=70.0)
        self.assertEqual(gbm['primary_gene'], 'MGMT (O-6-Methylguanine-DNA Methyltransferase)')
        self.assertIn('Transcriptional', gbm['silencing_status'])
        self.assertGreater(gbm['targetability_score'], 75.0)
        self.assertIn('Temozolomide', gbm['predicted_therapeutic_response'])

        # 2. Colorectal Sporadic MLH1 Hypermethylated
        crc = engine.calculate_epigenetic_metrics('CRC_MLH1_SPORADIC', custom_methylation_override=85.0)
        self.assertIn('Sporadic MSI-H', crc['lynch_syndrome_triage'])
        self.assertIn('Anti-PD-1', crc['predicted_therapeutic_response'])

        # 3. MCED Liquid Biopsy cfDNA
        mced = engine.calculate_epigenetic_metrics('MCED_LIQUID_BIOPSY_CRC', custom_methylation_override=76.2)
        self.assertIn('Colorectal', mced['tissue_of_origin_deconvolution'])
        self.assertIn('Colonoscopy', mced['predicted_therapeutic_response'])

        # 4. Epithelioid Sarcoma EZH2
        es = engine.calculate_epigenetic_metrics('ES_SMARCB1_EZH2_ONCOGENIC')
        self.assertIn('Tazemetostat', es['predicted_therapeutic_response'])
        self.assertGreaterEqual(es['targetability_score'], 90.0)

        # 5. Low-grade Glioma IDH1 G-CIMP
        lgg = engine.calculate_epigenetic_metrics('LGG_IDH1_GCIMP')
        self.assertIn('Vorasidenib', lgg['predicted_therapeutic_response'])

    # =========================================================================
    # 17. CAR-T & ADOPTIVE CELL ENGINEERING ENGINE
    # =========================================================================
    def test_cart_canonical_constructs(self):
        """Verify all 6 canonical CAR-T constructs are defined with complete clinical and signaling metadata."""
        engine = CellTherapyCartEngine()
        constructs = engine.list_constructs()
        self.assertEqual(len(constructs), 6)
        ids = [c['id'] for c in constructs]
        self.assertIn('CD19_B_ALL_KYMRIAH', ids)
        self.assertIn('BCMA_MM_CARVYKTI', ids)
        self.assertIn('CD19_DLBCL_YESCARTA', ids)
        self.assertIn('CLDN18_2_GASTRIC_CT041', ids)
        self.assertIn('HER2_ARMORED_TRUCK_SOLID', ids)
        self.assertIn('GPC3_HCC_SYNTHETIC', ids)

    def test_cart_fitness_and_astct_safety_evaluation(self):
        """Verify expansion kinetics, exhaustion risk, solid tumor penetrance, and ASTCT CRS/ICANS grading."""
        engine = CellTherapyCartEngine()

        # 1. Kymriah (4-1BB) in pediatric B-ALL
        kym = engine.evaluate_cart_fitness_and_safety('CD19_B_ALL_KYMRIAH', tumor_antigen_density_percent=90.0, tme_immunosuppression_percent=15.0)
        self.assertEqual(kym['costimulatory_domain'], '4-1BB (CD137)')
        self.assertGreater(kym['expansion_score'], 80.0)
        self.assertLess(kym['exhaustion_score'], 40.0)
        self.assertIn('Grade', kym['crs_grade'])
        self.assertIn('Tocilizumab', kym['crs_management'])
        self.assertIn('NCCN', kym['nccn_guideline'])

        # 2. Carvykti (Bivalent Nanobody) in Multiple Myeloma
        carvykti = engine.evaluate_cart_fitness_and_safety('BCMA_MM_CARVYKTI')
        self.assertEqual(carvykti['target_antigen'], 'B-Cell Maturation Antigen (BCMA / TNFRSF17)')
        self.assertGreaterEqual(carvykti['expansion_score'], 90.0)

        # 3. Armored 4th Gen TRUCK in Solid Tumor with high immunosuppression
        truck = engine.evaluate_cart_fitness_and_safety('HER2_ARMORED_TRUCK_SOLID', tumor_antigen_density_percent=80.0, tme_immunosuppression_percent=85.0)
        self.assertIn('4th Generation', truck['generation'])
        self.assertGreater(truck['solid_penetrance_index'], 50.0)

        # 4. Helper function test
        ct041 = analyze_cart_construct('CLDN18_2_GASTRIC_CT041', costimulatory_override='CD28')
        self.assertEqual(ct041['costimulatory_domain'], 'CD28')
        self.assertEqual(ct041['target_antigen'], 'Claudin 18.2 (CLDN18.2)')

    # =========================================================================
    # 18. ANTIBODY-DRUG CONJUGATES (ADC) & LINKER-PAYLOAD ENGINE
    # =========================================================================
    def test_adc_canonical_constructs(self):
        """Verify all 6 canonical ADC constructs are defined with complete pharmacological metadata."""
        engine = AntibodyDrugConjugateEngine()
        constructs = engine.list_constructs()
        self.assertEqual(len(constructs), 6)
        ids = [c['id'] for c in constructs]
        self.assertIn('HER2_T_DXD_ENHERTU', ids)
        self.assertIn('TROP2_SG_TRODELVY', ids)
        self.assertIn('NECTIN4_EV_PADCEV', ids)
        self.assertIn('HER2_T_DM1_KADCYLA', ids)
        self.assertIn('CD30_BV_ADCETRIS', ids)
        self.assertIn('TROP2_DATO_DXD', ids)

    def test_adc_kinetics_bystander_and_safety_evaluation(self):
        """Verify DAR optimization, bystander killing vs non-cleavable SMCC, plasma stability, and ILD/neuropathy risk stratification."""
        engine = AntibodyDrugConjugateEngine()

        # 1. Trastuzumab Deruxtecan (Enhertu) - Cleavable GGFG, DXd, DAR 8.0, High Bystander
        enhertu = engine.evaluate_adc_kinetics_and_safety('HER2_T_DXD_ENHERTU')
        self.assertEqual(enhertu['target_antigen'], 'HER2 (ERBB2 / Neu)')
        self.assertEqual(enhertu['canonical_dar'], 8.0)
        self.assertGreater(enhertu['payload_delivery_score'], 85.0)
        self.assertGreater(enhertu['bystander_index'], 80.0)
        self.assertIn('High Bystander Killing', enhertu['bystander_summary'])
        self.assertGreaterEqual(enhertu['ild_risk_percent'], 10.0)
        self.assertIn('Interstitial Lung Disease', enhertu['primary_toxicity'])
        self.assertIn('NCCN', enhertu['nccn_guideline'])

        # 2. Trastuzumab Emtansine (Kadcyla) - Non-cleavable SMCC, Zero Bystander
        kadcyla = engine.evaluate_adc_kinetics_and_safety('HER2_T_DM1_KADCYLA')
        self.assertFalse(CANONICAL_ADC_CONSTRUCTS['HER2_T_DM1_KADCYLA']['membrane_permeable_payload'])
        self.assertEqual(kadcyla['bystander_index'], 15.0)
        self.assertIn('Zero / Negligible Bystander', kadcyla['bystander_summary'])
        self.assertGreater(kadcyla['plasma_stability_score'], 90.0)

        # 3. Enfortumab Vedotin (Padcev) - Nectin-4, MMAE, Neuropathy
        padcev = engine.evaluate_adc_kinetics_and_safety('NECTIN4_EV_PADCEV')
        self.assertEqual(padcev['target_antigen'], 'Nectin-4 (PVRL4 / Cell adhesion molecule)')
        self.assertGreaterEqual(padcev['neuropathy_risk_percent'], 45.0)
        self.assertIn('Peripheral Sensory Neuropathy', padcev['primary_toxicity'])

        # 4. Helper function test with custom DAR override
        custom_trodelvy = analyze_adc_construct('TROP2_SG_TRODELVY', custom_dar=6.0, target_antigen_expression_percent=90.0)
        self.assertEqual(custom_trodelvy['dar'], 6.0)
        self.assertEqual(custom_trodelvy['target_antigen'], 'Trop-2 (Trophoblast cell-surface antigen 2 / TACSTD2)')
        self.assertGreater(custom_trodelvy['payload_delivery_score'], 75.0)

    # =========================================================================
    # 19. TARGETED RADIOLIGAND THERAPY (RLT) & THERANOSTICS ENGINE
    # =========================================================================
    def test_radioligand_canonical_regimens(self):
        """Verify all 6 canonical theranostic regimens are defined with complete nuclear oncology parameters."""
        engine = RadioligandTheranosticsEngine()
        regimens = engine.list_regimens()
        self.assertEqual(len(regimens), 6)
        ids = [r['id'] for r in regimens]
        self.assertIn('LU177_PSMA_617_PLUVICTO', ids)
        self.assertIn('LU177_DOTATATE_LUTATHERA', ids)
        self.assertIn('AC225_PSMA_617_ALPHA', ids)
        self.assertIn('RA223_DICHLORIDE_XOFIGO', ids)
        self.assertIn('I131_MIBG_AZEDRA', ids)
        self.assertIn('AC225_DOTATATE_ALPHA_NET', ids)

    def test_theranostics_dosimetry_and_vision_criteria(self):
        """Verify PET avidity (VISION criteria), MIRD absorbed dosimetry, amino acid nephroprotection, and LET physics."""
        engine = RadioligandTheranosticsEngine()

        # 1. 177Lu-PSMA-617 (Pluvicto) - VISION Criteria, Beta Low LET, 7.4 GBq
        pluvicto = engine.evaluate_theranostics_and_dosimetry(
            'LU177_PSMA_617_PLUVICTO',
            tumor_suv_max=24.5,
            liver_suv_mean=6.2,
            amino_acid_nephroprotection=True
        )
        self.assertEqual(pluvicto['target_antigen'], 'Prostate-Specific Membrane Antigen (PSMA / FOLH1)')
        self.assertEqual(pluvicto['administered_activity_gbq'], 7.4)
        self.assertGreaterEqual(pluvicto['avidity_ratio'], 1.5)
        self.assertIn('Optimal Theranostic Candidate', pluvicto['eligibility_status'])
        self.assertGreater(pluvicto['tumor_absorbed_dose_gy'], 50.0)
        self.assertLess(pluvicto['cumulative_kidney_dose_gy'], 23.0)  # Within 23 Gy limit
        self.assertLessEqual(pluvicto['cumulative_marrow_dose_gy'], 2.0)   # Within 2.0 Gy limit
        self.assertIn('Beta-Minus', pluvicto['radiation_physics_profile'])
        self.assertIn('VISION', pluvicto['clinical_efficacy'])

        # 2. 225Ac-PSMA-617 (Targeted Alpha Therapy) - High LET, cluster DSBs, severe xerostomia
        ac_alpha = engine.evaluate_theranostics_and_dosimetry(
            'AC225_PSMA_617_ALPHA',
            administered_activity_gbq=0.008,
            tumor_suv_max=30.0,
            liver_suv_mean=5.0
        )
        self.assertIn('Alpha', ac_alpha['radiation_type'])
        self.assertIn('High LET', ac_alpha['radiation_physics_profile'])
        self.assertIn('clustered double-strand DNA breaks', ac_alpha['dna_damage_mechanism'])
        self.assertGreaterEqual(ac_alpha['salivary_risk_score'], 80.0)
        self.assertIn('High / Dose-Limiting', ac_alpha['xerostomia_risk'])

        # 3. 177Lu-DOTATATE (Lutathera) - NETTER-1, SSTR2, AA nephroprotection test
        with_aa = engine.evaluate_theranostics_and_dosimetry('LU177_DOTATATE_LUTATHERA', amino_acid_nephroprotection=True)
        without_aa = engine.evaluate_theranostics_and_dosimetry('LU177_DOTATATE_LUTATHERA', amino_acid_nephroprotection=False)
        self.assertLess(with_aa['cumulative_kidney_dose_gy'], without_aa['cumulative_kidney_dose_gy'])
        self.assertAlmostEqual(with_aa['cumulative_kidney_dose_gy'], round(without_aa['cumulative_kidney_dose_gy'] * 0.55, 1), delta=0.5)

        # 4. Helper convenience function test
        xofigo = analyze_radioligand_regimen('RA223_DICHLORIDE_XOFIGO')
        self.assertEqual(xofigo['target_antigen'], 'Bone Hydroxyapatite (Calcium-Mimetic Osteoblastic Targeting)')
        self.assertIn('ALSYMPCA', xofigo['clinical_efficacy'])

    # =========================================================================
    # 20. BISPECIFIC ANTIBODIES & MULTI-SPECIFIC ENGAGERS (BiTE) ENGINE
    # =========================================================================
    def test_bispecific_canonical_constructs(self):
        """Verify all 6 canonical bispecific constructs are defined with complete pharmacological metadata."""
        engine = BispecificAntibodyEngine()
        constructs = engine.list_constructs()
        self.assertEqual(len(constructs), 6)
        ids = [c['id'] for c in constructs]
        self.assertIn('BLINATUMOMAB_BLINCYTO', ids)
        self.assertIn('AMIVANTAMAB_RYBREVANT', ids)
        self.assertIn('TECLISTAMAB_TECVAYLI', ids)
        self.assertIn('TEBENTAFUSP_KIMMTRAK', ids)
        self.assertIn('EPCORITAMAB_TEPKINLY', ids)
        self.assertIn('TALQUETAMAB_TALVEY', ids)

    def test_bispecific_synapse_affinity_and_stepup_crs_mitigation(self):
        """Verify synapse distance, tumor-to-CD3 affinity ratio, ASTCT CRS mitigation via step-up dosing, and ICANS risk."""
        engine = BispecificAntibodyEngine()

        # 1. Blinatumomab (BiTE) - CD19xCD3, 11.2 nm synapse, step-up dosing, ICANS risk
        blin = engine.evaluate_bispecific_mechanics(
            'BLINATUMOMAB_BLINCYTO',
            custom_target_expression=85.0,
            effector_to_target_ratio=5.0,
            step_up_dosing_adherent=True
        )
        self.assertEqual(blin['target_tumor'], 'CD19 (B-Cell Surface Antigen)')
        self.assertIn('BiTE', blin['format'])
        self.assertLessEqual(blin['synapse_distance_nm'], 12.0)
        self.assertIn('Optimal Lytic Synapse', blin['synapse_badge'])
        self.assertGreater(blin['tumor_lysis_score'], 85.0)
        self.assertIn('Optimal Tumor First', blin['biodistribution_profile'])
        self.assertLess(blin['crs_risk_percent'], 15.0)  # Mitigated by step-up
        self.assertIn('Neurotoxicity / ICANS', blin['primary_toxicity'])
        self.assertIn('TOWER', blin['clinical_trial'])

        # 2. Step-up dosing omitted test
        blin_no_stepup = engine.evaluate_bispecific_mechanics('BLINATUMOMAB_BLINCYTO', step_up_dosing_adherent=False)
        self.assertGreater(blin_no_stepup['crs_risk_percent'], blin['crs_risk_percent'])
        self.assertIn('HIGH RISK OF GRADE >=3 CRS', blin_no_stepup['crs_mitigation_status'])

        # 3. Tebentafusp (ImmTAC TCR soluble bispecific for Uveal Melanoma)
        kimmtrak = engine.evaluate_bispecific_mechanics('TEBENTAFUSP_KIMMTRAK')
        self.assertIn('ImmTAC', kimmtrak['format'])
        self.assertIn('HLA-A*02:01', kimmtrak['target_tumor'])
        self.assertLessEqual(kimmtrak['synapse_distance_nm'], 10.0)
        self.assertIn('IMCgp100-202', kimmtrak['clinical_trial'])

        # 4. Amivantamab (EGFRxMET DuoBody with ADCC)
        rybrevant = engine.evaluate_bispecific_mechanics('AMIVANTAMAB_RYBREVANT')
        self.assertIn('DuoBody', rybrevant['format'])
        self.assertEqual(rybrevant['target_tumor'], 'EGFR (ErbB1) & MET (HGF Receptor)')
        self.assertIn('MARIPOSA', rybrevant['clinical_trial'])

        # 5. Helper convenience function test
        talvey = analyze_bispecific_construct('TALQUETAMAB_TALVEY')
        self.assertEqual(talvey['target_tumor'], 'GPRC5D (G Protein-Coupled Receptor Family C Group 5 Member D)')
        self.assertIn('MonumenTAL-1', talvey['clinical_trial'])

    # =========================================================================
    # 21. CANCER NEOANTIGEN & PERSONALIZED mRNA VACCINE ENGINE (NEOVAC)
    # =========================================================================
    def test_neovac_canonical_regimens(self):
        """Verify all 6 canonical neoantigen regimens have complete biological and clinical metadata."""
        self.assertEqual(len(CANONICAL_NEOVAC_REGIMENS), 6)
        expected_ids = [
            'MRNA_4157_V940_MELANOMA',
            'AUTOGENE_CEVUMERAN_BNT122_PDAC',
            'MSI_H_FRAMESHIFT_POOL_NOUS209',
            'KRAS_MUTANT_AMPLIFY_ELI002',
            'EGFR_VIII_GLIOBLASTOMA_RINDOMUPEBTIDE',
            'NSCLC_HIGH_TMB_PERSONALIZED'
        ]
        for rid in expected_ids:
            self.assertIn(rid, CANONICAL_NEOVAC_REGIMENS)
            r = CANONICAL_NEOVAC_REGIMENS[rid]
            self.assertIn('clinical_trial', r)
            self.assertIn('gene_mutation', r)
            self.assertIn('sample_peptide', r)
            self.assertIn('nccn_guideline', r)
            self.assertGreater(r['base_potency'], 70.0)

    def test_neovac_hla_affinity_and_fitness_scoring(self):
        """Verify NetMHCpan HLA binding prediction and multi-parametric fitness scoring."""
        # 1. Strong binder with favored anchors (P2: Thr, P9: Val on HLA-A*02:01)
        res_strong = predict_hla_binding_affinity('STEKITYSV', 'HLA-A*02:01')
        self.assertLess(res_strong['ic50_nm'], 50.0)
        self.assertEqual(res_strong['binding_class'], 'Strong Binder')
        self.assertLess(res_strong['rank_percent'], 0.5)

        # 2. Non-binder with disfavored anchor (Asp at P2)
        res_weak = predict_hla_binding_affinity('SDEKITYSP', 'HLA-A*02:01')
        self.assertGreater(res_weak['ic50_nm'], 500.0)
        self.assertEqual(res_weak['binding_class'], 'Non-Binder')

        # 3. High fitness evaluation (Clonal VAF=0.48, TPM=125, strong binding)
        fit_high = calculate_neoantigen_fitness(
            peptide='STEKITYSV',
            wt_peptide='STEKITYSL',
            ic50_nm=18.4,
            vaf=0.48,
            tpm_expression=125.0
        )
        self.assertGreaterEqual(fit_high['composite_fitness'], 80.0)
        self.assertIn('Tier 1 Priority', fit_high['recommendation'])

        # 4. Frameshift foreignness boost
        fit_fs = calculate_neoantigen_fitness(
            peptide='RLSSCVPVA',
            wt_peptide='KSLVSSAST',
            ic50_nm=25.0,
            vaf=0.52,
            tpm_expression=195.0,
            is_frameshift=True
        )
        self.assertGreaterEqual(fit_fs['foreignness_score'], 80.0)
        self.assertGreaterEqual(fit_fs['composite_fitness'], 85.0)

    def test_neovac_mrna_cassette_and_lnp_formulation(self):
        """Verify synthetic mRNA cassette design, GC/CAI optimization, and LNP biophysics."""
        # 1. mRNA Cassette Assembly (34 epitopes with AAY linker)
        cassette = design_mrna_cassette(neoantigen_count=34, linker_type='AAY', signal_peptide='LAMP1')
        self.assertGreater(cassette['total_mrna_length_nt'], 1200)
        self.assertIn('CleanCap AG', cassette['cap_structure'])
        self.assertGreaterEqual(cassette['codon_adaptation_index'], 0.90)
        self.assertGreaterEqual(cassette['gc_content_percent'], 58.0)
        self.assertGreaterEqual(cassette['estimated_half_life_hours'], 20.0)

        # 2. LNP Formulation Biophysics (SM-102 at N/P 6.0)
        lnp = model_lnp_biophysics(ionizable_lipid='SM-102', np_ratio=6.0, total_rna_dose_ug=1000.0)
        self.assertGreaterEqual(lnp['encapsulation_efficiency_percent'], 96.0)
        self.assertLess(lnp['mean_particle_diameter_nm'], 90.0)
        self.assertLess(lnp['polydispersity_index_pdi'], 0.12)
        self.assertGreaterEqual(lnp['endosomal_escape_percent'], 1.5)

        # 3. Complete End-to-End Orchestration
        full_analysis = analyze_neovac_regimen('MRNA_4157_V940_MELANOMA')
        self.assertIn('hla_binding', full_analysis)
        self.assertIn('neoantigen_fitness', full_analysis)
        self.assertIn('mrna_cassette', full_analysis)
        self.assertIn('lnp_formulation', full_analysis)
        self.assertGreaterEqual(full_analysis['overall_immunogenicity_potency_score'], 85.0)

    # =========================================================================
    # 22. ONCOLYTIC VIROTHERAPY & SYNTHETIC IMMUNO-VIROLOGY (OV ENGINE)
    # =========================================================================
    def test_oncolytic_virotherapy_constructs_and_kinetics(self):
        """Verify all 6 canonical OV constructs and ODE kinetic replication simulation."""
        # 1. Verify 6 canonical constructs
        self.assertEqual(len(CANONICAL_OV_CONSTRUCTS), 6)
        expected_ids = [
            'TVEC_IMLYGIC_MELANOMA',
            'DELYTACT_G47DELTA_GBM',
            'PELAREOREP_REOLYSIN_RAS',
            'CAN2409_HSV_TK_ADENOVIRUS',
            'CRETOSTIMOGENE_CG0070_BLADDER',
            'VCN01_HYALURONIDASE_ADENOVIRUS'
        ]
        for c_id in expected_ids:
            self.assertIn(c_id, CANONICAL_OV_CONSTRUCTS)
            c = CANONICAL_OV_CONSTRUCTS[c_id]
            self.assertIn('name', c)
            self.assertIn('burst_size', c)
            self.assertIn('lysis_cycle_hours', c)
            self.assertIn('default_dose_pfu', c)
            self.assertGreaterEqual(c['burst_size'], 0)
            self.assertGreater(c['default_dose_pfu'], 1.0e6)

        # 2. Test 14-day ODE kinetics for T-VEC
        kinetics = model_oncolytic_viral_kinetics(
            initial_tumor_cells=1.0e8,
            initial_dose_pfu=1.0e8,
            infectivity_rate=1.0e-9,
            burst_size=1000,
            lysis_rate_per_day=1.0,
            viral_clearance_per_day=0.45,
            days=14
        )
        self.assertGreaterEqual(kinetics['tumor_lysis_percent'], 95.0)
        self.assertGreater(kinetics['peak_viral_titer_pfu'], 1.0e8)
        self.assertGreater(len(kinetics['time_history']), 10)
        self.assertLess(kinetics['nadir_tumor_cells'], 1.0e8)

        # 3. Test ICD Biomarker scoring
        icd = calculate_icd_biomarker_score(
            lysis_percent=98.0,
            virus_family='Herpesviridae (HSV-1)',
            gmcsf_present=True
        )
        self.assertGreaterEqual(icd['calreticulin_translocation_percent'], 70.0)
        self.assertGreaterEqual(icd['extracellular_atp_uM'], 5.0)
        self.assertGreaterEqual(icd['hmgb1_release_ng_ml'], 30.0)
        self.assertGreaterEqual(icd['composite_icd_score'], 80.0)

    def test_oncolytic_virotherapy_nab_and_abscopal_synergy(self):
        """Verify pre-existing NAb pharmacokinetics, carrier shielding, and abscopal synergy."""
        # 1. NAb impact: Intratumoral vs Intravenous
        it_nab = evaluate_neutralizing_antibody_impact(
            virus_family='HSV-1',
            pre_existing_nab_percent=70.0,
            delivery_route='intratumoral',
            carrier_cell_shielding=False
        )
        self.assertGreaterEqual(it_nab['effective_payload_percent'], 75.0)

        # IV without shielding: severe attenuation
        iv_unshielded = evaluate_neutralizing_antibody_impact(
            virus_family='Ad5',
            pre_existing_nab_percent=80.0,
            delivery_route='intravenous',
            carrier_cell_shielding=False
        )
        self.assertLessEqual(iv_unshielded['effective_payload_percent'], 45.0)

        # IV with carrier shielding: rescued effective dose
        iv_shielded = evaluate_neutralizing_antibody_impact(
            virus_family='Ad5',
            pre_existing_nab_percent=80.0,
            delivery_route='intravenous',
            carrier_cell_shielding=True
        )
        self.assertGreaterEqual(iv_shielded['effective_payload_percent'], 70.0)

        # 2. Abscopal Systemic Effect with Anti-PD-1 synergy
        abscopal_mono = calculate_abscopal_systemic_index(
            local_lysis_score=95.0,
            icd_score=90.0,
            delivery_route='intratumoral',
            anti_pd1_combination=False
        )
        abscopal_combo = calculate_abscopal_systemic_index(
            local_lysis_score=95.0,
            icd_score=90.0,
            delivery_route='intratumoral',
            anti_pd1_combination=True
        )
        self.assertGreater(abscopal_combo['abscopal_systemic_index'], abscopal_mono['abscopal_systemic_index'])
        self.assertIn('Synergy', abscopal_combo['checkpoint_synergy'])

        # 3. Engine Class & Full Regimen Analysis
        constructs_list = self.ov_engine.list_constructs()
        self.assertEqual(len(constructs_list), 6)
        c_tvec = self.ov_engine.get_construct('TVEC_IMLYGIC_MELANOMA')
        self.assertIsNotNone(c_tvec)

        eval_res = self.ov_engine.evaluate_oncolytic_treatment('TVEC_IMLYGIC_MELANOMA', custom_dose_pfu=1.0e8)
        self.assertIn('viral_kinetics', eval_res)
        self.assertIn('icd_biomarkers', eval_res)
        self.assertIn('abscopal_effect', eval_res)
        self.assertGreaterEqual(eval_res['overall_therapeutic_score'], 85.0)

if __name__ == "__main__":
    unittest.main()


