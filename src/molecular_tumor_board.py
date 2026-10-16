"""
Molecular Tumor Board (MTB) & Precision Oncology Decision Support Engine.
Implements AMP/ASCO/CAP biomarker evidence tiering (Tier I to IV),
actionable alteration matching, multi-disciplinary consensus synthesis,
and clinical trial matching (NCCN & ESMO Precision Oncology guidelines).
"""

from typing import Dict, List, Any, Optional

# AMP/ASCO/CAP Biomarker Evidence Tiers:
# Tier I: Biomarkers of Strong Clinical Significance (FDA-approved therapy or NCCN Cat 1)
# Tier II: Biomarkers of Potential Clinical Significance (Off-label, investigational, or Phase 2/3 trials)
# Tier III: Biomarkers of Unknown Clinical Significance (VUS)
# Tier IV: Benign or Likely Benign Alterations

CANONICAL_MTB_CASES: Dict[str, Dict[str, Any]] = {
    'CASE_01_LUNG_MET_EGFR': {
        'case_id': 'MTB-2026-001',
        'title': 'İlerlemiş Akciğer Adenokarsinomu (EGFR L858R + De Novo MET Amplifikasyonu)',
        'demographics': {'age': 62, 'gender': 'Kadın', 'smoking': 'Hiç içmemiş (0 Paket-Yıl)', 'ecog_ps': 1},
        'primary_diagnosis': 'Akciğer Adenokarsinomu (Evre IVB - Akciğer içi ve plevral yayılım)',
        'histopathology': 'TTF-1 pozitif, p40 negatif invaziv adenokarsinom',
        'genomic_alterations': [
            {'gene': 'EGFR', 'alteration': 'p.L858R (c.2573T>G, Exon 21)', 'vaf_percent': 38.5, 'tier': 'Tier I', 'implication': '1-3. Kuşak TKI Duyarlılığı'},
            {'gene': 'MET', 'alteration': 'MET Gene Copy Number Gain (CN >= 6.2)', 'vaf_percent': None, 'tier': 'Tier I/II', 'implication': 'Primer TKI Direnci Baypası'},
            {'gene': 'TP53', 'alteration': 'p.R273H', 'vaf_percent': 42.0, 'tier': 'Tier II', 'implication': 'Olumsuz Prognostik Gösterge'}
        ],
        'biomarkers': {'tmb_mut_per_mb': 5.2, 'msi_status': 'MSS', 'pdl1_tps_percent': 15.0},
        'mtb_consensus': {
            'molecular_pathologist': 'EGFR L858R ile eş zamanlı primer MET amplifikasyonu mevcuttur. Standart osimertinib monoterapisi durumunda 4-6 ay içinde erken progresyon riski yüksektir.',
            'medical_oncologist': 'Birinci basamakta Osimertinib (80 mg) + Savolitinib (300-600 mg) kombinasyonu (SAVANNAH / TATTON protokolü) veya Platin+Pemetreksed kemoterapisi eşliğinde TKI değerlendirilmelidir.',
            'clinical_geneticist': 'Hasta sigara içmemiş genç-orta yaş kadındır. Germline BRCA/Lynch paneli temizdir; somatik sürücü mutasyonlar hedeflenmelidir.',
            'guideline_recommendation': 'NCCN Kategori 2A / ESMO Tier I-B: EGFR-TKI + MET İnhibitörü Kombinasyonu',
            'relevant_trials': [{'nct': 'NCT03778229', 'name': 'SAVANNAH (Osimertinib + Savolitinib)', 'phase': 'Phase 2'}]
        }
    },
    'CASE_02_COLON_KRAS_G12C': {
        'case_id': 'MTB-2026-002',
        'title': 'Metastatik Kolorektal Kanser (KRAS G12C, Sol Kolon, MSS)',
        'demographics': {'age': 54, 'gender': 'Erkek', 'smoking': '15 Paket-Yıl', 'ecog_ps': 0},
        'primary_diagnosis': 'Kolon Adenokarsinomu (Evre IV - Karaciğer oligometastatik)',
        'histopathology': 'Orta diferansiye sol kolon adenokarsinomu, CDX-2 pozitif',
        'genomic_alterations': [
            {'gene': 'KRAS', 'alteration': 'p.G12C (c.34G>T, Exon 2)', 'vaf_percent': 31.2, 'tier': 'Tier I', 'implication': 'Anti-EGFR Tek Ajan Direnci, G12C Hedeflenebilir'},
            {'gene': 'APC', 'alteration': 'p.R1450* (Stop-gain)', 'vaf_percent': 34.0, 'tier': 'Tier III', 'implication': 'Wnt Yolağı Aktivasyonu'}
        ],
        'biomarkers': {'tmb_mut_per_mb': 6.8, 'msi_status': 'MSS (pMMR)', 'pdl1_tps_percent': 0.0},
        'mtb_consensus': {
            'molecular_pathologist': 'KRAS G12C kovalent cebi hedeflenebilir yapıdadır; ancak kolorektal kanserde EGFR aracılı hızlı feedback aktivasyonu nedeniyle tek ajan Sotorasib/Adagrasib yanıt oranı düşüktür (%10-15).',
            'medical_oncologist': 'CodeBreak-300 çalışması doğrultusunda Sotorasib + Panitumumab (Anti-EGFR) ikili blokajı veya FOLFIRI + Bevacizumab sonrası ikinci basamakta kovalent kombinasyon önerilir.',
            'clinical_geneticist': 'MSS fenotipinde Lynch şüphesi yoktur. Aile taraması gerektirmez.',
            'guideline_recommendation': 'NCCN Kategori 1: Sotorasib + Panitumumab veya Adagrasib + Setuksimab Kombinasyonu',
            'relevant_trials': [{'nct': 'NCT05192486', 'name': 'CodeBreak-300 (Sotorasib + Panitumumab)', 'phase': 'Phase 3'}]
        }
    },
    'CASE_03_TNBC_BRCA1_PDL1': {
        'case_id': 'MTB-2026-003',
        'title': 'Üçlü Negatif Meme Kanseri (gBRCA1+, PD-L1 CPS 18, TMB 12)',
        'demographics': {'age': 48, 'gender': 'Kadın', 'smoking': 'Hiç içmemiş', 'ecog_ps': 0},
        'primary_diagnosis': 'İnvaziv Duktal Meme Karsinomu (Evre III - T2N2M0, ER-, PR-, HER2- 0)',
        'histopathology': 'Grade 3, yüksek Ki-67 (%75), lenfovasküler invazyon pozitif',
        'genomic_alterations': [
            {'gene': 'BRCA1', 'alteration': 'c.68_69delAG (p.Glu23fs)', 'vaf_percent': 50.1, 'tier': 'Tier I (Germline)', 'implication': 'HRD Tamir Kusuru & PARP Duyarlılığı'},
            {'gene': 'TP53', 'alteration': 'p.R175H', 'vaf_percent': 48.0, 'tier': 'Tier II', 'implication': 'Apoptoz Direnci'}
        ],
        'biomarkers': {'tmb_mut_per_mb': 12.4, 'msi_status': 'MSS', 'pdl1_tps_percent': 18.0},
        'mtb_consensus': {
            'molecular_pathologist': 'Homolog rekombinasyon eksikliği (HRD) skoru yüksektir (>42). PD-L1 CPS >= 10 ve TMB-High eş zamanlı immünjenik ortam sağlar.',
            'medical_oncologist': 'Neoadjuvan dönemde KEYNOTE-522 protokolü (Karboplatin+Paklitaksel+Doksorubisin+Siklofosfamid + Pembrolizumab). Cerrahi sonrasında kalıntı hastalık varsa OlympiA protokolü ile 1 yıl adjuvan Olaparib.',
            'clinical_geneticist': 'Kalıtsal BRCA1 patojenik varyantı nedeniyle hastanın 1. derece akrabalarına (kız kardeş, çocukları) genetik test ve presemptomatik izlem (yıllık Meme MRG ve profilaktik BSO danışmanlığı).',
            'guideline_recommendation': 'NCCN Kategori 1 (KEYNOTE-522 & OlympiA): İmmünoterapi + Kemo ardından Adjuvan PARP İnhibitörü',
            'relevant_trials': [{'nct': 'NCT03036488', 'name': 'KEYNOTE-522 (Pembrolizumab in TNBC)', 'phase': 'Phase 3'}]
        }
    },
    'CASE_04_PROSTATE_HRD_CDK12': {
        'case_id': 'MTB-2026-004',
        'title': 'Metastatik Kastrasyona Dirençli Prostat Kanseri (mCRPC, CDK12+)',
        'demographics': {'age': 70, 'gender': 'Erkek', 'smoking': '20 Paket-Yıl', 'ecog_ps': 1},
        'primary_diagnosis': 'Prostat Adenokarsinomu (Evre IVB - Kemik ve Pelvik Lenf Nodu Metastazı)',
        'histopathology': 'Gleason Skoru 4+5 = 9 (Grup 5), Akut idrar retansiyonu ve yaygın osteoblastik lezyonlar',
        'genomic_alterations': [
            {'gene': 'CDK12', 'alteration': 'p.E462* (Çift allelik inaktivasyon)', 'vaf_percent': 44.0, 'tier': 'Tier I', 'implication': 'Tandem Duplikasyon & HRD Benzeri Fenotip'},
            {'gene': 'AR', 'alteration': 'AR Copy Number Gain (Amplifikasyon)', 'vaf_percent': None, 'tier': 'Tier II', 'implication': 'Birinci Kuşak Androjen Blokajı Direnci'}
        ],
        'biomarkers': {'tmb_mut_per_mb': 11.0, 'msi_status': 'MSS', 'pdl1_tps_percent': 5.0},
        'mtb_consensus': {
            'molecular_pathologist': 'CDK12 kaybı, genom genelinde fokal tandem duplikasyonlara ve neoantijen yükünde artışa yol açar. DNA tamir defekti taşır.',
            'medical_oncologist': 'TALAPRO-2 ve PROpel çalışmaları temelinde Enzalutamid veya Abirateron ile kombine PARP inhibitörü (Talazoparib veya Olaparib). Ayrıca CDK12 mutasyonlu olgular immünoterapiye (Anti-PD-1) de kısmi duyarlılık gösterebilir.',
            'clinical_geneticist': 'CDK12 germline taşıyıcılığı nadir de olsa dışlanmalı, germline kanda doğrulanmalıdır.',
            'guideline_recommendation': 'NCCN Kategori 1: PARP İnhibitörü + Yeni Nesil Androjen Reseptör İnhibitörü Kombinasyonu',
            'relevant_trials': [{'nct': 'NCT03395197', 'name': 'TALAPRO-2 (Talazoparib + Enzalutamide)', 'phase': 'Phase 3'}]
        }
    },
    'CASE_05_GBM_MGMT_METHYLATED': {
        'case_id': 'MTB-2026-005',
        'title': 'Glioblastoma Multiforme (GBM Evre IV, IDH-wt, MGMT Promotörü Metile)',
        'demographics': {'age': 58, 'gender': 'Kadın', 'smoking': '5 Paket-Yıl', 'ecog_ps': 1},
        'primary_diagnosis': 'Glioblastoma (WHO Grade 4 - Sol Temporal Lob)',
        'histopathology': 'Mikrovasküler proliferasyon, psödo-palizatlanan nekroz, GFAP pozitif',
        'genomic_alterations': [
            {'gene': 'IDH1', 'alteration': 'Wild-Type (R132 Negatif)', 'vaf_percent': None, 'tier': 'Tier I', 'implication': 'Klasik Primer Glioblastoma Fenotipi'},
            {'gene': 'MGMT', 'alteration': 'Promoter Hypermethylation (>45%)', 'vaf_percent': None, 'tier': 'Tier I', 'implication': 'Alkilleyici Ajan (Temozolomid) Duyarlılığı'},
            {'gene': 'TERT', 'alteration': 'p.C228T (Promoter)', 'vaf_percent': 46.5, 'tier': 'Tier II', 'implication': 'Telomeraz Aktivasyonu'}
        ],
        'biomarkers': {'tmb_mut_per_mb': 3.1, 'msi_status': 'MSS', 'pdl1_tps_percent': 0.0},
        'mtb_consensus': {
            'molecular_pathologist': 'MGMT metilasyonu, DNA tamir enzimi olan MGMT\'nin susturulduğunu kanıtlar. Bu durum alkilleyici kemoterapiye (Temozolomid) belirgin duyarlılık öngörür.',
            'medical_oncologist': 'Stupp Protokolü: Maksimal cerrahi rezeksiyon sonrası 6 hafta eşzamanlı Radyoterapi (60 Gy) + Temozolomid (75 mg/m²), ardından 6 siklus adjuvan Temozolomid ve TTFields (Optune) cihazı.',
            'clinical_geneticist': 'Kalıtsal sendrom (Li-Fraumeni vb.) düşündürecek aile öyküsü yoktur.',
            'guideline_recommendation': 'NCCN Kategori 1: Radyoterapi + Temozolomid + Tümör Tedavi Alanları (TTFields)',
            'relevant_trials': [{'nct': 'NCT00916409', 'name': 'EF-14 (TTFields + Temozolomide)', 'phase': 'Phase 3'}]
        }
    },
    'CASE_06_MELANOMA_BRAF_V600E': {
        'case_id': 'MTB-2026-006',
        'title': 'Metastatik Malign Melanom (BRAF V600E, LDH Yüksek, Evre IV)',
        'demographics': {'age': 52, 'gender': 'Erkek', 'smoking': '10 Paket-Yıl', 'ecog_ps': 0},
        'primary_diagnosis': 'Kutanöz Malign Melanom (Evre IV M1c - Karaciğer ve Akciğer Tutulumu)',
        'histopathology': 'S100+, Melan-A+, HMB-45+ iğsi hücreli melanom',
        'genomic_alterations': [
            {'gene': 'BRAF', 'alteration': 'p.V600E (c.1799T>A, Exon 15)', 'vaf_percent': 42.0, 'tier': 'Tier I', 'implication': 'BRAF+MEK İnhibitörü Hedeflenebilir'},
            {'gene': 'CDKN2A', 'alteration': 'Homozygous Deletion', 'vaf_percent': None, 'tier': 'Tier II', 'implication': 'Hücre Döngüsü Kontrol Kaybı'}
        ],
        'biomarkers': {'tmb_mut_per_mb': 24.5, 'msi_status': 'MSS', 'pdl1_tps_percent': 45.0},
        'mtb_consensus': {
            'molecular_pathologist': 'BRAF V600E ile birlikte TMB-High (>20 mut/Mb) ve yüksek PD-L1 mevcuttur. Hem hedefe yönelik TKI hem de immünoterapi potansiyeli çok yüksektir.',
            'medical_oncologist': 'DREAMseq çalışması uyarınca: Yüksek tümör yükü ve LDH yüksekliğinde acil yanıt için Dabrafenib + Trametinib ile indüksiyon, veya uzun vadeli kalıcı küratif yanıt için Nivolumab + Ipilimumab ikili immünoterapisi.',
            'clinical_geneticist': 'CDKN2A germline mutasyonları melanom ve pankreas kanseri yatkınlığı yapabilir; germline kontrolü tavsiye edilir.',
            'guideline_recommendation': 'NCCN Kategori 1: Nivolumab + Ipilimumab veya Dabrafenib + Trametinib',
            'relevant_trials': [{'nct': 'NCT02224781', 'name': 'DREAMseq (Sequencing Targeted vs IO)', 'phase': 'Phase 3'}]
        }
    }
}


class MolecularTumorBoardEngine:
    """
    Precision Oncology Molecular Tumor Board (MTB) Decision Support Engine.
    Evaluates complex patient cases, assigns AMP/ASCO/CAP evidence tiers,
    synthesizes multi-disciplinary consensus, and matches clinical trials.
    """

    @classmethod
    def get_case(cls, case_key: str) -> Optional[Dict[str, Any]]:
        """Retrieves cataloged MTB clinical case by key."""
        return CANONICAL_MTB_CASES.get(case_key)

    @classmethod
    def list_all_cases(cls) -> List[Dict[str, Any]]:
        """Lists summary metadata of all preloaded MTB cases."""
        return [
            {
                'key': k,
                'case_id': v['case_id'],
                'title': v['title'],
                'primary_diagnosis': v['primary_diagnosis'],
                'primary_alteration': f"{v['genomic_alterations'][0]['gene']} {v['genomic_alterations'][0]['alteration']}",
                'guideline': v['mtb_consensus']['guideline_recommendation']
            }
            for k, v in CANONICAL_MTB_CASES.items()
        ]

    @classmethod
    def evaluate_custom_case(cls, cancer_type: str, alterations: List[Dict[str, Any]],
                              tmb: float = 5.0, msi: str = 'MSS', pdl1: float = 0.0) -> Dict[str, Any]:
        """
        Dynamically analyzes custom patient genomic profile and synthesizes
        an automated Molecular Tumor Board (MTB) consensus report.
        """
        tier1_targets = []
        tier2_targets = []
        recommendations = []
        trials = []

        for alt in alterations:
            gene = alt.get('gene', '').upper()
            var = alt.get('alteration', '')
            vaf = alt.get('vaf_percent', 0.0)

            if gene == 'EGFR':
                tier1_targets.append({'gene': 'EGFR', 'alteration': var, 'vaf': vaf, 'therapy': 'Osimertinib (Tagrisso) 80mg oral'})
                recommendations.append("EGFR aktive edici mutasyon için 3. kuşak TKI (Osimertinib) standart birinci basamak tedavidir.")
            elif gene == 'KRAS':
                if 'G12C' in var.upper():
                    tier1_targets.append({'gene': 'KRAS', 'alteration': var, 'vaf': vaf, 'therapy': 'Sotorasib veya Adagrasib + Anti-EGFR (Kolon) / Monoterapi (Akciğer)'})
                    recommendations.append("KRAS G12C kovalent inhibitörü endikedir.")
                else:
                    tier2_targets.append({'gene': 'KRAS', 'alteration': var, 'vaf': vaf, 'therapy': 'MEK/SHP2 İnhibitör Denemeleri'})
                    recommendations.append(f"KRAS {var} klasik kodon mutasyonu; anti-EGFR monoklonal antikorlarına direnç gösterir.")
            elif gene in ['BRCA1', 'BRCA2']:
                tier1_targets.append({'gene': gene, 'alteration': var, 'vaf': vaf, 'therapy': 'PARP İnhibitörü (Olaparib, Rucaparib, Talazoparib)'})
                recommendations.append(f"{gene} kaybı sentetik letalite yoluyla PARP inhibitörlerine yüksek duyarlılık sağlar.")
            elif gene == 'BRAF' and 'V600' in var.upper():
                tier1_targets.append({'gene': 'BRAF', 'alteration': var, 'vaf': vaf, 'therapy': 'Dabrafenib + Trametinib veya Enkorafenib + Setuksimab'})
                recommendations.append("BRAF V600E mutasyonu için kombine BRAF + MEK blokajı gereklidir.")
            elif gene == 'MET' and 'AMP' in var.upper():
                tier2_targets.append({'gene': 'MET', 'alteration': var, 'vaf': vaf, 'therapy': 'Savolitinib / Tepotinib'})
                recommendations.append("MET amplifikasyonu TKI baypas direnci yaratır; MET inhibitörü eklenmelidir.")
            elif gene == 'TP53':
                tier2_targets.append({'gene': 'TP53', 'alteration': var, 'vaf': vaf, 'therapy': 'Standart Kemoterapi / WEE1 İnhibitörleri'})

        # Biomarker evaluations
        is_tmb_high = tmb >= 10.0
        is_msi_high = 'MSI-H' in msi.upper() or 'DMMR' in msi.upper()
        if is_msi_high or is_tmb_high or pdl1 >= 50.0:
            recommendations.append(f"İmmün kontrol noktası blokajı (Pembrolizumab / Nivolumab) endikedir (TMB: {tmb}, MSI: {msi}, PD-L1: %{pdl1}).")

        consensus_summary = " | ".join(recommendations) if recommendations else "Standart histolojik rehberlere uygun tedavi planlanmalıdır."

        return {
            'cancer_type': cancer_type,
            'biomarker_profile': {'tmb': tmb, 'msi': msi, 'pdl1_percent': pdl1, 'is_tmb_high': is_tmb_high, 'is_msi_high': is_msi_high},
            'tier1_actionable_targets': tier1_targets,
            'tier2_investigational_targets': tier2_targets,
            'mtb_consensus_recommendation': consensus_summary,
            'evidence_strength': 'Tier I (FDA / NCCN Cat 1)' if tier1_targets else ('Tier II (Clinical Trials)' if tier2_targets else 'Tier III/IV')
        }


if __name__ == '__main__':
    print("=" * 70)
    print("MOLECULAR TUMOR BOARD (MTB) PRECISION ONCOLOGY ENGINE")
    print("=" * 70)
    for c in MolecularTumorBoardEngine.list_all_cases():
        print(f"\n[{c['case_id']}] {c['title']}")
        print(f"  Tanı: {c['primary_diagnosis']}")
        print(f"  Kılavuz Önerisi: {c['guideline']}")
