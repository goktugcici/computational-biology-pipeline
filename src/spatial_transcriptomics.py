"""
Spatial Transcriptomics & Single-Cell Tumor Architecture Engine.
Models 10x Genomics Visium / Xenium / CosMx spatial tissue neighborhoods,
Tertiary Lymphoid Structures (TLS), Cancer-Associated Fibroblast (CAF) desmoplastic barriers,
Tumor-Associated Macrophage (TAM) polarization (M1 vs. M2), and ligand-receptor spatial interactions.
"""

from typing import Dict, List, Any, Optional

CANONICAL_SPATIAL_SCENARIOS: Dict[str, Dict[str, Any]] = {
    'NSCLC_TLS_MATURE': {
        'id': 'NSCLC_TLS_MATURE',
        'cancer_type': 'Küçük Hücreli Dışı Akciğer Kanseri (NSCLC)',
        'title': 'Olgun Tersiyer Lenfoid Yapılar (Mature TLS) & Yüksek İmmünojenite',
        'tissue_architecture': 'Tümör periferinde organize germinal merkezli B ve T hücre yuvaları; CXCL13+ ve CD20+ zengin.',
        'cellular_composition': {
            'tumor_epithelial': 40,
            'cd8_t_cells': 25,
            'b_cells_plasma': 20,
            'caf_fibroblasts': 10,
            'tam_macrophages': 5
        },
        'tls_maturation_stage': 'Evre 3: Olgun Germinal Merkezli Sekonder TLS',
        'tls_score': 88.0,
        'stromal_barrier_score': 22.0,
        'ligand_receptor_pairs': ['CXCL13 - CXCR5', 'CCL19/21 - CCR7', 'CD40 - CD40L', 'PD-L1 - PD-1'],
        'clinical_implication': 'İmmün kontrol noktası blokajına (Anti-PD-1/PD-L1) olağanüstü yüksek ve kalıcı yanıt; uzun genel sağkalım ile korele.',
        'guidance': 'Birinci basamakta Pembrolizumab veya Kemo-İmmünoterapi ile derin patolojik remisyon (pCR) hedeflenmelidir.'
    },
    'PDAC_DESMOPLASTIC_CAF': {
        'id': 'PDAC_DESMOPLASTIC_CAF',
        'cancer_type': 'Pankreatik Duktal Adenokarsinom (PDAC)',
        'title': 'Desmoplastik Fibröz Stroma & İmmün Dışlanmış CAF Bariyeri',
        'tissue_architecture': 'Yoğun kollajen I/III, alfa-SMA+ miyofibroblastik CAF tabakası; T hücreleri stroma çeperinde hapsedilmiş.',
        'cellular_composition': {
            'tumor_epithelial': 25,
            'cd8_t_cells': 5,
            'b_cells_plasma': 3,
            'caf_fibroblasts': 52,
            'tam_macrophages': 15
        },
        'tls_maturation_stage': 'Evre 0: TLS Yok (İmmün Boşluk / Fibröz Çöl)',
        'tls_score': 8.0,
        'stromal_barrier_score': 92.0,
        'ligand_receptor_pairs': ['TGFB1 - TGFBR2', 'CXCL12 - CXCR4', 'FAP - Integrin', 'PDGFRA - PDGFB'],
        'clinical_implication': 'Aşırı stromal interstisyel sıvı basıncı ve hipoksi; standart kemoterapi ve kontrol noktası inhibitörlerine primer direnç.',
        'guidance': 'FOLFIRINOX veya Gemstabin+Nab-Paklitaksel öncesi stroma gevşetici stratejiler veya FAP hedefli faz çalışmaları.'
    },
    'GBM_HYPOXIC_TAM': {
        'id': 'GBM_HYPOXIC_TAM',
        'cancer_type': 'Glioblastoma Multiforme (GBM - IDH-Wildtype)',
        'title': 'Hipoksik Nekroz & M2 Polarize İmmünsüpresif Makrofaj Hakimiyeti',
        'tissue_architecture': 'Psödo-palizatlı nekroz sınırında yüksek HIF-1alfa, yoğun VEGF salınımı ve CD206+/CD163+ M2 TAM infiltrasyonu.',
        'cellular_composition': {
            'tumor_epithelial': 55,
            'cd8_t_cells': 4,
            'b_cells_plasma': 1,
            'caf_fibroblasts': 5,
            'tam_macrophages': 35
        },
        'tls_maturation_stage': 'Evre 0: TLS Yok (İmmünosüpresif Mikroçevre)',
        'tls_score': 5.0,
        'stromal_barrier_score': 74.0,
        'ligand_receptor_pairs': ['VEGFA - VEGFR2', 'CSF1 - CSF1R', 'CD47 - SIRPA', 'ANGPT2 - TIE2'],
        'clinical_implication': 'Makrofajlar tümör büyümesini kamçılar ve CD47 "don\'t eat me" sinyali ile fagositozu kilitler.',
        'guidance': 'Standart Stupp protokolü (Radyoterapi + Temozolomid); ilerleyen olgularda Anti-VEGF (Bevacizumab) veya CSF-1R inhibitörleri.'
    },
    'MELANOMA_INVASIVE_FRONT': {
        'id': 'MELANOMA_INVASIVE_FRONT',
        'cancer_type': 'Kutanöz Malign Melanom',
        'title': 'İnvaziv Tümör Sınırında CD8+ T Hücre İnfiltrasyonu & Adaptif PD-L1',
        'tissue_architecture': 'Tümör-stroma temas yüzeyinde yoğun sitotoksik T hücre infiltrasyonu, adaptif IFN-gama salınımı ve PD-L1 ko-lokalizasyonu.',
        'cellular_composition': {
            'tumor_epithelial': 45,
            'cd8_t_cells': 32,
            'b_cells_plasma': 8,
            'caf_fibroblasts': 8,
            'tam_macrophages': 7
        },
        'tls_maturation_stage': 'Evre 2: Primer Foliküler TLS (Organize Agregatlar)',
        'tls_score': 64.0,
        'stromal_barrier_score': 30.0,
        'ligand_receptor_pairs': ['IFNG - IFNGR1', 'PD-L1 - PD-1', 'LAG-3 - MHC-II', 'TIGIT - PVR'],
        'clinical_implication': 'Tümör hücresi immün baskı altındadır; PD-1 ve LAG-3 blokajına yüksek duyarlılık gösterir.',
        'guidance': 'Nivolumab + Relatlimab (Anti-LAG-3) veya Nivolumab + Ipilimumab ikili immünoterapisi Kategori 1 seçenektir.'
    }
}


class SpatialTranscriptomicsEngine:
    """
    Simulates spatial microarchitecture, cellular deconvolution,
    TLS maturation index, and desmoplastic stroma barriers.
    """

    @classmethod
    def list_scenarios(cls) -> List[Dict[str, Any]]:
        return list(CANONICAL_SPATIAL_SCENARIOS.values())

    @classmethod
    def get_scenario(cls, scenario_id: str) -> Optional[Dict[str, Any]]:
        return CANONICAL_SPATIAL_SCENARIOS.get(scenario_id)

    @classmethod
    def calculate_spatial_metrics(
        cls,
        cd8_density: float,
        b_cell_density: float,
        caf_density: float,
        m2_tam_ratio: float,
        cxcl13_positive: bool = False
    ) -> Dict[str, Any]:
        """
        Calculates TLS Maturation Score, Stroma Desmoplasia Barrier,
        and provides spatial translational oncologic recommendations.
        """
        # TLS Score calculation (0 - 100)
        base_tls = (cd8_density * 1.2) + (b_cell_density * 2.0)
        if cxcl13_positive:
            base_tls += 30.0
        # M2 TAM and CAF dampens TLS formation
        base_tls -= (caf_density * 0.3) + (m2_tam_ratio * 0.4)
        tls_score = max(0.0, min(100.0, base_tls))

        # Stroma Barrier calculation (0 - 100)
        stroma_barrier = (caf_density * 1.4) + (m2_tam_ratio * 0.6) - (cd8_density * 0.2)
        stroma_score = max(0.0, min(100.0, stroma_barrier))

        # Classification
        if tls_score >= 70.0:
            tls_stage = "Evre 3: Olgun Germinal Merkezli Sekonder TLS (Yüksek İmmün Yanıt)"
            icb_prediction = "Çok Yüksek İmmünoterapi Yanıtı (ORR >%65)"
        elif tls_score >= 40.0:
            tls_stage = "Evre 2: Erken Foliküler Agregat TLS (Orta İmmün Yanıt)"
            icb_prediction = "Orta Düzey İmmünoterapi Yanıtı (ORR ~%35-45)"
        elif tls_score >= 20.0:
            tls_stage = "Evre 1: Olgunlaşmamış Lenfoid Küme (Immature Cluster)"
            icb_prediction = "Düşük/Gecikmiş İmmün Yanıt"
        else:
            tls_stage = "Evre 0: TLS Yok (İmmün Çöl / Dışlanmış Doku)"
            icb_prediction = "Kontrol Noktası Monoterapisine Primer Dirençli"

        if stroma_score >= 70.0:
            stroma_status = "Aşırı Desmoplazi & Yoğun Fibrotik Bariyer (T Hücreleri Dışlanmış)"
            recommendation = "Anti-TGF-beta, FAP hedefli ajanlar veya stroma geçirgenliğini artıran sitotoksik kemoterapi kombinasyonu gereklidir."
        elif stroma_score >= 40.0:
            stroma_status = "Orta Düzey Stromal Yoğunluk"
            recommendation = "Standart kombinasyon rejimleri ve invaziv sınır değerlendirmesi."
        else:
            stroma_status = "Düşük Stromal Direnç (Serbest Hücresel İnfiltrasyon)"
            recommendation = "Doğrudan tümör parankimine erişim mevcuttur; hedefe yönelik ajanlar ve immünoterapi en yüksek penetrasyona ulaşır."

        return {
            'tls_score': round(tls_score, 1),
            'tls_stage': tls_stage,
            'icb_prediction': icb_prediction,
            'stroma_barrier_score': round(stroma_score, 1),
            'stroma_status': stroma_status,
            'clinical_recommendation': recommendation
        }
