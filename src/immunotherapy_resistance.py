"""
Immunotherapy Resistance & Immune Microenvironment (TIME) Predictor Engine.
Models primary, adaptive, and acquired resistance to immune checkpoint blockade (ICB),
including HLA-A/B/C loss, B2M loss, STK11/KEAP1 co-mutations, PTEN loss,
JAK1/JAK2 loss-of-function, MDSC / Treg infiltration, and VEGF/TGF-beta immunosuppression.
"""

from typing import Dict, List, Any, Optional

CANONICAL_ICB_SCENARIOS: Dict[str, Dict[str, Any]] = {
    'NSCLC_STK11_KEAP1': {
        'id': 'NSCLC_STK11_KEAP1',
        'cancer_type': 'Küçük Hücreli Dışı Akciğer Kanseri (NSCLC)',
        'title': 'KRAS Mutasyonu + STK11/KEAP1 Eş-Mutasyonu (Soğuk İmmün Çevre)',
        'demographics': {'age': 64, 'gender': 'Erkek', 'smoking': '35 Paket-Yıl'},
        'biomarkers': {'tmb': 14.2, 'msi': 'MSS', 'pdl1_tps': 1.0, 'hla_loss': False},
        'mutations': ['KRAS p.G12C', 'STK11 (LKB1) Truncating', 'KEAP1 Nonsense'],
        'time_phenotype': 'İmmün Çöl (Immune Desert / Cold Tumor)',
        'mechanism': 'STK11/LKB1 inaktivasyonu STING yolağını susturur; T-hücre infiltrasyonunu ve tip-I interferon salınımını bloke eder. TMB yüksek olsa dahi Anti-PD-1 monoterapisine dirençlidir.',
        'clinical_strategy': 'Anti-PD-1/PD-L1 monoterapisi yerine Platin Bazlı Kemoterapi + Anti-CTLA-4 (Tremelimumab/Ipilimumab) + Anti-PD-1 kombinasyonu (CheckMate-9LA veya POSEIDON protokolü) veya KRAS G12C kovalent inhibitörü.',
        'evidence_level': 'NCCN Kategori 1 / ESMO ESCAT Tier I-A'
    },
    'MELANOMA_JAK1_B2M': {
        'id': 'MELANOMA_JAK1_B2M',
        'cancer_type': 'Malign Melanom (Kutanöz)',
        'title': 'Kazanılmış İmmün Kaçış (B2M Delesyonu & JAK1 İnaktivasyonu)',
        'demographics': {'age': 55, 'gender': 'Kadın', 'smoking': '0 Paket-Yıl'},
        'biomarkers': {'tmb': 28.5, 'msi': 'MSS', 'pdl1_tps': 80.0, 'hla_loss': True},
        'mutations': ['BRAF p.V600E', 'B2M (Beta-2 Mikroglobulin) Exon 1 Delesyonu', 'JAK1 p.Q503*'],
        'time_phenotype': 'Kazanılmış Antijen Sunum Yetmezliği (Acquired Resistance)',
        'mechanism': 'B2M kaybı yüzeyde MHC-I antijen sunumunu tamamen ortadan kaldırır. CD8+ T-hücreleri tümör neoantijenlerini tanıyamaz. JAK1 kaybı ise interferon-gama reseptör sinyalini ve PD-L1 indüksiyonunu felç eder.',
        'clinical_strategy': 'T-hücre bağımsız yollar: NK hücre aktivatörleri, sitotoksik bispesifik veya hedefe yönelik BRAF+MEK inhibitörleri (Dabrafenib+Trametinib) veya TIL (Tümör İnfiltre Eden Lenfosit) hücresel tedavisi.',
        'evidence_level': 'ESMO Tier I-B'
    },
    'COLON_DMMR_MSI_HIGH': {
        'id': 'COLON_DMMR_MSI_HIGH',
        'cancer_type': 'Kolorektal Adenokarsinom',
        'title': 'Sıcak Tümör & Güçlü İmmünojenik Yanıt (MSI-High / dMMR)',
        'demographics': {'age': 51, 'gender': 'Erkek', 'smoking': '5 Paket-Yıl'},
        'biomarkers': {'tmb': 42.0, 'msi': 'MSI-H', 'pdl1_tps': 65.0, 'hla_loss': False},
        'mutations': ['MLH1 Promoter Hipermetilasyonu', 'BRAF p.V600E', 'TGFBR2 Frameshift'],
        'time_phenotype': 'İmmün İltihaplı (Immune Inflamed / Hot Tumor)',
        'mechanism': 'Uyumsuzluk onarım kusuru (dMMR) binlerce mikrosatellit çerçeve kayması ve devasa neoantijen yükü üretir. Tümör mikroçevresi CD8+ ve Th1 hücreleriyle doludur.',
        'clinical_strategy': 'Birinci basamakta Pembrolizumab monoterapisi (KEYNOTE-177) veya Nivolumab + Ipilimumab (CheckMate-142). Cerrahi öncesi neoadjuvan kullanımda %60-70 tam patolojik yanıt (pCR) bildirilmektedir (NICHE-2 protokolü).',
        'evidence_level': 'NCCN Kategori 1 / FDA Birinci Basamak Onaylı'
    },
    'BREAST_TNBC_VEGF_MDSC': {
        'id': 'BREAST_TNBC_VEGF_MDSC',
        'cancer_type': 'Üçlü Negatif Meme Kanseri (TNBC)',
        'title': 'İmmün Dışlanmış (Immune Excluded) - Yüksek VEGF & TGF-β / MDSC Stroma',
        'demographics': {'age': 46, 'gender': 'Kadın', 'smoking': '0 Paket-Yıl'},
        'biomarkers': {'tmb': 6.5, 'msi': 'MSS', 'pdl1_tps': 15.0, 'hla_loss': False},
        'mutations': ['TP53 p.R248Q', 'PIK3CA p.H1047R', 'PTEN Delesyon'],
        'time_phenotype': 'İmmün Dışlanmış (Immune Excluded / Stroma Zengin)',
        'mechanism': 'T hücreleri stromada sıkışmıştır; tümör parankimine giremez. Aşırı VEGF ve TGF-beta ekspresyonu miyeloid kökenli baskılayıcı hücreleri (MDSC) ve Treg hücrelerini çekerek fiziko-kimyasal bariyer oluşturur.',
        'clinical_strategy': 'Anti-VEGF/Anti-Anjiyogenez veya Kemoterapi (Paklitaksel/Karboplatin) ile stroma gevşetildikten sonra Anti-PD-1 (Pembrolizumab - KEYNOTE-522 veya KEYNOTE-355 protokolü).',
        'evidence_level': 'NCCN Kategori 1'
    }
}


class ImmunotherapyResistanceEngine:
    """
    Evaluates Tumor Immune Microenvironment (TIME), HLA loss,
    STK11/KEAP1/B2M driven primary/acquired resistance, and predicts ICB efficacy.
    """

    @classmethod
    def list_scenarios(cls) -> List[Dict[str, Any]]:
        return list(CANONICAL_ICB_SCENARIOS.values())

    @classmethod
    def get_scenario(cls, scenario_id: str) -> Optional[Dict[str, Any]]:
        return CANONICAL_ICB_SCENARIOS.get(scenario_id)

    @classmethod
    def evaluate_tumor_immune_status(
        cls,
        cancer_type: str,
        tmb: float,
        msi_status: str,
        pdl1_percent: float,
        stk11_mut: bool = False,
        keap1_mut: bool = False,
        b2m_loss: bool = False,
        jak_mut: bool = False,
        pten_loss: bool = False,
        vegf_high: bool = False
    ) -> Dict[str, Any]:
        """
        Calculates composite ICB Sensitivity Score (0-100), resistance mechanisms,
        and provides tailored biomarker-guided therapeutic guidance.
        """
        # Baseline probability from classic biomarkers
        icb_score = 30.0

        if msi_status.upper() in ['MSI-H', 'MSI-HIGH', 'DMMR']:
            icb_score += 45.0
        elif tmb >= 20.0:
            icb_score += 35.0
        elif tmb >= 10.0:
            icb_score += 20.0

        if pdl1_percent >= 50.0:
            icb_score += 25.0
        elif pdl1_percent >= 1.0:
            icb_score += 15.0

        resistance_reasons: List[str] = []
        bypass_therapies: List[str] = []

        # Detrimental negative drivers
        if stk11_mut:
            icb_score -= 30.0
            resistance_reasons.append("STK11 (LKB1) inaktivasyonu STING yolağını baskılar; 'soğuk tümör' fenotipi oluşturur.")
            bypass_therapies.append("Anti-CTLA-4 (Ipilimumab/Tremelimumab) + Kemoterapi kombinasyonu")

        if keap1_mut:
            icb_score -= 20.0
            resistance_reasons.append("KEAP1 mutasyonu NRF2 aşırı aktivasyonu ve immün direnci indükler.")

        if b2m_loss:
            icb_score -= 45.0
            resistance_reasons.append("B2M kaybı / delesyonu: Hücre yüzeyinde MHC-I antijen sunumu çökmüştür. CD8+ T-hücreleri tümörü göremez.")
            bypass_therapies.append("NK-hücre aktivatörleri, hedefe yönelik tirozin kinaz inhibitörleri veya T-hücre bağımsız sitotoksisite")

        if jak_mut:
            icb_score -= 35.0
            resistance_reasons.append("JAK1/JAK2 inaktivasyonu: İnterferon-gama sinyali yanıtsızdır; PD-L1 yukarı regüle olamaz.")

        if pten_loss:
            icb_score -= 15.0
            resistance_reasons.append("PTEN kaybı immünosüpresif sitokin salınımını tetikler ve T-hücre infiltrasyonunu kısıtlar.")

        if vegf_high:
            icb_score -= 15.0
            resistance_reasons.append("Yüksek VEGF & TGF-beta: İmmün dışlanmış (immune-excluded) yoğun stroma ve MDSC birikimi.")
            bypass_therapies.append("Anti-VEGF (Bevacizumab/Lenvatinib) ile damarsal ve stromal normalizasyon")

        # Bound score between 5.0 and 99.0
        final_score = max(5.0, min(99.0, icb_score))

        if final_score >= 70.0:
            phenotype = "İmmün İltihaplı (Hot / Inflamed Tumor)"
            response_category = "YÜKSEK DUYARLILIK (Güçlü ICB Yanıtı Beklenir)"
            primary_rec = "Anti-PD-1 / Anti-PD-L1 Monoterapisi veya Standart İmmünoterapi Kombinasyonu"
        elif final_score >= 40.0:
            phenotype = "İmmün Dışlanmış / Orta Derece Reaktif (Immune Excluded)"
            response_category = "ORTA DÜZEY YANIT (Kombinasyon Tedavisi Gerekli)"
            primary_rec = "İmmünoterapi + Kemoterapi veya Anti-VEGF / Hedefli Ajan İkili Protokolü"
        else:
            phenotype = "İmmün Çöl / Dirençli Soğuk Tümör (Cold Tumor / Resistant)"
            response_category = "DÜŞÜK YANIT / PRİMER DİRENÇLİ"
            primary_rec = "Tek ajan Anti-PD-1 önerilmez; STING agonistleri, ikili kontrol noktası veya hedefe yönelik baypas"

        return {
            'cancer_type': cancer_type,
            'icb_sensitivity_score': round(final_score, 1),
            'tumor_phenotype': phenotype,
            'response_category': response_category,
            'resistance_reasons': resistance_reasons,
            'recommended_strategy': primary_rec,
            'bypass_therapies': bypass_therapies,
            'tmb': tmb,
            'msi_status': msi_status,
            'pdl1_percent': pdl1_percent
        }
