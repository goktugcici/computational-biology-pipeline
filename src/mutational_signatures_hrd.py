"""
COSMIC Mutational Signatures & HRD (Homologous Recombination Deficiency) Engine.
Models 96-trinucleotide mutational contexts (SBS1, SBS2/13 APOBEC, SBS3 HRD, SBS4 Smoking,
SBS7 UV, SBS10 POLE, SBS31 Platinum), calculates genomic scarring indices
(HRD Score = LOH + TAI + LST), and predicts PARP inhibitor and platinum sensitivity.
"""

from typing import Dict, List, Any, Optional

CANONICAL_MUTATIONAL_SIGNATURES: Dict[str, Dict[str, Any]] = {
    'SBS1_DEAMINATION': {
        'id': 'SBS1_DEAMINATION',
        'title': 'SBS1: Spontan 5-Metilsitozin Deaminasyonu (Biyolojik Yaş Saati)',
        'etiology': 'Yaşlanmaya bağlı spontan sitozin deaminasyonu (CpG dinükleotitlerinde mC > T dönüşümü).',
        'mutation_type': 'C>T transisyonları (özellikle N[C>T]G bağlamında)',
        'associated_cancers': ['Tüm insan kanserleri (yaşla korele evrensel imza)', 'Kolorektal', 'Mide', 'Gliom'],
        'typical_prevalence': 'Neredeyse tüm kanser türlerinde arka plan saati olarak %10-40 oranında bulunur.',
        'clinical_actionability': 'Biyolojik saattir; tek başına hedefe yönelik tedavi hedefi oluşturmaz.',
        'relative_contributions': {'C>A': 5, 'C>G': 5, 'C>T': 75, 'T>A': 5, 'T>C': 5, 'T>G': 5}
    },
    'SBS2_13_APOBEC': {
        'id': 'SBS2_13_APOBEC',
        'title': 'SBS2 & SBS13: APOBEC Sitidin Deaminaz Enzim Hiperaktivasyonu',
        'etiology': 'APOBEC3A ve APOBEC3B sitidin deaminaz enzimlerinin aşırı aktivasyonu ve tek iplikli DNA hasarı (kataegis).',
        'mutation_type': 'C>T ve C>G mutasyonları (özellikle T[C>T]N ve T[C>G]A bağlamında)',
        'associated_cancers': ['Mesane (Ürotelyal Karsinom - %70+)', 'Meme Kanseri (özellikle HER2+)', 'Serviks Kanseri', 'Baş-Boyun'],
        'typical_prevalence': 'Ürotelyal karsinom ve meme kanserlerinde genomik mutasyonların %50-80\'ini oluşturabilir.',
        'clinical_actionability': 'Yüksek TMB oluşturur; Anti-PD-1/PD-L1 kontrol noktası blokajına (İmmünoterapi) artmış duyarlılık ile ilişkilidir.',
        'relative_contributions': {'C>A': 10, 'C>G': 35, 'C>T': 45, 'T>A': 4, 'T>C': 3, 'T>G': 3}
    },
    'SBS3_HRD_BRCA': {
        'id': 'SBS3_HRD_BRCA',
        'title': 'SBS3: Homolog Rekombinasyon Eksikliği (HRD / BRCA1/2 Kusuru)',
        'etiology': 'Çift iplik DNA kırıklarının hatasız homolog rekombinasyon (HRR) ile onarılamaması; mikrohomoloji aracılı uç birleşmesi (MMEJ) kullanımı.',
        'mutation_type': 'Tüm baz değişimleri boyunca geniş spektrumlu dağınık mutasyonlar ve >3 bp mikrohomoloji delesyonları.',
        'associated_cancers': ['Yüksek Dereceli Seröz Over Kanseri (%50)', 'Üçlü Negatif Meme Kanseri (%20)', 'Prostat (mCRPC)', 'Pankreas'],
        'typical_prevalence': 'BRCA1, BRCA2, PALB2, RAD51C alterasyonu olan tümörlerde hakim imzadır.',
        'clinical_actionability': 'PARP İnhibitörlerine (Olaparib, Rucaparib, Niraparib, Talazoparib) ve Platin Türevlerine (Sisplatin, Karboplatin) AŞIRI DUYARLILIK.',
        'relative_contributions': {'C>A': 20, 'C>G': 15, 'C>T': 25, 'T>A': 15, 'T>C': 15, 'T>G': 10}
    },
    'SBS4_SMOKING': {
        'id': 'SBS4_SMOKING',
        'title': 'SBS4: Tütün Dumanı Kanserojenleri (Benzo[a]piren Adüktleri)',
        'etiology': 'Sigara dumanındaki polisiklik aromatik hidrokarbonların (özellikle Benzo[a]piren) DNA guanin bazlarına kovalent bağlanması ve transkripsiyonel iplik sapması.',
        'mutation_type': 'Karakteristik C>A transversiyonları (özellikle C[C>A]A ve G[C>A]T bağlamında).',
        'associated_cancers': ['Akciğer Adenokarsinomu', 'Akciğer Skuamöz Hücreli Karsinomu', 'Larenks', 'Özofagus'],
        'typical_prevalence': 'Sigara içicisi akciğer kanserlerinde TMB\'nin %60-90\'ını oluşturur; paket-yıl ile güçlü doğrusal korelasyon gösterir.',
        'clinical_actionability': 'Yüksek TMB (>10 mut/Mb) sağlar; Anti-PD-1 monoterapisine (Pembrolizumab) mükemmel yanıt olasılığı oluşturur.',
        'relative_contributions': {'C>A': 65, 'C>G': 10, 'C>T': 12, 'T>A': 5, 'T>C': 5, 'T>G': 3}
    },
    'SBS7_UV_LIGHT': {
        'id': 'SBS7_UV_LIGHT',
        'title': 'SBS7a/b/c/d: Ultraviyole (UV) Güneş Işığı Hasarı (Pirimidin Dimerleri)',
        'etiology': 'Güneş ışığı UV radyasyonunun bitişik pirimidin bazları arasında siklobütan pirimidin dimerleri ve 6-4 fotopürünler oluşturması.',
        'mutation_type': 'Karakteristik C>T transisyonları ve çift nükleotit CC>TT dinükleotit değişimleri (özellikle dinükleotit pirimidinlerde).',
        'associated_cancers': ['Kutanöz Malign Melanom (%90+)', 'Skuamöz Hücreli Cilt Kanseri', 'Bazal Hücreli Karsinom'],
        'typical_prevalence': 'Güneşe maruz kalan cilt tümörlerinde yüzlerce mut/Mb ile dünyadaki en yüksek TMB değerlerini yaratır.',
        'clinical_actionability': 'Yüksek neoantijen yükü nedeniyle İmmün Kontrol Noktası Blokajına (Nivolumab + Ipilimumab) en duyarlı kanser grubudur.',
        'relative_contributions': {'C>A': 4, 'C>G': 2, 'C>T': 85, 'T>A': 3, 'T>C': 4, 'T>G': 2}
    },
    'SBS10_POLE_PROOFREADING': {
        'id': 'SBS10_POLE_PROOFREADING',
        'title': 'SBS10a/b: DNA Polimeraz Epsilon (POLE) Proofreading Ekzonükleaz Kusuru',
        'etiology': 'DNA replikasyon polimerazı POLE ekzonükleaz domaininde (örn. P286R, V411L) mutasyon sonucu replikasyon hata düzeltme kabiliyetinin kaybı.',
        'mutation_type': 'Devasa C>A ve T>G transversiyonları (özellikle T[C>A]T bağlamında ultra-hipermutatör fenotip).',
        'associated_cancers': ['Endometriyum Kanseri (%7-10)', 'Kolorektal Kanser (%1-2)'],
        'typical_prevalence': 'Genomda 100 ile 500+ mutasyon/Mb seviyesinde aşırı hipermutasyon (Ultra-Hypermutator).',
        'clinical_actionability': 'Mükemmel prognoz ve Anti-PD-1 immünoterapisine olağanüstü yüksek ve kalıcı yanıt oranı.',
        'relative_contributions': {'C>A': 48, 'C>G': 10, 'C>T': 22, 'T>A': 5, 'T>C': 5, 'T>G': 10}
    }
}


class MutationalSignaturesHrdEngine:
    """
    Bioinformatics engine for COSMIC Mutational Signatures deconvolution,
    Genomic Scar Index (LOH + TAI + LST = HRD Score), and PARP/Platinum sensitivity estimation.
    """

    @classmethod
    def list_signatures(cls) -> List[Dict[str, Any]]:
        return list(CANONICAL_MUTATIONAL_SIGNATURES.values())

    @classmethod
    def get_signature(cls, signature_id: str) -> Optional[Dict[str, Any]]:
        return CANONICAL_MUTATIONAL_SIGNATURES.get(signature_id)

    @classmethod
    def calculate_hrd_genomic_scar(
        cls,
        loh_score: int,
        tai_score: int,
        lst_score: int,
        brca_status: str = 'WT',
        cancer_type: str = 'Over Kanseri'
    ) -> Dict[str, Any]:
        """
        Calculates composite HRD score:
        HRD Score = LOH (Loss of Heterozygosity) + TAI (Telomeric Allelic Imbalance) + LST (Large-scale State Transitions).
        FDA threshold (Myriad myChoice CDx standard): HRD Score >= 42 or deleterious BRCA1/2 mutation indicates HRD-positive.
        """
        # Clamp inputs
        loh = max(0, min(40, int(loh_score)))
        tai = max(0, min(40, int(tai_score)))
        lst = max(0, min(40, int(lst_score)))

        composite_hrd = loh + tai + lst
        brca_pathogenic = brca_status.upper() in ['MUT', 'MUTANT', 'PATHOGENIC', 'BRCA1_MUT', 'BRCA2_MUT']

        # Clinical determination
        hrd_positive = (composite_hrd >= 42) or brca_pathogenic

        if brca_pathogenic and composite_hrd >= 42:
            clinical_status = "HRD Pozitif (BRCA Mutant + Yüksek Genomik Skar)"
            confidence = "Çok Yüksek (%95 Güven)"
            parp_recommendation = "Olaparib, Niraparib, Rucaparib veya Talazoparib Kategori 1 Endikasyon"
            platinum_recommendation = "Platin bazlı kemoterapiye (Karboplatin / Sisplatin) belirgin üstün yanıt beklenir."
        elif brca_pathogenic:
            clinical_status = "HRD Pozitif (BRCA Mutant / Fonksiyon Kaybı)"
            confidence = "Yüksek (Moleküler Tanımlı)"
            parp_recommendation = "PARP inhibitörleri ile idame veya monoterapi Kategori 1 önerilir."
            platinum_recommendation = "Platin kemoterapisine yüksek duyarlılık beklenir."
        elif composite_hrd >= 42:
            clinical_status = "HRD Pozitif (BRCA Wild-Type / Epigenetik veya 'BRCAness' Genomik Skar)"
            confidence = "Yüksek (Genomik Skar Eşiği Aşılmış: ≥42)"
            parp_recommendation = "Niraparib (PRIMA / NOVA) veya Olaparib + Bevacizumab (PAOLA-1) idamesi uygundur."
            platinum_recommendation = "Platin kemoterapisine iyi yanıt olasılığı."
        else:
            clinical_status = "HRD Negatif (Homolog Onarım Yeterli / HR-Proficient)"
            confidence = "Normal (Düşük Skar: <42)"
            parp_recommendation = "PARP inhibitörü monoterapisinden beklenen klinik fayda sınırlıdır; alternatif biyobelirteçler araştırılmalıdır."
            platinum_recommendation = "Standart platin duyarlılığı (özel aşırı duyarlılık beklenmez)."

        return {
            'cancer_type': cancer_type,
            'loh_score': loh,
            'tai_score': tai,
            'lst_score': lst,
            'composite_hrd_score': composite_hrd,
            'hrd_positive': hrd_positive,
            'brca_status': brca_status,
            'clinical_status': clinical_status,
            'confidence': confidence,
            'parp_recommendation': parp_recommendation,
            'platinum_recommendation': platinum_recommendation,
            'resistance_warning': "Direnç riski: Tedavi sürecinde BRCA1/2 geninde okuma çerçevesini onaran ikincil geri dönüş (secondary reversion) mutasyonları gelişebilir."
        }

    @classmethod
    def decompose_patient_spectrum(
        cls,
        c_to_a: float,
        c_to_g: float,
        c_to_t: float,
        t_to_a: float,
        t_to_c: float,
        t_to_g: float
    ) -> Dict[str, Any]:
        """
        Deconvolutes and matches single-base substitution (SBS) spectrum to COSMIC etiologies.
        """
        total = c_to_a + c_to_g + c_to_t + t_to_a + t_to_c + t_to_g
        if total <= 0:
            total = 100.0

        p_ca = (c_to_a / total) * 100.0
        p_cg = (c_to_g / total) * 100.0
        p_ct = (c_to_t / total) * 100.0
        p_ta = (t_to_a / total) * 100.0
        p_tc = (t_to_c / total) * 100.0
        p_tg = (t_to_g / total) * 100.0

        # Heuristic best signature match
        if p_ca >= 45.0:
            dominant_sig = 'SBS4_SMOKING'
            reason = "Aşırı C>A transversiyonları tütün dumanı polisiklik aromatik hidrokarbonları (benzo[a]piren) ile uyumludur."
        elif p_ct >= 70.0:
            dominant_sig = 'SBS7_UV_LIGHT'
            reason = "Baskın C>T pirimidin dimer transisyonları güneş / ultraviyole (UV) radyasyon hasarını düşündürmektedir."
        elif p_cg >= 25.0 and p_ct >= 35.0:
            dominant_sig = 'SBS2_13_APOBEC'
            reason = "Eşzamanlı C>G ve C>T artışı APOBEC3 sitidin deaminaz atağını (kataegis) yansıtmaktadır."
        elif p_ca >= 30.0 and p_tg >= 8.0:
            dominant_sig = 'SBS10_POLE_PROOFREADING'
            reason = "Yüksek C>A ve T>G oranı DNA polimeraz epsilon (POLE) proofreading ekzonükleaz kusuruna işaret eder."
        elif p_ct >= 55.0:
            dominant_sig = 'SBS1_DEAMINATION'
            reason = "CpG dinükleotitlerinde 5-metilsitozin deaminasyonuna bağlı biyolojik yaşlanma imzasıdır."
        else:
            dominant_sig = 'SBS3_HRD_BRCA'
            reason = "Tüm baz değişimlerine yayılan geniş ve dengeli mutasyon dağılımı bozulmuş homolog onarım (HRD) ile uyumludur."

        sig_data = CANONICAL_MUTATIONAL_SIGNATURES.get(dominant_sig, {})

        return {
            'spectrum_percentages': {
                'C>A': round(p_ca, 1),
                'C>G': round(p_cg, 1),
                'C>T': round(p_ct, 1),
                'T>A': round(p_ta, 1),
                'T>C': round(p_tc, 1),
                'T>G': round(p_tg, 1)
            },
            'dominant_signature': dominant_sig,
            'dominant_title': sig_data.get('title', ''),
            'dominant_etiology': sig_data.get('etiology', ''),
            'match_reason': reason
        }
