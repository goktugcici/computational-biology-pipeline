"""
Tumor Mutational Burden (TMB) and Microsatellite Instability (MSI)
Immunotherapy Response Prediction Engine for Precision Oncology.
"""

from typing import Dict, List, Any, Optional

CANONICAL_TUMOR_COHORTS: Dict[str, Dict[str, Any]] = {
    'MELANOMA_CUTANEOUS': {
        'name': 'Cutaneous Melanoma (UV-induced)',
        'typical_tmb': 18.4,
        'msi_status': 'MSS',
        'pdl1_tps': 45,
        'primary_regimen': 'Nivolumab + Ipilimumab (Anti-PD-1 + Anti-CTLA-4)',
        'description': 'UV ışınımı kaynaklı yüksek mutasyon yükü ve bol neoantijen üretimi.'
    },
    'LUNG_NSCLC_SMOKER': {
        'name': 'Non-Small Cell Lung Cancer (NSCLC)',
        'typical_tmb': 12.8,
        'msi_status': 'MSS',
        'pdl1_tps': 60,
        'primary_regimen': 'Pembrolizumab Monoterapisi (Keynote-042)',
        'description': 'Tütün karsinojenlerine bağlı transversiyon mutasyonları ve TMB-High.'
    },
    'COLORECTAL_MSI_H': {
        'name': 'Colorectal Carcinoma (Lynch / dMMR)',
        'typical_tmb': 28.5,
        'msi_status': 'MSI-H (dMMR)',
        'pdl1_tps': 30,
        'primary_regimen': 'Pembrolizumab veya Dostarlimab (Doku-Agnostik Onay)',
        'description': 'MLH1/MSH2 defekti kaynaklı mikrosatellit instabilitesi ve hipermutasyon.'
    },
    'PANCREATIC_MSS': {
        'name': 'Pancreatic Ductal Adenocarcinoma',
        'typical_tmb': 2.1,
        'msi_status': 'MSS (pMMR)',
        'pdl1_tps': 0,
        'primary_regimen': 'Kemoterapi (FOLFIRINOX) - İmmünoterapi Dirençli',
        'description': 'İmmünolojik olarak "soğuk" tümör mikroçevresi ve düşük TMB.'
    }
}


class ImmunotherapyPredictor:
    """
    Evaluates patient molecular biomarkers (TMB, MSI, PD-L1)
    to predict therapeutic efficacy of immune checkpoint blockade.
    """

    @classmethod
    def evaluate_response(cls, tmb_mut_per_mb: float, msi_status: str = 'MSS', pdl1_tps_percent: float = 0.0) -> Dict[str, Any]:
        """
        Calculates composite Immunotherapy Efficacy Score (0 to 100%)
        based on FDA criteria (Keynote-158 TMB >= 10 mut/Mb, MSI-H tissue-agnostic approval).
        """
        # 1. TMB Score (0 to 45 pts)
        is_tmb_high = tmb_mut_per_mb >= 10.0
        if tmb_mut_per_mb >= 20.0:
            tmb_score = 45.0
        elif tmb_mut_per_mb >= 10.0:
            tmb_score = 35.0 + (tmb_mut_per_mb - 10.0) * 1.0
        else:
            tmb_score = (tmb_mut_per_mb / 10.0) * 20.0

        # 2. MSI Score (0 to 35 pts)
        is_msi_high = 'MSI-H' in msi_status.upper() or 'DMMR' in msi_status.upper()
        msi_score = 35.0 if is_msi_high else 5.0

        # 3. PD-L1 Score (0 to 20 pts)
        pdl1_clamped = min(100.0, max(0.0, pdl1_tps_percent))
        pdl1_score = (pdl1_clamped / 100.0) * 20.0

        # Total Composite Probability
        total_prob = round(tmb_score + msi_score + pdl1_score, 1)
        total_prob = min(98.0, max(5.0, total_prob))

        # Clinical Classification
        if total_prob >= 70.0:
            classification = 'Exceptional Response Likely (Yüksek Yanıt Beklentisi)'
            badge_color = '#10b981' # Emerald
            recommended_regimen = 'Pembrolizumab (Keytruda) veya Nivolumab Monoterapisi'
        elif total_prob >= 45.0:
            classification = 'Moderate Response Potential (Ilımlı Yanıt Potansiyeli)'
            badge_color = '#38bdf8' # Sky blue
            recommended_regimen = 'İkili İmmün Kontrol Blokajı (Nivolumab + Ipilimumab) veya Kemo-İmmünoterapi'
        else:
            classification = 'Immune Desert / Resistant ("Soğuk" Tümör / İmmünoterapi Dirençli)'
            badge_color = '#ef4444' # Red
            recommended_regimen = 'Standart Kemoterapi / Hedefe Yönelik Tedavi / ADC (Antikor-İlaç Konjugatı)'

        return {
            'tmb_mut_per_mb': round(tmb_mut_per_mb, 1),
            'is_tmb_high': is_tmb_high,
            'msi_status': msi_status.upper(),
            'is_msi_high': is_msi_high,
            'pdl1_tps_percent': round(pdl1_clamped, 1),
            'response_probability_percent': total_prob,
            'classification': classification,
            'badge_color': badge_color,
            'recommended_regimen': recommended_regimen,
            'biomarker_breakdown': {
                'tmb_contribution': round(tmb_score, 1),
                'msi_contribution': round(msi_score, 1),
                'pdl1_contribution': round(pdl1_score, 1)
            }
        }

    @classmethod
    def evaluate_cohort(cls, cohort_key: str) -> Dict[str, Any]:
        """Evaluates a canonical cancer cohort."""
        if cohort_key not in CANONICAL_TUMOR_COHORTS:
            cohort_key = 'LUNG_NSCLC_SMOKER'
        c = CANONICAL_TUMOR_COHORTS[cohort_key]
        res = cls.evaluate_response(c['typical_tmb'], c['msi_status'], c['pdl1_tps'])
        res['cohort_name'] = c['name']
        res['cohort_description'] = c['description']
        return res


if __name__ == '__main__':
    predictor = ImmunotherapyPredictor()
    eval_res = predictor.evaluate_cohort('COLORECTAL_MSI_H')
    print(f"Cohort: {eval_res['cohort_name']}")
    print(f"TMB: {eval_res['tmb_mut_per_mb']} mut/Mb | MSI: {eval_res['msi_status']} | PD-L1: {eval_res['pdl1_tps_percent']}%")
    print(f"Response Probability: {eval_res['response_probability_percent']}% -> {eval_res['classification']}")
    print(f"Recommended Regimen: {eval_res['recommended_regimen']}")
