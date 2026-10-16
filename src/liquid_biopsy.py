"""
Liquid Biopsy & Circulating Tumor DNA (ctDNA) Minimal Residual Disease (MRD) Engine.
Monitors longitudinal Variant Allele Frequency (VAF %), molecular treatment response,
clonal dynamics, and acquired resistance mutation emergence from plasma cell-free DNA.
"""

from typing import Dict, List, Any, Optional

ACQUIRED_RESISTANCE_MUTATIONS: Dict[str, Dict[str, Any]] = {
    'EGFR_C797S': {
        'target_therapy': 'Osimertinib (3rd-gen TKI)',
        'gene': 'EGFR',
        'mutation': 'p.Cys797Ser (c.2389T>A)',
        'mechanism': 'Kovalent bağlanma cebindeki sistein 797 kalıntısının serin ile değişmesi; Osimertinibin kovalent bağ kurmasını engeller.',
        'actionable_bypass': 'Allosterik 4. jenerasyon EGFR TKI (örn. BLU-945, BDTX-1535) veya antikor-ilaç konjugatı (Patritumab Deruxtecan - HER3-DXd).'
    },
    'MET_AMPLIFICATION': {
        'target_therapy': 'Osimertinib / Gefitinib',
        'gene': 'MET',
        'mutation': 'MET Gene Copy Number Gain (Copy Number >= 5)',
        'mechanism': 'EGFR yolağını baypas ederek ERBB3 fosforilasyonunu ve downstream PI3K/AKT sinyalini devam ettirir.',
        'actionable_bypass': 'Kombinasyon: Osimertinib + MET inhibitörü Savolitinib (SAVANNAH Faz 2 çalışması) veya Tepotinib.'
    },
    'KRAS_Y96D': {
        'target_therapy': 'Sotorasib / Adagrasib (KRAS G12C İnhibitörleri)',
        'gene': 'KRAS',
        'mutation': 'p.Tyr96Asp (Switch-II cebi direnci)',
        'mechanism': 'Switch-II allosterik cebindeki Y96 mutasyonu, kovalent inhibitörlerin bağlanmasını sterik olarak bloke eder.',
        'actionable_bypass': 'Pan-KRAS inhibitörleri (RMC-6236, Darovasertib) veya SOS1/SHP2 inhibitör kombinasyonları.'
    },
    'BRCA_REVERSION': {
        'target_therapy': 'Olaparib / Talazoparib (PARP İnhibitörleri)',
        'gene': 'BRCA1 / BRCA2',
        'mutation': 'Secondary in-frame restore / Reversion mutation',
        'mechanism': 'İkincil mutasyon okuma çerçevesini yeniden düzelterek homolog rekombinasyon (HR) kapasitesini geri kazandırır.',
        'actionable_bypass': 'ATR inhibitörleri (Ceralasertib), WEE1 inhibitörleri (Adavosertib) veya Platin bazlı ikili kemoterapi.'
    }
}


class LiquidBiopsyEngine:
    """
    Evaluates plasma ctDNA liquid biopsy samples over time.
    Calculates molecular response, minimal residual disease (MRD) status,
    and identifies acquired resistance clones.
    """

    @classmethod
    def evaluate_longitudinal_ctdna(cls, baseline_vaf_percent: float,
                                    followup_vaf_percent: float,
                                    primary_gene: str = "EGFR",
                                    detected_new_mutations: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Calculates molecular response classification and relapse lead time
        based on longitudinal ctDNA VAF changes.
        """
        detected_new = detected_new_mutations or []

        # Relative change in VAF
        if baseline_vaf_percent > 0:
            delta_vaf_percent = round(((followup_vaf_percent - baseline_vaf_percent) / baseline_vaf_percent) * 100.0, 1)
        else:
            delta_vaf_percent = 100.0 if followup_vaf_percent > 0 else 0.0

        # Molecular Response Classification (RECIST molecular proxy)
        if followup_vaf_percent < 0.05:
            response_class = "Molecular Complete Response (mCR) - ctDNA Temizliği"
            response_code = "MCR"
            color = "#10b981"
            mrd_status = "MRD Negatif (Derin Moleküler Yanıt)"
            relapse_risk = "Düşük (Uzun PFS Beklentisi)"
            clinical_action = "Mevcut hedefe yönelik veya adjuvan tedaviye devam ediniz. 3 ay sonra rutin ctDNA takibi önerilir."
        elif delta_vaf_percent <= -50.0:
            response_class = "Molecular Partial Response (mPR)"
            response_code = "MPR"
            color = "#38bdf8"
            mrd_status = "MRD Düşük Pozitif (Tedaviye Yanıt Var)"
            relapse_risk = "Orta Düzey"
            clinical_action = "Tümör klonu belirgin şekilde baskılanmaktadır. Tedaviye devam ediniz."
        elif delta_vaf_percent <= 25.0:
            response_class = "Molecular Stable Disease (mSD)"
            response_code = "MSD"
            color = "#f59e0b"
            mrd_status = "MRD Pozitif (Persistan Klonalite)"
            relapse_risk = "Orta/Yüksek"
            clinical_action = "Klon stabil seyretmektedir. Kısa aralıklı (6-8 hafta) ctDNA ve radyolojik izlem önerilir."
        else:
            response_class = "Molecular Progressive Disease (mPD) - Moleküler Nüks"
            response_code = "MPD"
            color = "#ef4444"
            mrd_status = "MRD Yüksek Pozitif (Moleküler Nüks / Direnç)"
            relapse_risk = "Çok Yüksek (Radyolojik Nüksten 3-6 Ay Önce Saptandı)"
            clinical_action = "Moleküler progresyon saptandı! Radyolojik nüks beklenmeksizin direnç paneli taranmalı ve rejim değişikliği planlanmalıdır."

        # Check for resistance mutations
        resistance_findings = []
        for mut_key in detected_new:
            clean_key = mut_key.strip().upper().replace(" ", "_")
            if clean_key in ACQUIRED_RESISTANCE_MUTATIONS:
                info = ACQUIRED_RESISTANCE_MUTATIONS[clean_key]
                resistance_findings.append({
                    'mutation_key': clean_key,
                    'gene': info['gene'],
                    'mutation': info['mutation'],
                    'mechanism': info['mechanism'],
                    'actionable_bypass': info['actionable_bypass']
                })

        return {
            'target_gene': primary_gene.upper(),
            'baseline_vaf_percent': baseline_vaf_percent,
            'followup_vaf_percent': followup_vaf_percent,
            'vaf_change_percent': delta_vaf_percent,
            'delta_vaf_percent': delta_vaf_percent,
            'response_classification': response_class,
            'response_code': response_code,
            'color_code': color,
            'mrd_status': mrd_status,
            'relapse_risk': relapse_risk,
            'relapse_lead_time_months': '3 ila 6 Ay (Radyolojik Progresyon Öncesi Erken Uyarı)' if response_code == 'MPD' else 'N/A',
            'clinical_action': clinical_action,
            'acquired_resistance_detected': len(resistance_findings) > 0,
            'resistance_clones': resistance_findings,
            'resistance_findings': resistance_findings
        }


if __name__ == '__main__':
    print("=" * 70)
    print("LIQUID BIOPSY & ctDNA MINIMAL RESIDUAL DISEASE (MRD) ENGINE")
    print("=" * 70)

    # Scenario 1: Molecular Complete Response (ctDNA Clearance)
    res_clear = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(baseline_vaf_percent=6.5, followup_vaf_percent=0.02, primary_gene="EGFR")
    print(f"\n[Scenario 1: Post-Op ctDNA Clearance]")
    print(f"  VAF Change: {res_clear['baseline_vaf_percent']}% -> {res_clear['followup_vaf_percent']}% (Delta: {res_clear['vaf_change_percent']}%)")
    print(f"  Response: {res_clear['response_classification']} -> {res_clear['mrd_status']}")
    print(f"  Action: {res_clear['clinical_action']}")

    # Scenario 2: Molecular Progression with Acquired EGFR C797S
    res_prog = LiquidBiopsyEngine.evaluate_longitudinal_ctdna(
        baseline_vaf_percent=0.8,
        followup_vaf_percent=4.2,
        primary_gene="EGFR",
        detected_new_mutations=['EGFR_C797S', 'MET_AMPLIFICATION']
    )
    print(f"\n[Scenario 2: Molecular Relapse & Resistance Emergence]")
    print(f"  VAF Change: {res_prog['baseline_vaf_percent']}% -> {res_prog['followup_vaf_percent']}% (Delta: +{res_prog['vaf_change_percent']}%)")
    print(f"  Response: {res_prog['response_classification']} ({res_prog['relapse_risk']})")
    print(f"  Lead Time: {res_prog['relapse_lead_time_months']}")
    print(f"  Acquired Clones Detected ({len(res_prog['resistance_clones'])}):")
    for rc in res_prog['resistance_clones']:
        print(f"    * {rc['mutation']} [{rc['gene']}]: {rc['mechanism']}")
        print(f"      -> Actionable Bypass: {rc['actionable_bypass']}")
