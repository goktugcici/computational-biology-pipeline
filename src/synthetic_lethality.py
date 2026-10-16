"""
Synthetic Lethality & CRISPR DepMap Cancer Dependency Engine.
Models genetic vulnerabilities, synthetic lethal gene pairs, and precision targeting
derived from Broad Institute Cancer Dependency Map (DepMap Chronos / Project Score).
"""

from typing import Dict, List, Any, Optional

SYNTHETIC_LETHAL_CATALOG: Dict[str, Dict[str, Any]] = {
    'BRCA1_DEFICIENCY': {
        'biomarker_gene': 'BRCA1',
        'biomarker_state': 'Loss-of-Function / Pathogenic Inactivating Mutation (HRD+)',
        'synthetic_lethal_target': 'PARP1',
        'target_full_name': 'Poly(ADP-Ribose) Polymerase 1',
        'target_uniprot': 'P09874',
        'depmap_chronos_score': -1.24,
        'depmap_essentiality': 'Extreme Dependency in HRD+ Tumors',
        'vulnerability_score': 98.0,
        'mechanistic_rationale': 'Homolog rekombinasyon (HR) onarım defekti olan hücreler, tek zincir kırıklarının onarımında PARP aracılı baz eksizyon tamirine (BER) bağımlı hale gelir. PARP inhibisyonu çift zincir kırıklarına ve replikasyon çatalı çökmesine yol açarak tümör hücresini seçici olarak öldürür.',
        'clinical_therapies': [
            {'drug': 'Olaparib (Lynparza)', 'fda_status': 'FDA Onaylı', 'indication': 'Over, Meme, Prostat, Pankreas CA'},
            {'drug': 'Talazoparib (Talzenna)', 'fda_status': 'FDA Onaylı', 'indication': 'gBRCA+ HER2- Meme CA'},
            {'drug': 'Niraparib (Zejula)', 'fda_status': 'FDA Onaylı', 'indication': 'HRD+ Over Kanseri'}
        ],
        'primary_malignancies': ['Over Kanseri', 'Meme Kanseri', 'Prostat Kanseri', 'Pankreas Kanseri']
    },
    'BRCA2_DEFICIENCY': {
        'biomarker_gene': 'BRCA2',
        'biomarker_state': 'Loss-of-Function / Homologous Recombination Deficiency',
        'synthetic_lethal_target': 'PARP1',
        'target_full_name': 'Poly(ADP-Ribose) Polymerase 1',
        'target_uniprot': 'P09874',
        'depmap_chronos_score': -1.18,
        'depmap_essentiality': 'Extreme Dependency',
        'vulnerability_score': 96.0,
        'mechanistic_rationale': 'RAD51 yükleme defekti nedeniyle çift zincir DNA kırıklarını onaramayan hücreler, PARP inhibisyonuna aşırı duyarlıdır.',
        'clinical_therapies': [
            {'drug': 'Olaparib (Lynparza)', 'fda_status': 'FDA Onaylı', 'indication': 'Prostat, Meme ve Pankreas CA'},
            {'drug': 'Rucaparib (Rubraca)', 'fda_status': 'FDA Onaylı', 'indication': 'BRCA+ Over ve Prostat CA'}
        ],
        'primary_malignancies': ['Prostat Kanseri', 'Over Kanseri', 'Meme Kanseri']
    },
    'MSI_H_DMMR': {
        'biomarker_gene': 'MLH1 / MSH2 / MSH6 / PMS2',
        'biomarker_state': 'Mismatch Repair Deficiency (dMMR / MSI-H)',
        'synthetic_lethal_target': 'WRN',
        'target_full_name': 'Werner Syndrome RecQ Like Helicase',
        'target_uniprot': 'Q14191',
        'depmap_chronos_score': -1.45,
        'depmap_essentiality': 'Exclusive Selective Dependency in MSI-H',
        'vulnerability_score': 95.0,
        'mechanistic_rationale': 'Mikrosatellit instabilitesi (MSI-H) tümörlerinde non-B DNA (TA dinükleotid saç tokası) yapıları birikir. WRN helikaz bu yapıları açamadığında kromozom kırılmaları meydana gelir. MSS (sağlam) hücrelerde WRN vazgeçilmez değildir; yalnızca MSI-H tümörler için seçici sentetik letaldir.',
        'clinical_therapies': [
            {'drug': 'HRO761 (Novartis)', 'fda_status': 'Faz 1/2 Klinik Çalışma', 'indication': 'Gelişmiş MSI-H/dMMR Kolorektal ve Endometriyum CA'},
            {'drug': 'RO7589831 (Roche)', 'fda_status': 'Faz 1 Klinik Çalışma', 'indication': 'MSI-H Katı Tümörler'}
        ],
        'primary_malignancies': ['Kolorektal Kanser (MSI-H)', 'Endometriyum Kanseri', 'Mide Kanseri']
    },
    'MTAP_DELETION': {
        'biomarker_gene': 'MTAP',
        'biomarker_state': 'Homozygous Deletion (Co-deleted with CDKN2A at 9p21.3)',
        'synthetic_lethal_target': 'PRMT5 / MAT2A',
        'target_full_name': 'Protein Arginine Methyltransferase 5 / MAT2A',
        'target_uniprot': 'O14744',
        'depmap_chronos_score': -0.92,
        'depmap_essentiality': 'High Essentiality in MTAP-Null',
        'vulnerability_score': 88.0,
        'mechanistic_rationale': 'MTAP delesyonu sonucu biriken MTA (metiltiyo-adenozin), PRMT5 enzimini kısmen inhibe eder. Tümör hücreleri, SAM metabolizmasını sağlayan MAT2A ve arjinin metilasyonunu sürdüren PRMT5 inhibisyonuna karşı seçici olarak aşırı hassaslaşır.',
        'clinical_therapies': [
            {'drug': 'MRTX1719 (Mirati / BMS)', 'fda_status': 'Faz 1/2 Klinik Çalışma', 'indication': 'MTAP-delesyonlu Akciğer ve Pankreas CA'},
            {'drug': 'AG-270 (MAT2A İnhibitörü)', 'fda_status': 'Faz 2 Klinik Çalışma', 'indication': 'MTAP-Null Katı Tümörler'}
        ],
        'primary_malignancies': ['Küçük Hücreli Dışı Akciğer Kanseri', 'Pankreas Kanseri', 'Mezotelyoma', 'Glioblastoma']
    },
    'ARID1A_MUTATION': {
        'biomarker_gene': 'ARID1A',
        'biomarker_state': 'Inactivating Truncating Mutation (SWI/SNF Complex)',
        'synthetic_lethal_target': 'EZH2',
        'target_full_name': 'Enhancer of Zeste 2 Polycomb Repressive Complex 2 Subunit',
        'target_uniprot': 'Q15910',
        'depmap_chronos_score': -0.85,
        'depmap_essentiality': 'Strong Lineage Dependency',
        'vulnerability_score': 84.0,
        'mechanistic_rationale': 'ARID1A kaybı, epigenetik bivalent gen bölgelerinde SWI/SNF ve PRC2 arasındaki dengeyi bozar. EZH2 katalitik histon metilasyonunun inhibisyonu, ARID1A eksikliği olan hücrelerde diferansiyasyonu ve apoptozu indükler.',
        'clinical_therapies': [
            {'drug': 'Tazemetostat (Tazverik)', 'fda_status': 'FDA Onaylı', 'indication': 'Epiteloid Sarkom ve Foliküler Lenfoma (ARID1A ve EZH2 mutant)'}
        ],
        'primary_malignancies': ['Over Şeffaf Hücreli Karsinomu', 'Mide Kanseri', 'Endometriyum Kanseri', 'Mesane Kanseri']
    },
    'KRAS_ACTIVATION': {
        'biomarker_gene': 'KRAS',
        'biomarker_state': 'Oncogenic Codon 12/13/61 Mutation (G12C / G12D / G12V)',
        'synthetic_lethal_target': 'SHP2 (PTPN11) & SOS1',
        'target_full_name': 'Tyrosine-Protein Phosphatase Non-Receptor Type 11 / SOS1',
        'target_uniprot': 'Q06124',
        'depmap_chronos_score': -0.78,
        'depmap_essentiality': 'Adaptive Resistance Feedback Dependency',
        'vulnerability_score': 82.0,
        'mechanistic_rationale': 'Doğrudan KRAS(G12C) inhibitörleri (Sotorasib) uygulandığında tümör hızlıca RTK adaptif geri bildirim döngüsü ile KRAS’ı yeniden aktive eder. SHP2 veya SOS1 inhibisyonu ile bu geri besleme kesildiğinde KRAS blokajı derinleşir ve sentetik letal yanıt elde edilir.',
        'clinical_therapies': [
            {'drug': 'TNO155 (SHP2i) + Sotorasib', 'fda_status': 'Faz 1/2 Klinik Çalışma', 'indication': 'KRAS G12C NSCLC ve Kolorektal CA'},
            {'drug': 'BI-1701963 (SOS1i) + Trametinib', 'fda_status': 'Faz 1/2 Klinik Çalışma', 'indication': 'Pan-KRAS Mutant Solid Tümörler'}
        ],
        'primary_malignancies': ['Akciğer Kanseri (NSCLC)', 'Kolorektal Kanser', 'Pankreas Adenokarsinomu']
    }
}


class SyntheticLethalityEngine:
    """
    Evaluates genetic alterations to identify synthetic lethal target dependencies,
    DepMap Chronos essentiality scores, and tailored targeted therapies.
    """

    @classmethod
    def query_vulnerabilities(cls, gene_or_biomarker: str) -> Dict[str, Any]:
        """
        Finds synthetic lethal partner genes and DepMap dependency metrics
        for a given cancer driver or biomarker.
        """
        key = gene_or_biomarker.strip().upper()

        # Map common gene symbols to catalog keys
        symbol_map = {
            'BRCA1': 'BRCA1_DEFICIENCY',
            'BRCA2': 'BRCA2_DEFICIENCY',
            'MSI': 'MSI_H_DMMR',
            'MSI-H': 'MSI_H_DMMR',
            'DMMR': 'MSI_H_DMMR',
            'MLH1': 'MSI_H_DMMR',
            'MSH2': 'MSI_H_DMMR',
            'MTAP': 'MTAP_DELETION',
            'CDKN2A': 'MTAP_DELETION',
            'ARID1A': 'ARID1A_MUTATION',
            'KRAS': 'KRAS_ACTIVATION'
        }

        matched_key = symbol_map.get(key, key if key in SYNTHETIC_LETHAL_CATALOG else None)

        if not matched_key:
            return {
                'query': gene_or_biomarker,
                'status': 'no_cataloged_synthetic_lethality',
                'message': f"'{gene_or_biomarker}' geni için klinik veya DepMap onaylı doğrudan sentetik letal hedef bulunamadı.",
                'available_biomarkers': list(symbol_map.keys())
            }

        data = SYNTHETIC_LETHAL_CATALOG[matched_key]
        return {
            'query': gene_or_biomarker,
            'status': 'synthetic_lethal_identified',
            'biomarker_gene': data['biomarker_gene'],
            'biomarker_state': data['biomarker_state'],
            'synthetic_lethal_target': data['synthetic_lethal_target'],
            'target_full_name': data['target_full_name'],
            'target_uniprot': data['target_uniprot'],
            'depmap_chronos_score': data['depmap_chronos_score'],
            'depmap_essentiality': data['depmap_essentiality'],
            'vulnerability_score': data['vulnerability_score'],
            'mechanistic_rationale': data['mechanistic_rationale'],
            'clinical_therapies': data['clinical_therapies'],
            'primary_malignancies': data['primary_malignancies']
        }

    @classmethod
    def list_all_vulnerabilities(cls) -> List[Dict[str, Any]]:
        """Returns all cataloged synthetic lethal partnerships."""
        return [
            {
                'key': k,
                'biomarker': v['biomarker_gene'],
                'target': v['synthetic_lethal_target'],
                'depmap_chronos': v['depmap_chronos_score'],
                'vulnerability': v['vulnerability_score'],
                'top_drug': v['clinical_therapies'][0]['drug']
            }
            for k, v in SYNTHETIC_LETHAL_CATALOG.items()
        ]


if __name__ == '__main__':
    print("=" * 70)
    print("SYNTHETIC LETHALITY & DEPMAP CANCER DEPENDENCY ENGINE")
    print("=" * 70)

    for query in ['BRCA1', 'MSI-H', 'MTAP', 'KRAS']:
        res = SyntheticLethalityEngine.query_vulnerabilities(query)
        print(f"\n[Biomarker: {res['biomarker_gene']}]")
        print(f"  Synthetic Lethal Target: {res['synthetic_lethal_target']} ({res['target_full_name']})")
        print(f"  DepMap Chronos Score: {res['depmap_chronos_score']} -> Vulnerability: %{res['vulnerability_score']}")
        print(f"  Targeted Therapy: {res['clinical_therapies'][0]['drug']} ({res['clinical_therapies'][0]['fda_status']})")
