"""
Cancer Pathway & Signaling Cascade Engine.
Models oncogenic signal transduction networks, targeted drug perturbations,
downstream transcriptional flux, and feedback resistance mechanisms.
"""

from typing import Dict, List, Any, Optional

CANONICAL_PATHWAYS: Dict[str, Dict[str, Any]] = {
    'MAPK_ERK': {
        'id': 'MAPK_ERK',
        'name': 'RTK / RAS / RAF / MEK / ERK (MAPK Kaskadı)',
        'description': 'Hücre çoğalması, farklılaşması ve sağkalımını yöneten en temel mitogenik sinyal yolağıdır. Katı tümörlerin %30-40\'ında aşırı aktiftir.',
        'nodes': [
            {'id': 'EGFR', 'name': 'EGFR (RTK)', 'role': 'Reseptör Tirozin Kinaz', 'drugs': ['Osimertinib', 'Erlotinib', 'Gefitinib', 'Cetuximab']},
            {'id': 'GRB2_SOS', 'name': 'GRB2 / SOS1', 'role': 'Adaptör & GEF Kompleksi', 'drugs': ['BI-3406 (SOS1i)']},
            {'id': 'KRAS', 'name': 'KRAS (GTPase)', 'role': 'Moleküler Anahtar (GTP/GDP)', 'drugs': ['Sotorasib (G12C)', 'Adagrasib (G12C)', 'MRTX1133 (G12D)']},
            {'id': 'BRAF', 'name': 'BRAF (Ser/Thr Kinaz)', 'role': 'RAF Kinaz Ailesi', 'drugs': ['Dabrafenib', 'Encorafenib', 'Vemurafenib']},
            {'id': 'MEK1_2', 'name': 'MEK1 / MEK2', 'role': 'İkili Özgüllüklü Kinaz', 'drugs': ['Trametinib', 'Binimetinib', 'Cobimetinib']},
            {'id': 'ERK1_2', 'name': 'ERK1 / ERK2', 'role': 'Efektör MAPK Kinaz', 'drugs': ['Ulixertinib (ERK1/2i)']},
            {'id': 'MYC_ELK1', 'name': 'MYC / ELK1 / FOS', 'role': 'Transkripsiyon Faktörleri', 'drugs': []}
        ],
        'edges': [
            {'from': 'EGFR', 'to': 'GRB2_SOS', 'type': 'activation'},
            {'from': 'GRB2_SOS', 'to': 'KRAS', 'type': 'activation'},
            {'from': 'KRAS', 'to': 'BRAF', 'type': 'activation'},
            {'from': 'BRAF', 'to': 'MEK1_2', 'type': 'activation'},
            {'from': 'MEK1_2', 'to': 'ERK1_2', 'type': 'activation'},
            {'from': 'ERK1_2', 'to': 'MYC_ELK1', 'type': 'activation'},
            {'from': 'ERK1_2', 'to': 'EGFR', 'type': 'negative_feedback'}
        ],
        'primary_cancers': ['Akciğer Adenokarsinomu', 'Kolorektal Kanser', 'Malign Melanom', 'Pankreas Duktal Karsinomu'],
        'resistance_mechanisms': [
            'KRAS G12C tek ajan kullanımında EGFR aracılı hızlı feedback aktivasyonu (CodeBreak-300)',
            'BRAF V600E kolon kanserinde EGFR feedback yeniden aktivasyonu (BEACON)',
            'MET amplifikasyonu veya HER2 aşırı ekspresyonu ile kaskad baypası'
        ]
    },
    'PI3K_AKT_MTOR': {
        'id': 'PI3K_AKT_MTOR',
        'name': 'PI3K / AKT / PTEN / mTOR Sağkalım Yolağı',
        'description': 'Hücresel büyüme, protein translasyonu, glukoz metabolizması ve apoptoz direncini kontrol eden majör anabolik sinyal eksenidir.',
        'nodes': [
            {'id': 'RTK_IRS', 'name': 'RTK / IRS1', 'role': 'Membran Reseptör Kompleksi', 'drugs': ['Trastuzumab', 'Pertuzumab']},
            {'id': 'PI3K', 'name': 'PI3K (p110α/p85)', 'role': 'Lipid Kinaz (PIP2 -> PIP3)', 'drugs': ['Alpelisib (p110α)', 'Inavolisib']},
            {'id': 'PTEN', 'name': 'PTEN', 'role': 'Tümör Baskılayıcı Fosfataz (PIP3 -> PIP2)', 'drugs': []},
            {'id': 'AKT', 'name': 'AKT1 / AKT2', 'role': 'Serin/Treonin Kinaz', 'drugs': ['Capivasertib', 'Ipatasertib']},
            {'id': 'TSC_RHEB', 'name': 'TSC1/2 & RHEB', 'role': 'GTPase Düzenleyici', 'drugs': []},
            {'id': 'MTORC1', 'name': 'mTORC1 (Raptor)', 'role': 'Besin & Enerji Sensörü', 'drugs': ['Everolimus', 'Temsirolimus', 'Sirolimus']},
            {'id': 'S6K_4EBP1', 'name': 'p70S6K / 4E-BP1', 'role': 'Translasyon Başlatıcı', 'drugs': []}
        ],
        'edges': [
            {'from': 'RTK_IRS', 'to': 'PI3K', 'type': 'activation'},
            {'from': 'PI3K', 'to': 'AKT', 'type': 'activation'},
            {'from': 'PTEN', 'to': 'AKT', 'type': 'inhibition'},
            {'from': 'AKT', 'to': 'TSC_RHEB', 'type': 'inhibition'},
            {'from': 'TSC_RHEB', 'to': 'MTORC1', 'type': 'activation'},
            {'from': 'MTORC1', 'to': 'S6K_4EBP1', 'type': 'activation'}
        ],
        'primary_cancers': ['Meme Karsinomu (HR+/HER2-)', 'Endometrium Karsinomu', 'Prostat Karsinomu', 'Glioblastoma'],
        'resistance_mechanisms': [
            'PTEN kaybı durumunda PI3K alfa inhibitörlerine primer direnç',
            'mTORC1 inhibisyonu sonrası IRS-1 aracılı AKT fosforilasyonu geri tepmesi',
            'HER3 ve FGFR aşırı ekspresyonu'
        ]
    },
    'DDR_HRR': {
        'id': 'DDR_HRR',
        'name': 'DNA Hasar Yanıtı (DDR) & Homolog Rekombinasyon (HRR)',
        'description': 'Çift iplikçikli DNA kırıklarının (DSB) hatasız onarımını sağlayan genom koruyucu sistemdir. Sentetik letalite tedavisinin temelidir.',
        'nodes': [
            {'id': 'DSB_LESION', 'name': 'DNA Çift İplik Kırığı (DSB)', 'role': 'Genomik Hasar', 'drugs': ['Platin Kemoterapisi (Sisplatin/Karboplatin)']},
            {'id': 'ATM_ATR', 'name': 'ATM / ATR Kinazlar', 'role': 'Hasar Sensör Kinazları', 'drugs': ['Ceralasertib (ATRi)']},
            {'id': 'CHEK1_2', 'name': 'CHK1 / CHK2', 'role': 'Kontrol Noktası Kinazları', 'drugs': []},
            {'id': 'BRCA1_2', 'name': 'BRCA1 / BRCA2', 'role': 'HRR Tamir İskelet Proteinleri', 'drugs': []},
            {'id': 'RAD51', 'name': 'RAD51', 'role': 'Rekombinaz Filament Formasyonu', 'drugs': []},
            {'id': 'PARP1', 'name': 'PARP1 / PARP2', 'role': 'Tek İplik Tamir Enzimi', 'drugs': ['Olaparib', 'Talazoparib', 'Niraparib', 'Rucaparib']}
        ],
        'edges': [
            {'from': 'DSB_LESION', 'to': 'ATM_ATR', 'type': 'activation'},
            {'from': 'ATM_ATR', 'to': 'CHEK1_2', 'type': 'activation'},
            {'from': 'CHEK1_2', 'to': 'BRCA1_2', 'type': 'activation'},
            {'from': 'BRCA1_2', 'to': 'RAD51', 'type': 'activation'},
            {'from': 'PARP1', 'to': 'DSB_LESION', 'type': 'inhibition'}
        ],
        'primary_cancers': ['Yumurtalık (Over) Kanseri', 'Üçlü Negatif Meme Kanseri', 'mCRPC Prostat', 'Pankreas Adeno'],
        'resistance_mechanisms': [
            'BRCA1/2 reversiyon mutasyonları (açık okuma çerçevesinin yeniden kurulması)',
            'PARP1 mutasyonları (tuzaklanma afinitesinin azalması)',
            '53BP1 kaybı ile alternatif uç birleştirme yoluyla HR yetkinliğinin kısmi restorasyonu'
        ]
    },
    'CELL_CYCLE_CDK': {
        'id': 'CELL_CYCLE_CDK',
        'name': 'Hücre Döngüsü G1/S Kontrolü (CDK4/6 - RB1 Aksı)',
        'description': 'Hücrenin kontrolsüz bölünmesini önleyen G1/S fazı kısıtlama noktasıdır. Tümörlerin çoğunda p16 kaybı veya Cyclin D1 artışı görülür.',
        'nodes': [
            {'id': 'MITOGEN_SIGNAL', 'name': 'Mitojenik Sinyal (ER/RTK)', 'role': 'Giriş Uyaranı', 'drugs': ['Fulvestrant', 'Letrozol']},
            {'id': 'CYCLIN_D1', 'name': 'Cyclin D1 (CCND1)', 'role': 'Döngü Düzenleyici Siklin', 'drugs': []},
            {'id': 'CDKN2A_P16', 'name': 'p16INK4A (CDKN2A)', 'role': 'Doğal CDK4/6 İnhibitörü', 'drugs': []},
            {'id': 'CDK4_6', 'name': 'CDK4 / CDK6', 'role': 'Siklin Bağımlı Kinazlar', 'drugs': ['Palbociclib', 'Ribociclib', 'Abemaciclib']},
            {'id': 'RB1', 'name': 'Retinoblastoma (pRB)', 'role': 'Tümör Baskılayıcı (Fosforillenmemiş)', 'drugs': []},
            {'id': 'E2F', 'name': 'E2F Transkripsiyon Faktörleri', 'role': 'S-Fazı Genlerini Başlatır', 'drugs': []}
        ],
        'edges': [
            {'from': 'MITOGEN_SIGNAL', 'to': 'CYCLIN_D1', 'type': 'activation'},
            {'from': 'CYCLIN_D1', 'to': 'CDK4_6', 'type': 'activation'},
            {'from': 'CDKN2A_P16', 'to': 'CDK4_6', 'type': 'inhibition'},
            {'from': 'CDK4_6', 'to': 'RB1', 'type': 'inhibition'},
            {'from': 'RB1', 'to': 'E2F', 'type': 'inhibition'}
        ],
        'primary_cancers': ['HR+/HER2- Metastatik Meme Kanseri', 'Malign Melanom', 'Akciğer Skuamöz Karsinomu'],
        'resistance_mechanisms': [
            'RB1 homozigot delesyonu veya inaktivasyonu (CDK4/6i\'ye mutlak direnç)',
            'CCNE1 (Cyclin E1) amplifikasyonu ile CDK2 baypası',
            'FGFR1 veya AKT aşırı aktivasyonu'
        ]
    },
    'IMMUNE_CHECKPOINT': {
        'id': 'IMMUNE_CHECKPOINT',
        'name': 'İmmün Kontrol Noktası & T-Hücre Anergisi (PD-1 / PD-L1)',
        'description': 'Tümör mikroçevresinde sitotoksik CD8+ T hücrelerinin tükenmesini ve immün kaçışı sağlayan ligand-reseptör eksenidir.',
        'nodes': [
            {'id': 'TCR_MHC', 'name': 'TCR : Neoantijen/MHC-I', 'role': 'T-Hücre Tanıma Kompleksi', 'drugs': []},
            {'id': 'TMB_NEOAG', 'name': 'TMB & Neoantijen Yükü', 'role': 'İmmünojenik Uyarım', 'drugs': []},
            {'id': 'PD1', 'name': 'PD-1 (CD279 - T-Hücresi)', 'role': 'İnhibitör Reseptör', 'drugs': ['Pembrolizumab', 'Nivolumab', 'Cemiplimab']},
            {'id': 'PDL1', 'name': 'PD-L1 (CD274 - Tümör)', 'role': 'İmmün Kaçış Ligandı', 'drugs': ['Atezolizumab', 'Durvalumab', 'Avelumab']},
            {'id': 'SHP2', 'name': 'SHP-2 Fosfataz', 'role': 'TCR Sinyal Söndürücü', 'drugs': []},
            {'id': 'CYTOTOXIC_KILL', 'name': 'Tümör Hücresi Lizisi & Apoptoz', 'role': 'İmmün Sitotoksisite', 'drugs': []}
        ],
        'edges': [
            {'from': 'TMB_NEOAG', 'to': 'TCR_MHC', 'type': 'activation'},
            {'from': 'TCR_MHC', 'to': 'CYTOTOXIC_KILL', 'type': 'activation'},
            {'from': 'PDL1', 'to': 'PD1', 'type': 'activation'},
            {'from': 'PD1', 'to': 'SHP2', 'type': 'activation'},
            {'from': 'SHP2', 'to': 'CYTOTOXIC_KILL', 'type': 'inhibition'}
        ],
        'primary_cancers': ['Küçük Hücreli Dışı Akciğer Kanseri', 'Malign Melanom', 'MSI-H Kolorektal', 'Mesane Ürotelyal'],
        'resistance_mechanisms': [
            'B2M (Beta-2 Mikroglobulin) mutasyonu ile MHC-I yüzey ekspresyon kaybı',
            'JAK1/JAK2 inaktivasyonu ile interferon-gama yanıtsızlığı',
            'LAG-3 veya TIGIT alternatif kontrol noktalarının aşırı ekspresyonu'
        ]
    }
}


class CancerPathwayEngine:
    """
    Cancer Signaling Pathway Simulation & Perturbation Engine.
    Simulates pathway flux, oncogenic driver alterations, and targeted drug effects.
    """

    def __init__(self):
        self.pathways = CANONICAL_PATHWAYS

    def list_pathways(self) -> List[Dict[str, Any]]:
        """Returns metadata for all available canonical signaling pathways."""
        results = []
        for pid, pdata in self.pathways.items():
            results.append({
                'id': pid,
                'name': pdata['name'],
                'description': pdata['description'],
                'primary_cancers': pdata['primary_cancers'],
                'node_count': len(pdata['nodes']),
                'edge_count': len(pdata['edges'])
            })
        return results

    def get_pathway(self, pathway_id: str) -> Optional[Dict[str, Any]]:
        """Returns full network definition for a given pathway ID."""
        return self.pathways.get(pathway_id)

    def simulate_pathway_flux(
        self,
        pathway_id: str,
        alterations: Optional[List[str]] = None,
        drugs: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Simulates the biochemical flux and phenotypic outcome of a pathway.
        Args:
            pathway_id: Identifier of the pathway (e.g., 'MAPK_ERK')
            alterations: List of present mutations (e.g., ['KRAS_G12C', 'MET_AMP'])
            drugs: List of applied targeted therapies (e.g., ['Sotorasib'])
        Returns:
            Dictionary with pathway activity (0-100), proliferation score, apoptosis score,
            drug efficacy status, and resistance warnings.
        """
        pathway = self.pathways.get(pathway_id)
        if not pathway:
            raise ValueError(f"Unknown pathway: {pathway_id}")

        alterations = [a.upper() for a in (alterations or [])]
        drugs = [d.upper() for d in (drugs or [])]

        # Base basal activity
        base_activity = 40.0
        activation_bonus = 0.0
        inhibition_effect = 0.0
        resistance_flags = []
        clinical_notes = []

        if pathway_id == 'MAPK_ERK':
            # Driver mutations
            if any('KRAS' in a for a in alterations):
                activation_bonus += 45.0
                clinical_notes.append("Onkojenik KRAS mutasyonu MAPK yolağını konstitütif olarak aktifler.")
            if any('BRAF' in a for a in alterations):
                activation_bonus += 50.0
                clinical_notes.append("BRAF V600E monomerik sinyali azami düzeye çıkarır.")
            if any('EGFR' in a for a in alterations):
                activation_bonus += 35.0

            # Drug treatments
            has_kras_i = any('SOTORASIB' in d or 'ADAGRASIB' in d for d in drugs)
            has_braf_i = any('DABRAFENIB' in d or 'ENCORAFENIB' in d for d in drugs)
            has_mek_i = any('TRAMETINIB' in d or 'BINIMETINIB' in d for d in drugs)
            has_anti_egfr = any('CETUXIMAB' in d or 'PANITUMUMAB' in d or 'OSIMERTINIB' in d for d in drugs)

            if has_kras_i:
                if has_anti_egfr:
                    inhibition_effect += 75.0
                    clinical_notes.append("Kovalent KRAS G12C + Anti-EGFR ikili blokajı sinerjik baskılama sağlar (CodeBreak-300).")
                else:
                    inhibition_effect += 40.0
                    resistance_flags.append("Anti-EGFR eklenmediğinde EGFR aracılı hızlı feedback aktivasyonu riski mevcuttur.")

            if has_braf_i and has_mek_i:
                inhibition_effect += 80.0
                clinical_notes.append("BRAF + MEK ikili kaskad inhibisyonu paradoksal aktivasyonu önler.")
            elif has_braf_i:
                inhibition_effect += 45.0
                resistance_flags.append("Tek ajan BRAF inhibisyonu kolon kanserinde yetersizdir; MEK veya EGFR blokajı gereklidir.")

            if any('MET' in a for a in alterations):
                inhibition_effect *= 0.4
                resistance_flags.append("Eş zamanlı MET amplifikasyonu TKI direncine yol açar; MET inhibitörü (Capmatinib/Savolitinib) gereklidir.")

        elif pathway_id == 'PI3K_AKT_MTOR':
            if any('PIK3CA' in a for a in alterations):
                activation_bonus += 45.0
            if any('PTEN' in a for a in alterations):
                activation_bonus += 40.0

            has_pi3k_i = any('ALPELISIB' in d or 'INAVOLISIB' in d for d in drugs)
            has_akt_i = any('CAPIVASERTIB' in d for d in drugs)
            has_mtor_i = any('EVEROLIMUS' in d for d in drugs)

            if has_pi3k_i:
                if any('PTEN' in a for a in alterations):
                    inhibition_effect += 20.0
                    resistance_flags.append("PTEN kaybı varlığında tek başına PI3K alfa inhibitörüne primer direnç gelişir.")
                else:
                    inhibition_effect += 65.0
            if has_akt_i:
                inhibition_effect += 70.0
            if has_mtor_i:
                inhibition_effect += 55.0

        elif pathway_id == 'DDR_HRR':
            is_hrd = any('BRCA' in a or 'CDK12' in a or 'PALB2' in a for a in alterations)
            if is_hrd:
                activation_bonus -= 30.0  # Defective repair capacity
                clinical_notes.append("Homolog rekombinasyon tamir eksikliği (HRD) saptandı; sentetik letalite penceresi açıktır.")

            has_parp_i = any('OLAPARIB' in d or 'TALAZOPARIB' in d or 'NIRAPARIB' in d for d in drugs)
            has_platinum = any('SISPLATIN' in d or 'KARBOPLATIN' in d or 'PLATIN' in d for d in drugs)

            if is_hrd and has_parp_i:
                inhibition_effect += 85.0
                clinical_notes.append("HRD + PARP İnhibitörü sentetik letalite mekanizması ile çift iplik kırığı krizine yol açar.")
            elif has_parp_i:
                inhibition_effect += 25.0
                resistance_flags.append("HR yetkin (BRCA-wt) tümörlerde PARP inhibitörü monoterapisi sınırlı etkiye sahiptir.")

        elif pathway_id == 'CELL_CYCLE_CDK':
            if any('CDKN2A' in a or 'CCND1' in a for a in alterations):
                activation_bonus += 45.0
            is_rb1_lost = any('RB1' in a for a in alterations)
            if is_rb1_lost:
                activation_bonus += 50.0

            has_cdk46_i = any('PALBOCICLIB' in d or 'RIBOCICLIB' in d or 'ABEMACICLIB' in d for d in drugs)
            if has_cdk46_i:
                if is_rb1_lost:
                    inhibition_effect += 5.0
                    resistance_flags.append("RB1 mutasyonu/kaybı CDK4/6 inhibitörlerine mutlak direnç oluşturur.")
                else:
                    inhibition_effect += 75.0
                    clinical_notes.append("CDK4/6 inhibisyonu ile pRB hipofosforile tutularak E2F salınımı bloke edilir.")

        elif pathway_id == 'IMMUNE_CHECKPOINT':
            is_hot_tumor = any('MSI' in a or 'TMB' in a for a in alterations)
            if is_hot_tumor:
                activation_bonus += 40.0

            has_io = any('PEMBROLIZUMAB' in d or 'NIVOLUMAB' in d or 'ATEZOLIZUMAB' in d for d in drugs)
            if has_io:
                if is_hot_tumor:
                    inhibition_effect += 85.0
                    clinical_notes.append("Yüksek neoantijen yükü + PD-1/PD-L1 blokajı ile güçlü sitotoksik T-hücre reaktivasyonu.")
                else:
                    inhibition_effect += 45.0

        # Calculate composite pathway score
        raw_flux = base_activity + activation_bonus - inhibition_effect
        final_flux = max(5.0, min(100.0, raw_flux))
        proliferation_score = max(5.0, min(100.0, final_flux * 0.95))
        apoptosis_score = max(5.0, min(100.0, (100.0 - final_flux) * 0.90 + (inhibition_effect * 0.2)))

        return {
            'pathway_id': pathway_id,
            'pathway_name': pathway['name'],
            'final_flux': round(final_flux, 1),
            'proliferation_score': round(proliferation_score, 1),
            'apoptosis_score': round(apoptosis_score, 1),
            'alterations_applied': alterations,
            'drugs_applied': drugs,
            'resistance_flags': resistance_flags,
            'clinical_notes': clinical_notes,
            'status': 'BASKILANDI' if final_flux < 35.0 else ('KISMEN KONTROLDE' if final_flux < 60.0 else 'AŞIRI AKTİF / DİRENÇLİ')
        }
