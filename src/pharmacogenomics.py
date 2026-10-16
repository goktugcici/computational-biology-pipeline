"""
Pharmacogenomics (PGx) Safety & Clinical Drug-Gene Toxicity Engine.
Implements CPIC (Clinical Pharmacogenetics Implementation Consortium) Level A/B
and PharmGKB clinical guidelines for precision oncology and supportive therapies.
"""

from typing import Dict, List, Any, Optional, Tuple

# CPIC Level A/B Guidelines and Diplotype Definitions
CPIC_GENE_GUIDELINES: Dict[str, Dict[str, Any]] = {
    'DPYD': {
        'gene_symbol': 'DPYD',
        'full_name': 'Dihydropyrimidine Dehydrogenase',
        'chromosome': '1p21.3',
        'primary_drugs': ['Fluorouracil (5-FU)', 'Capecitabine (Xeloda)', 'Tegafur'],
        'therapeutic_area': 'Gastrointestinal & Breast Oncology',
        'metabolic_pathway': 'Hepatic catabolism of fluoropyrimidines (inactivates >80% of dose)',
        'alleles': {
            '*1': {'activity_score': 1.0, 'description': 'Normal function (Wild-Type)'},
            '*2A': {'activity_score': 0.0, 'rs_id': 'rs3918290', 'hgvs': 'c.1905+1G>A', 'description': 'No function (Splice site variant, lethal toxicity risk)'},
            '*13': {'activity_score': 0.0, 'rs_id': 'rs55886062', 'hgvs': 'c.1679T>G (p.I560S)', 'description': 'No function (Severe myelosuppression risk)'},
            'c.2846A>T': {'activity_score': 0.5, 'rs_id': 'rs67376798', 'hgvs': 'p.D949V', 'description': 'Decreased function (Moderate to severe toxicity risk)'},
            'HapB3': {'activity_score': 0.5, 'rs_id': 'rs75017182', 'hgvs': 'c.1129-5923C>G', 'description': 'Decreased function (Deep intronic variant)'}
        },
        'phenotypes': {
            'NM': {'name': 'Normal Metabolizer', 'min_as': 2.0, 'max_as': 2.0, 'risk': 'Standard Risk', 'color': '#10b981'},
            'IM': {'name': 'Intermediate Metabolizer', 'min_as': 1.0, 'max_as': 1.5, 'risk': 'Moderate/High Toxicity Risk', 'color': '#f59e0b'},
            'PM': {'name': 'Poor Metabolizer', 'min_as': 0.0, 'max_as': 0.5, 'risk': 'Fatal / Extreme Toxicity Risk (Contraindicated)', 'color': '#ef4444'}
        }
    },
    'UGT1A1': {
        'gene_symbol': 'UGT1A1',
        'full_name': 'UDP Glucuronosyltransferase Family 1 Member A1',
        'chromosome': '2q37.1',
        'primary_drugs': ['Irinotecan (Camptosar, FOLFIRI, FOLFIRINOX)'],
        'therapeutic_area': 'Colorectal & Pancreatic Oncology',
        'metabolic_pathway': 'Glucuronidation of active toxic metabolite SN-38 to inactive SN-38G',
        'alleles': {
            '*1': {'function': 'Normal', 'description': '(TA)6 repeat promoter allele'},
            '*28': {'function': 'Decreased', 'rs_id': 'rs8175347', 'description': '(TA)7 repeat promoter allele (Gilbert Syndrome allele)'},
            '*6': {'function': 'Decreased', 'rs_id': 'rs4148323', 'description': 'c.211G>A (p.G71R) common in East Asian populations'}
        },
        'phenotypes': {
            'NM': {'name': 'Normal Metabolizer (*1/*1)', 'risk': 'Standard Risk', 'color': '#10b981'},
            'IM': {'name': 'Intermediate Metabolizer (*1/*28 or *1/*6)', 'risk': 'Moderate Neutropenia Risk', 'color': '#f59e0b'},
            'PM': {'name': 'Poor Metabolizer (*28/*28 or *6/*6)', 'risk': 'Severe Grade 4 Neutropenia & Diarrhea Risk', 'color': '#ef4444'}
        }
    },
    'TPMT': {
        'gene_symbol': 'TPMT',
        'full_name': 'Thiopurine S-Methyltransferase',
        'chromosome': '6p22.3',
        'primary_drugs': ['6-Mercaptopurine (6-MP)', 'Azathioprine', 'Thioguanine'],
        'therapeutic_area': 'Leukemia (ALL) & Immunosuppression',
        'metabolic_pathway': 'S-methylation of thiopurines, diverting from cytotoxic TGN accumulation',
        'alleles': {
            '*1': {'function': 'Normal', 'description': 'Wild-type'},
            '*2': {'function': 'No function', 'rs_id': 'rs1800462', 'hgvs': 'c.238G>C (p.A80P)'},
            '*3A': {'function': 'No function', 'rs_id': 'rs1800460 / rs1142345', 'hgvs': 'c.460G>A + c.719A>G'},
            '*3C': {'function': 'No function', 'rs_id': 'rs1142345', 'hgvs': 'c.719A>G (p.Y240C)'}
        },
        'phenotypes': {
            'NM': {'name': 'Normal Metabolizer', 'risk': 'Standard Risk', 'color': '#10b981'},
            'IM': {'name': 'Intermediate Metabolizer', 'risk': 'High Myelosuppression Risk', 'color': '#f59e0b'},
            'PM': {'name': 'Poor Metabolizer', 'risk': 'Life-Threatening Pancytopenia', 'color': '#ef4444'}
        }
    },
    'CYP2D6': {
        'gene_symbol': 'CYP2D6',
        'full_name': 'Cytochrome P450 Family 2 Subfamily D Member 6',
        'chromosome': '22q13.2',
        'primary_drugs': ['Tamoxifen (Nolvadex)', 'Ondansetron', 'Tramadol', 'Codeine'],
        'therapeutic_area': 'ER+ Breast Oncology & Supportive Care',
        'metabolic_pathway': 'Bioactivation of Tamoxifen to potent anti-estrogen Endoxifen (100x higher affinity)',
        'alleles': {
            '*1': {'function': 'Normal', 'description': 'Wild-type'},
            '*2': {'function': 'Normal', 'description': 'Normal function'},
            '*4': {'function': 'No function', 'rs_id': 'rs3892097', 'description': 'Splice defect (most common in Europeans)'},
            '*10': {'function': 'Decreased', 'rs_id': 'rs1065852', 'description': 'Decreased function (common in East Asians)'},
            '*41': {'function': 'Decreased', 'rs_id': 'rs28371725', 'description': 'Decreased function'}
        },
        'phenotypes': {
            'UM': {'name': 'Ultrarapid Metabolizer', 'risk': 'Standard Efficacy (Increased prodrug conversion)', 'color': '#38bdf8'},
            'NM': {'name': 'Normal Metabolizer', 'risk': 'Standard Efficacy', 'color': '#10b981'},
            'IM': {'name': 'Intermediate Metabolizer', 'risk': 'Suboptimal Endoxifen Levels', 'color': '#f59e0b'},
            'PM': {'name': 'Poor Metabolizer', 'risk': 'Severe Efficacy Failure (Therapeutic Resistance)', 'color': '#ef4444'}
        }
    }
}


class PharmacogenomicsEngine:
    """
    Evaluates patient pharmacogenomic diplotypes to predict drug toxicity,
    metabolic clearance, and CPIC guideline dosing recommendations.
    """

    @classmethod
    def evaluate_dpyd(cls, allele_1: str = '*1', allele_2: str = '*1') -> Dict[str, Any]:
        """
        Calculates DPYD Activity Score and CPIC dosing guideline for Fluoropyrimidines (5-FU, Capecitabine).
        Activity Score range: 0.0 to 2.0.
        """
        alleles = CPIC_GENE_GUIDELINES['DPYD']['alleles']
        as_1 = alleles.get(allele_1, {}).get('activity_score', 1.0)
        as_2 = alleles.get(allele_2, {}).get('activity_score', 1.0)
        total_as = round(as_1 + as_2, 1)

        if total_as == 2.0:
            phenotype_code = 'NM'
            dose_recommendation = '100% Standart Doz. Başlangıç dozu azaltımı gerekmez.'
            action_code = 'STANDARD_DOSE'
            severity = 'LOW'
        elif total_as >= 1.0:
            phenotype_code = 'IM'
            dose_recommendation = '%50 Başlangıç Dozu Azaltımı Önerilir. Terapötik ilaç izlemi (TDM) veya toksisiteye göre titrasyon yapınız.'
            action_code = 'DOSE_REDUCTION_50'
            severity = 'HIGH'
        else:
            phenotype_code = 'PM'
            dose_recommendation = 'KONTRENDİKE! 5-FU ve Kapesitabin kullanımından tamamen kaçınınız. Ölümcül nötropeni ve mukozit riski.'
            action_code = 'CONTRAINDICATED'
            severity = 'CRITICAL'

        pheno = CPIC_GENE_GUIDELINES['DPYD']['phenotypes'][phenotype_code]

        return {
            'gene': 'DPYD',
            'diplotype': f"{allele_1}/{allele_2}",
            'activity_score': total_as,
            'phenotype': pheno['name'],
            'phenotype_code': phenotype_code,
            'risk_level': pheno['risk'],
            'toxicity_severity': severity,
            'color_code': pheno['color'],
            'affected_drugs': CPIC_GENE_GUIDELINES['DPYD']['primary_drugs'],
            'clinical_action': dose_recommendation,
            'action_code': action_code,
            'guideline_source': 'CPIC Level A / PharmGKB 2024 Guidelines'
        }

    @classmethod
    def evaluate_ugt1a1(cls, diplotype: str = '*1/*1') -> Dict[str, Any]:
        """
        Evaluates UGT1A1 promoter repeats for Irinotecan (Camptosar) toxicity risk.
        """
        d = diplotype.replace(" ", "")
        if d in ['*28/*28', '*6/*6', '*28/*6']:
            phenotype_code = 'PM'
            dose_recommendation = 'İrinotekan başlangıç dozunda en az %25-30 azaltım önerilir. Ağır Grade 4 nötropeni ve şiddetli diyare riski yüksektir.'
            severity = 'HIGH'
        elif any(x in d for x in ['*28', '*6']):
            phenotype_code = 'IM'
            dose_recommendation = 'Standart başlangıç dozu uygulanabilir; ancak nötropeni takibi haftalık tam kan sayımı ile yapılmalıdır.'
            severity = 'MODERATE'
        else:
            phenotype_code = 'NM'
            dose_recommendation = 'Standart irinotekan dozu (FOLFIRI / FOLFIRINOX protokolü) güvenle uygulanabilir.'
            severity = 'LOW'

        pheno = CPIC_GENE_GUIDELINES['UGT1A1']['phenotypes'][phenotype_code]

        return {
            'gene': 'UGT1A1',
            'diplotype': diplotype,
            'phenotype': pheno['name'],
            'phenotype_code': phenotype_code,
            'risk_level': pheno['risk'],
            'toxicity_severity': severity,
            'color_code': pheno['color'],
            'affected_drugs': CPIC_GENE_GUIDELINES['UGT1A1']['primary_drugs'],
            'clinical_action': dose_recommendation,
            'guideline_source': 'CPIC Level A Guideline for Irinotecan'
        }

    @classmethod
    def evaluate_cyp2d6_tamoxifen(cls, diplotype: str = '*1/*1') -> Dict[str, Any]:
        """
        Evaluates CYP2D6 bioactivation of Tamoxifen to active Endoxifen in ER+ Breast Cancer.
        """
        d = diplotype.replace(" ", "")
        if d in ['*4/*4', '*4/*5', '*5/*5', '*4/*10']:
            phenotype_code = 'PM'
            recommendation = 'Tamoksifen yetersiz aktif metabolit (Endoksifen) üretir. Postmenopozal hastalarda Aromataz İnhibitörüne (Anastrozol/Letrozol) geçilmesi önerilir.'
            severity = 'CRITICAL'
        elif any(x in d for x in ['*4', '*10', '*41']):
            phenotype_code = 'IM'
            recommendation = 'Ilımlı endoksifen seviyeleri. Tamoksifen dozu artırılabilir (örn. 40 mg/gün) veya aromataz inhibitörü değerlendirilmelidir.'
            severity = 'MODERATE'
        else:
            phenotype_code = 'NM'
            recommendation = 'Standart 20 mg/gün Tamoksifen tedavisi ile optimal endoksifen serum konsantrasyonu sağlanır.'
            severity = 'LOW'

        pheno = CPIC_GENE_GUIDELINES['CYP2D6']['phenotypes'][phenotype_code]

        return {
            'gene': 'CYP2D6',
            'diplotype': diplotype,
            'phenotype': pheno['name'],
            'phenotype_code': phenotype_code,
            'risk_level': pheno['risk'],
            'efficacy_impact': 'Subtherapeutic Resistance' if phenotype_code == 'PM' else 'Optimal Bioactivation',
            'color_code': pheno['color'],
            'affected_drugs': ['Tamoxifen (Nolvadex)'],
            'clinical_action': recommendation,
            'guideline_source': 'CPIC Level A Guideline for Tamoxifen'
        }

    @classmethod
    def profile_all_oncology_pgx(cls, dpyd_diplotype: Tuple[str, str] = ('*1', '*1'),
                                  ugt1a1_diplotype: str = '*1/*1',
                                  cyp2d6_diplotype: str = '*1/*1') -> Dict[str, Any]:
        """Profiles complete oncology pharmacogenomics panel for a patient."""
        return {
            'dpyd_fluoropyrimidine': cls.evaluate_dpyd(dpyd_diplotype[0], dpyd_diplotype[1]),
            'ugt1a1_irinotecan': cls.evaluate_ugt1a1(ugt1a1_diplotype),
            'cyp2d6_tamoxifen': cls.evaluate_cyp2d6_tamoxifen(cyp2d6_diplotype)
        }


if __name__ == '__main__':
    print("=" * 70)
    print("PHARMACOGENOMICS (PGx) CLINICAL TOXICITY ENGINE")
    print("=" * 70)

    # Test DPYD *1/*2A (Intermediate Metabolizer, high risk of lethal 5-FU toxicity)
    dpyd_res = PharmacogenomicsEngine.evaluate_dpyd('*1', '*2A')
    print(f"\n[DPYD Fluoropyrimidine Check: {dpyd_res['diplotype']}]")
    print(f"  Activity Score: {dpyd_res['activity_score']} -> {dpyd_res['phenotype']}")
    print(f"  Toxicity Risk: {dpyd_res['risk_level']} ({dpyd_res['toxicity_severity']})")
    print(f"  Clinical Action: {dpyd_res['clinical_action']}")

    # Test UGT1A1 *28/*28 (Irinotecan poor metabolizer)
    ugt_res = PharmacogenomicsEngine.evaluate_ugt1a1('*28/*28')
    print(f"\n[UGT1A1 Irinotecan Check: {ugt_res['diplotype']}]")
    print(f"  Phenotype: {ugt_res['phenotype']} -> {ugt_res['risk_level']}")
    print(f"  Clinical Action: {ugt_res['clinical_action']}")

    # Test CYP2D6 *4/*4 (Tamoxifen resistance)
    cyp_res = PharmacogenomicsEngine.evaluate_cyp2d6_tamoxifen('*4/*4')
    print(f"\n[CYP2D6 Tamoxifen Check: {cyp_res['diplotype']}]")
    print(f"  Phenotype: {cyp_res['phenotype']} -> {cyp_res['efficacy_impact']}")
    print(f"  Clinical Action: {cyp_res['clinical_action']}")
