"""
In Silico Drug Screening and Lipinski Drug-Likeness Scoring Engine.
Evaluates Lipinski Rule of 5, Veber criteria, quantitative drug-likeness (QED proxy),
and binding affinity proxy scores against disease-associated target binding pockets.
"""

from typing import Dict, List, Any, Optional

CANDIDATE_DRUG_LIBRARY: List[Dict[str, Any]] = [
    {
        'id': 'OSIMERTINIB',
        'name': 'Osimertinib (Tagrisso)',
        'target_gene': 'EGFR',
        'primary_indication': 'EGFR T790M / L858R NSCLC',
        'molecular_weight': 499.6,
        'xlogp': 3.7,
        'hbd_count': 1,
        'hba_count': 7,
        'rotatable_bonds': 9,
        'tpsa': 75.3,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'GEFITINIB',
        'name': 'Gefitinib (Iressa)',
        'target_gene': 'EGFR',
        'primary_indication': 'EGFR-mutant metastatic NSCLC',
        'molecular_weight': 446.9,
        'xlogp': 4.1,
        'hbd_count': 1,
        'hba_count': 6,
        'rotatable_bonds': 7,
        'tpsa': 68.7,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'SOTORASIB',
        'name': 'Sotorasib (Lumakras)',
        'target_gene': 'KRAS',
        'primary_indication': 'KRAS G12C Advanced NSCLC',
        'molecular_weight': 560.6,
        'xlogp': 2.8,
        'hbd_count': 0,
        'hba_count': 7,
        'rotatable_bonds': 5,
        'tpsa': 91.5,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'ADAGRASIB',
        'name': 'Adagrasib (Krazati)',
        'target_gene': 'KRAS',
        'primary_indication': 'KRAS G12C Colorectal / NSCLC',
        'molecular_weight': 604.1,
        'xlogp': 3.9,
        'hbd_count': 0,
        'hba_count': 6,
        'rotatable_bonds': 5,
        'tpsa': 65.2,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'DABRAFENIB',
        'name': 'Dabrafenib (Tafinlar)',
        'target_gene': 'BRAF',
        'primary_indication': 'BRAF V600E/K Melanoma & NSCLC',
        'molecular_weight': 519.6,
        'xlogp': 4.3,
        'hbd_count': 2,
        'hba_count': 6,
        'rotatable_bonds': 6,
        'tpsa': 105.0,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'VOXELOTOR',
        'name': 'Voxelotor (Oxbryta)',
        'target_gene': 'HBB',
        'primary_indication': 'Sickle Cell Disease (HbS Polymerization)',
        'molecular_weight': 337.4,
        'xlogp': 2.4,
        'hbd_count': 1,
        'hba_count': 5,
        'rotatable_bonds': 5,
        'tpsa': 62.1,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'IMATINIB',
        'name': 'Imatinib (Gleevec)',
        'target_gene': 'BCR-ABL / KIT',
        'primary_indication': 'Chronic Myeloid Leukemia (CML)',
        'molecular_weight': 493.6,
        'xlogp': 3.5,
        'hbd_count': 2,
        'hba_count': 7,
        'rotatable_bonds': 7,
        'tpsa': 86.3,
        'charge': 0,
        'fda_status': 'Approved'
    },
    {
        'id': 'METFORMIN',
        'name': 'Metformin',
        'target_gene': 'AMPK / PRKAA1',
        'primary_indication': 'Type 2 Diabetes Mellitus',
        'molecular_weight': 129.16,
        'xlogp': -1.4,
        'hbd_count': 3,
        'hba_count': 2,
        'rotatable_bonds': 0,
        'tpsa': 88.0,
        'charge': 1,
        'fda_status': 'Approved'
    },
    {
        'id': 'ASPIRIN',
        'name': 'Acetylsalicylic Acid (Aspirin)',
        'target_gene': 'PTGS1 / COX-1',
        'primary_indication': 'Analgesic, Antiplatelet',
        'molecular_weight': 180.16,
        'xlogp': 1.2,
        'hbd_count': 1,
        'hba_count': 3,
        'rotatable_bonds': 3,
        'tpsa': 63.6,
        'charge': 0,
        'fda_status': 'Approved'
    }
]

# Target binding pocket descriptors for affinity scoring
TARGET_POCKET_PROFILES: Dict[str, Dict[str, Any]] = {
    'EGFR': {
        'name': 'EGFR ATP-Binding Kinase Cleft (C797 / L858)',
        'optimal_mw': 475.0,
        'optimal_logp': 3.8,
        'preferred_hba': 7,
        'preferred_hbd': 1,
        'pocket_hydrophobicity': 0.75
    },
    'KRAS': {
        'name': 'KRAS Switch-II Allosteric Pocket (C12 / P-loop)',
        'optimal_mw': 570.0,
        'optimal_logp': 3.2,
        'preferred_hba': 7,
        'preferred_hbd': 0,
        'pocket_hydrophobicity': 0.70
    },
    'BRAF': {
        'name': 'BRAF Kinase V600E Activation Loop Cleft',
        'optimal_mw': 520.0,
        'optimal_logp': 4.2,
        'preferred_hba': 6,
        'preferred_hbd': 2,
        'pocket_hydrophobicity': 0.80
    },
    'HBB': {
        'name': 'Hemoglobin Tetramer Alpha-Cleft N-terminus',
        'optimal_mw': 350.0,
        'optimal_logp': 2.5,
        'preferred_hba': 5,
        'preferred_hbd': 1,
        'pocket_hydrophobicity': 0.55
    }
}


class DockingScorer:
    """
    Evaluates Lipinski rule compliance and computes in silico binding scores
    for pharmacology candidates against specified target pockets.
    """

    @staticmethod
    def evaluate_lipinski(drug: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculates compliance with Lipinski's Rule of 5 and Veber's extensions:
        1. MW <= 500 Da
        2. LogP <= 5.0
        3. H-Bond Donors <= 5
        4. H-Bond Acceptors <= 10
        Veber: Rotatable bonds <= 10, TPSA <= 140
        """
        mw = drug.get('molecular_weight', 0.0)
        logp = drug.get('xlogp', 0.0)
        hbd = drug.get('hbd_count', 0)
        hba = drug.get('hba_count', 0)
        rot_bonds = drug.get('rotatable_bonds', 0)
        tpsa = drug.get('tpsa', 0.0)

        violations = []
        if mw > 500:
            violations.append(f"MW > 500 ({mw} Da)")
        if logp > 5:
            violations.append(f"LogP > 5 ({logp})")
        if hbd > 5:
            violations.append(f"HBD > 5 ({hbd})")
        if hba > 10:
            violations.append(f"HBA > 10 ({hba})")

        veber_violations = []
        if rot_bonds > 10:
            veber_violations.append(f"Rotatable Bonds > 10 ({rot_bonds})")
        if tpsa > 140:
            veber_violations.append(f"TPSA > 140 ({tpsa} Å²)")

        lipinski_pass = len(violations) <= 1  # 1 violation is widely accepted in oral drug discovery
        
        # Calculate Drug-Likeness QED proxy score (0 to 100)
        score_mw = max(0, 100 - abs(mw - 380) * 0.25)
        score_logp = max(0, 100 - abs(logp - 2.5) * 18.0)
        score_hba = max(0, 100 - abs(hba - 5) * 12.0)
        score_hbd = max(0, 100 - abs(hbd - 2) * 15.0)
        score_tpsa = max(0, 100 - abs(tpsa - 75) * 0.6)
        qed_score = round((score_mw * 0.25 + score_logp * 0.25 + score_hba * 0.20 + score_hbd * 0.15 + score_tpsa * 0.15), 1)
        qed_score = min(100.0, max(10.0, qed_score))

        return {
            'drug_id': drug.get('id'),
            'name': drug.get('name'),
            'molecular_weight': mw,
            'xlogp': logp,
            'hbd_count': hbd,
            'hba_count': hba,
            'rotatable_bonds': rot_bonds,
            'tpsa': tpsa,
            'lipinski_violations': violations,
            'lipinski_violation_count': len(violations),
            'lipinski_pass': lipinski_pass,
            'veber_violations': veber_violations,
            'veber_pass': len(veber_violations) == 0,
            'drug_likeness_score': qed_score
        }

    @classmethod
    def calculate_binding_affinity_proxy(cls, drug: Dict[str, Any], target_gene: str) -> Dict[str, Any]:
        """
        Calculates an in silico binding affinity proxy score (kcal/mol proxy, lower/more negative is stronger).
        Uses shape matching, lipophilicity complimentarity, and hydrogen bonding interaction terms.
        """
        target = target_gene.upper()
        if target not in TARGET_POCKET_PROFILES:
            target = 'EGFR'  # default fallback
        
        pocket = TARGET_POCKET_PROFILES[target]
        lip = cls.evaluate_lipinski(drug)

        # Baseline pocket electrostatic & hydrophobic contact energy
        mw_ratio = min(drug['molecular_weight'], pocket['optimal_mw']) / max(drug['molecular_weight'], pocket['optimal_mw'])
        logp_diff = abs(drug['xlogp'] - pocket['optimal_logp'])
        logp_fitness = max(0.0, 1.0 - (logp_diff / 5.0))
        
        hba_fitness = 1.0 - (abs(drug['hba_count'] - pocket['preferred_hba']) / 10.0)
        hbd_fitness = 1.0 - (abs(drug['hbd_count'] - pocket['preferred_hbd']) / 6.0)

        contact_score = (mw_ratio * 0.40) + (logp_fitness * 0.35) + (max(0, hba_fitness) * 0.15) + (max(0, hbd_fitness) * 0.10)
        
        # Converted to simulated estimated ΔG (binding free energy, typically -6 to -12 kcal/mol for approved kinase inhibitors)
        simulated_delta_g = round(-4.5 - (contact_score * 6.8), 2)
        
        # Specific high-affinity boost if known targeted pair
        if drug.get('target_gene') == target:
            simulated_delta_g = round(simulated_delta_g - 1.4, 2)

        return {
            'target_gene': target,
            'target_pocket': pocket['name'],
            'drug_id': drug['id'],
            'drug_name': drug['name'],
            'estimated_delta_g_kcal_mol': simulated_delta_g,
            'contact_fitness_score': round(contact_score * 100, 1),
            'drug_likeness_score': lip['drug_likeness_score'],
            'lipinski_pass': lip['lipinski_pass']
        }

    @classmethod
    def screen_library_against_target(cls, target_gene: str) -> List[Dict[str, Any]]:
        """
        Screens all candidates in CANDIDATE_DRUG_LIBRARY against the target binding pocket
        and returns ranked candidates ordered by binding affinity (most negative ΔG first).
        """
        screen_results = []
        for drug in CANDIDATE_DRUG_LIBRARY:
            res = cls.calculate_binding_affinity_proxy(drug, target_gene)
            screen_results.append(res)
        
        # Rank by binding affinity (lower kcal/mol is more stable)
        screen_results.sort(key=lambda x: x['estimated_delta_g_kcal_mol'])
        for rank, item in enumerate(screen_results, start=1):
            item['screen_rank'] = rank
        return screen_results


if __name__ == '__main__':
    scorer = DockingScorer()
    ranked = scorer.screen_library_against_target('EGFR')
    print("--- In Silico Drug Screening Results (Target: EGFR) ---")
    for r in ranked[:4]:
        print(f"Rank #{r['screen_rank']}: {r['drug_name']} -> dG: {r['estimated_delta_g_kcal_mol']} kcal/mol (Fitness: {r['contact_fitness_score']}%)")
