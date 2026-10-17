"""
Bispecific Antibodies & Multi-Specific T-Cell Engagers (BiTE) Studio Module
Comprehensive computational biophysics, protein engineering, and clinical oncology engine for:
- Canonical FDA/EMA-approved and landmark investigational bispecific constructs:
  Blinatumomab (Blincyto), Amivantamab (Rybrevant), Teclistamab (Tecvayli),
  Tebentafusp (Kimmtrak), Epcoritamab (Tepkinly), Talquetamab (Talvey).
- Molecular architectural formats:
  BiTE (Tandem scFv without Fc, ~55 kDa, ~2h half-life),
  DuoBody / Full IgG (Controlled Fab-arm exchange, Fc-bearing, ~150 kDa, ~10-14d half-life),
  ImmTAC (High-affinity TCR soluble targeting domain x anti-CD3 scFv, intracellular antigen targeting via HLA-A*02:01).
- Immunological synapse distance (<15 nm optimal lysis vs steric shielding).
- Target vs. CD3 affinity ratio (Kd optimization to avoid peripheral T-cell trapping).
- Step-up dosing kinetics & ASTCT Cytokine Release Syndrome (CRS) risk stratification.
- On-target off-tumor toxicity modeling (dysgeusia/nail toxicity for GPRC5D, rash for EGFR/MET).
- Landmark Phase 3 clinical trials (MARIPOSA, MajesTEC-1, MonumenTAL-1, EPCORE, IMCgp100-202).
"""

from typing import Dict, Any, List, Optional


CANONICAL_BISPECIFIC_CONSTRUCTS: Dict[str, Dict[str, Any]] = {
    "BLINATUMOMAB_BLINCYTO": {
        "id": "BLINATUMOMAB_BLINCYTO",
        "name": "Blinatumomab (Blincyto)",
        "target_tumor": "CD19 (B-Cell Surface Antigen)",
        "target_effector": "CD3 epsilon (T-Cell Co-Receptor)",
        "format": "BiTE (Bispecific T-cell Engager - Tandem scFv, No Fc)",
        "molecular_weight_kda": 54.1,
        "half_life_hours": 2.1,  # Requires continuous 28-day IV infusion
        "delivery_mode": "Continuous 28-day IV infusion via portable pump",
        "synapse_distance_nm": 11.2,  # Compact distance for lytic pore formation
        "kd_target_nm": 1.49,
        "kd_cd3_nm": 100.0,
        "affinity_ratio": 67.1,  # 67x higher affinity to CD19 than CD3
        "indication": "Relapsed/Refractory B-cell ALL & Minimal Residual Disease (MRD) Positive B-ALL",
        "clinical_trial": "TOWER Phase 3 Trial: Prolonged median overall survival (7.7 vs 4.0 mos, HR=0.71, p=0.01) and doubled complete remission rates (44% vs 25%). ECOG-ACRIN E1910: 83% 3-year OS in MRD-negative consolidation.",
        "crs_risk_percent": 15.0,
        "icans_risk_percent": 24.0,
        "step_up_dosing": "9 mcg/day continuous IV on Days 1-7, then step-up to 28 mcg/day on Days 8-28. Mandatory dexamethasone premedication to prevent severe CRS.",
        "primary_toxicity": "Neurotoxicity / ICANS (Tremor, encephalopathy, aphasia, seizures) & Cytokine Release Syndrome.",
        "nccn_guideline": "NCCN Category 1 Preferred for adult and pediatric patients with CD19+ B-ALL in first or second complete remission with MRD >= 0.1% (MRD-positive) and relapsed/refractory B-ALL.",
        "base_potency": 94.0
    },
    "AMIVANTAMAB_RYBREVANT": {
        "id": "AMIVANTAMAB_RYBREVANT",
        "name": "Amivantamab (Rybrevant)",
        "target_tumor": "EGFR (ErbB1) & MET (HGF Receptor)",
        "target_effector": "Immune Effector Cells (Fc gamma RIIIa / CD16a via enhanced ADCC & trogocytosis)",
        "format": "DuoBody (Full IgG1 bispecific with low fucose for enhanced ADCC/trogocytosis)",
        "molecular_weight_kda": 148.0,
        "half_life_hours": 271.0,  # ~11.3 days
        "delivery_mode": "Weekly IV infusion for 4 weeks, then biweekly (with Lazertinib or Chemotherapy)",
        "synapse_distance_nm": 14.5,
        "kd_target_nm": 1.4,   # Dual EGFR/MET binding
        "kd_cd3_nm": 0.0,      # Effector Fc-directed, non-CD3
        "affinity_ratio": 95.0,
        "indication": "NSCLC with EGFR Exon 20 Insertion & EGFR-mutated NSCLC post-Osimertinib",
        "clinical_trial": "MARIPOSA Phase 3 Trial: Amivantamab + Lazertinib reduced risk of progression or death by 30% (PFS 23.7 vs 16.6 mos, HR=0.70) vs Osimertinib monotherapy in frontline EGFR-mutated NSCLC. PAPILLON: PFS HR=0.40 in exon 20 insertions.",
        "crs_risk_percent": 6.0,
        "icans_risk_percent": 1.0,
        "step_up_dosing": "Split dose on Cycle 1 Day 1 (350 mg) and Day 2 (1050 mg) to prevent severe infusion-related reactions (IRRs). Antihistamines, antipyretics, and glucocorticoids mandatory.",
        "primary_toxicity": "Dermatologic Toxicity (Acneiform rash, paronychia, pruritus) & Infusion-Related Reactions (IRRs ~65%).",
        "nccn_guideline": "NCCN Category 1 Preferred for EGFR Exon 20 insertion NSCLC, and Preferred Combination with Lazertinib for classical EGFR mutations (Ex19del/L858R) in advanced NSCLC.",
        "base_potency": 91.0
    },
    "TECLISTAMAB_TECVAYLI": {
        "id": "TECLISTAMAB_TECVAYLI",
        "name": "Teclistamab (Tecvayli)",
        "target_tumor": "BCMA (B-Cell Maturation Antigen)",
        "target_effector": "CD3 epsilon (T-Cell Co-Receptor)",
        "format": "DuoBody (Humanized IgG4-PAA bispecific T-cell engager)",
        "molecular_weight_kda": 146.0,
        "half_life_hours": 240.0,  # ~10 days
        "delivery_mode": "Subcutaneous (SC) injection with step-up dosing, then weekly 1.5 mg/kg",
        "synapse_distance_nm": 12.8,
        "kd_target_nm": 0.20,
        "kd_cd3_nm": 18.0,
        "affinity_ratio": 90.0,  # High selectivity for BCMA
        "indication": "Relapsed/Refractory Multiple Myeloma after >=4 prior lines (PI, IMiD, anti-CD38)",
        "clinical_trial": "MajesTEC-1 Pivotal Phase 1/2 Trial: 63.0% Overall Response Rate (ORR), with 39.4% achieving Complete Response (CR) or better. Median duration of response was 18.4 months in penta-refractory patients.",
        "crs_risk_percent": 72.0,  # Mostly Grade 1-2 with step-up dosing
        "icans_risk_percent": 14.5,
        "step_up_dosing": "Day 1: 0.06 mg/kg SC; Day 4: 0.3 mg/kg SC; Day 7: 1.5 mg/kg full treatment dose. Premedication with acetaminophen, diphenhydramine, and dexamethasone.",
        "primary_toxicity": "Cytokine Release Syndrome (CRS - 72%), Severe Infections (Hypogammaglobulinemia requiring IVIG), Neutropenia.",
        "nccn_guideline": "NCCN Category 1 for patients with relapsed or refractory multiple myeloma who have received at least 4 prior lines of therapy including a proteasome inhibitor, an immunomodulatory agent, and an anti-CD38 monoclonal antibody.",
        "base_potency": 93.0
    },
    "TEBENTAFUSP_KIMMTRAK": {
        "id": "TEBENTAFUSP_KIMMTRAK",
        "name": "Tebentafusp (Kimmtrak / IMCgp100)",
        "target_tumor": "gp100 Peptide-HLA-A*02:01 Complex (Intracellular Melanocytic Antigen)",
        "target_effector": "CD3 epsilon (T-Cell Co-Receptor)",
        "format": "ImmTAC (Immune Mobilizing Monoclonal TCRs Against Cancer - Soluble TCR x anti-CD3)",
        "molecular_weight_kda": 77.0,
        "half_life_hours": 7.5,
        "delivery_mode": "Weekly IV infusion with step-up dosing (Day 1: 20 mcg, Day 8: 30 mcg, Day 15+: 68 mcg)",
        "synapse_distance_nm": 9.8,  # Ultra-close TCR-like immunological synapse
        "kd_target_nm": 0.05,  # Picomolar affinity TCR domain (50 pM)
        "kd_cd3_nm": 2.5,
        "affinity_ratio": 50.0,
        "indication": "HLA-A*02:01-Positive Adult Patients with Unresectable or Metastatic Uveal Melanoma",
        "clinical_trial": "IMCgp100-202 Phase 3 Trial: Statistically significant improvement in overall survival (median OS 21.7 vs 16.0 mos, HR=0.51, p<0.001) vs investigator's choice (pembrolizumab, ipilimumab, or dacarbazine) - First therapy ever to improve OS in metastatic uveal melanoma.",
        "crs_risk_percent": 89.0,  # Grade 1-2 common upon first 3 doses
        "icans_risk_percent": 3.0,
        "step_up_dosing": "Day 1: 20 mcg IV; Day 8: 30 mcg IV; Day 15 and weekly thereafter: 68 mcg IV. Post-infusion overnight monitoring required for the first 3 doses.",
        "primary_toxicity": "Cytokine Release Syndrome (flushing, fever, hypotension) & Melanocyte Target-Mediated Cutaneous Toxicity (Rash, erythema, vitiligo, pruritus).",
        "nccn_guideline": "NCCN Category 1 Preferred for HLA-A*02:01-positive metastatic uveal melanoma. HLA-A*02:01 genotyping is mandatory prior to treatment initiation.",
        "base_potency": 95.0
    },
    "EPCORITAMAB_TEPKINLY": {
        "id": "EPCORITAMAB_TEPKINLY",
        "name": "Epcoritamab (Tepkinly / Epkinly)",
        "target_tumor": "CD20 (MS4A1 / Mature B-Cell Antigen)",
        "target_effector": "CD3 epsilon (T-Cell Co-Receptor)",
        "format": "DuoBody (Full-length IgG1-bispecific with Fc-silencing mutations)",
        "molecular_weight_kda": 147.0,
        "half_life_hours": 336.0,  # ~14 days
        "delivery_mode": "Subcutaneous (SC) injection with step-up dosing (0.16 mg, 0.8 mg, then 48 mg)",
        "synapse_distance_nm": 13.0,
        "kd_target_nm": 0.40,
        "kd_cd3_nm": 12.0,
        "affinity_ratio": 30.0,
        "indication": "Relapsed/Refractory Diffuse Large B-Cell Lymphoma (DLBCL) & Follicular Lymphoma",
        "clinical_trial": "EPCORE NHL-1 Phase 1/2 Trial: In heavily pretreated relapsed DLBCL (including CAR-T failures), overall response rate of 63.1% with 38.9% complete metabolic remissions. Median duration of response was 15.6 months.",
        "crs_risk_percent": 51.0,  # Low high-grade CRS due to SC formulation
        "icans_risk_percent": 6.1,
        "step_up_dosing": "Cycle 1 Day 1: 0.16 mg SC; Day 8: 0.8 mg SC; Day 15: 48 mg SC; Day 22: 48 mg SC. Subcutaneous delivery provides slower absorption kinetics, blunting peak IL-6 release.",
        "primary_toxicity": "Cytokine Release Syndrome (51%), Neutropenia, Injection site reactions.",
        "nccn_guideline": "NCCN Category 1 for third-line and subsequent DLBCL, high-grade B-cell lymphoma, and relapsed follicular lymphoma post-2 prior lines of systemic therapy.",
        "base_potency": 92.0
    },
    "TALQUETAMAB_TALVEY": {
        "id": "TALQUETAMAB_TALVEY",
        "name": "Talquetamab (Talvey)",
        "target_tumor": "GPRC5D (G Protein-Coupled Receptor Family C Group 5 Member D)",
        "target_effector": "CD3 epsilon (T-Cell Co-Receptor)",
        "format": "DuoBody (Humanized IgG4-bispecific with modified Fc)",
        "molecular_weight_kda": 146.0,
        "half_life_hours": 216.0,  # ~9 days
        "delivery_mode": "Subcutaneous (SC) injection: 0.4 mg/kg weekly or 0.8 mg/kg every 2 weeks",
        "synapse_distance_nm": 13.5,
        "kd_target_nm": 0.85,
        "kd_cd3_nm": 15.0,
        "affinity_ratio": 17.6,
        "indication": "Relapsed/Refractory Multiple Myeloma after >=4 prior lines (Including BCMA-exposed)",
        "clinical_trial": "MonumenTAL-1 Phase 1/2 Trial: 74.1% Overall Response Rate (ORR) in patients receiving 0.4 mg/kg weekly, and 71.7% in patients receiving 0.8 mg/kg biweekly. Highly active in BCMA-refractory and post-CAR-T patients.",
        "crs_risk_percent": 79.0,
        "icans_risk_percent": 11.0,
        "step_up_dosing": "Weekly Schedule: Day 1: 0.01 mg/kg; Day 4: 0.06 mg/kg; Day 7: 0.4 mg/kg. Biweekly: Day 1: 0.01 mg/kg; Day 4: 0.06 mg/kg; Day 7: 0.3 mg/kg; Day 10: 0.8 mg/kg.",
        "primary_toxicity": "GPRC5D On-Target Cutaneous Toxicity (Dysgeusia / taste loss ~70%, Nail dystrophy ~50%, Palmar-plantar desquamation), CRS, and Weight loss.",
        "nccn_guideline": "NCCN Category 1 for relapsed/refractory multiple myeloma with prior exposure to an immunomodulatory agent, proteasome inhibitor, and anti-CD38 antibody. Ideal non-cross-resistant salvage option following BCMA-directed CAR-T or ADC therapy.",
        "base_potency": 91.0
    }
}


class BispecificAntibodyEngine:
    """
    Bioinformatics, protein structural engineering, and clinical oncology engine for
    Bispecific Antibodies and Multi-Specific T-Cell Engagers (BiTE/DuoBody/ImmTAC).
    """

    def __init__(self):
        self.constructs = CANONICAL_BISPECIFIC_CONSTRUCTS

    def list_constructs(self) -> List[Dict[str, Any]]:
        """Returns all canonical bispecific constructs with complete pharmacological metadata."""
        return list(self.constructs.values())

    def get_construct(self, construct_id: str) -> Optional[Dict[str, Any]]:
        """Fetches a specific bispecific construct by ID."""
        return self.constructs.get(construct_id)

    def evaluate_bispecific_mechanics(
        self,
        construct_id: str,
        custom_target_expression: float = 75.0,
        effector_to_target_ratio: float = 5.0,
        step_up_dosing_adherent: bool = True
    ) -> Dict[str, Any]:
        """
        Simulates immunological synapse kinetics, cytolytic potency, peripheral T-cell trapping,
        and ASTCT Cytokine Release Syndrome (CRS) / ICANS neurotoxicity grading.
        """
        c = self.constructs.get(construct_id)
        if not c:
            raise ValueError(f"Unknown construct ID: '{construct_id}'. Valid IDs: {list(self.constructs.keys())}")

        target_exp = max(5.0, min(100.0, float(custom_target_expression)))
        et_ratio = max(0.5, min(20.0, float(effector_to_target_ratio)))

        # 1. Immunological Synapse Efficiency & Cytolytic Lysis Score
        # Optimal immunological synapse distance is ~10-15 nm. Distances >15 nm reduce perforin/granzyme delivery.
        synapse_dist = c["synapse_distance_nm"]
        if synapse_dist <= 12.0:
            synapse_efficiency = 1.0
            synapse_badge = "Optimal Lytic Synapse (<12 nm)"
        elif synapse_dist <= 15.0:
            synapse_efficiency = 0.90
            synapse_badge = "Functional Cytolytic Synapse (12-15 nm)"
        else:
            synapse_efficiency = 0.70
            synapse_badge = "Sterically Sub-optimal (>15 nm)"

        # Target affinity weighting and E:T ratio scaling
        potency_base = c["base_potency"]
        # E:T ratio saturates around 5:1 to 10:1
        et_factor = min(1.3, 0.5 + (et_ratio / 10.0) * 0.8)
        tumor_lysis_score = round(
            min(100.0, (potency_base * 0.55) + (target_exp * 0.25 * synapse_efficiency) + (et_factor * 20.0)),
            1
        )

        # 2. Peripheral T-Cell Trapping & Distribution
        # If Kd_CD3 is too tight (low Kd), bispecific binds T-cells in circulation/spleen before reaching the tumor.
        # High affinity ratio (Kd_CD3 >> Kd_target) favors tumor homing first.
        if c["kd_cd3_nm"] > 0:
            if c["affinity_ratio"] >= 50.0:
                biodistribution_profile = "Optimal Tumor First: High tumor-to-CD3 affinity ratio enables selective tumor binding before T-cell recruitment, minimizing peripheral sink trapping."
                biodistribution_score = 92.0
            elif c["affinity_ratio"] >= 20.0:
                biodistribution_profile = "Balanced Distribution: Moderate affinity ratio provides rapid crosslinking with mild peripheral sequestration."
                biodistribution_score = 80.0
            else:
                biodistribution_profile = "Sub-optimal / Peripheral Trapping Risk: Lower affinity ratio causes early binding to circulating T-cells, potentially dampening solid tumor penetration."
                biodistribution_score = 65.0
        else:
            # Non-CD3 effector (e.g. Amivantamab dual RTK)
            biodistribution_profile = "Dual Receptor Cross-Inhibition: Binds EGFR and MET on tumor cells simultaneously, inducing receptor internalization and Fc-directed immune killing (ADCC/trogocytosis)."
            biodistribution_score = 95.0

        # 3. CRS and Neurotoxicity (ICANS) Risk Stratification (ASTCT Criteria)
        # Step-up dosing reduces Grade 3/4 CRS by 75%
        raw_crs = c["crs_risk_percent"]
        raw_icans = c["icans_risk_percent"]

        if step_up_dosing_adherent:
            active_crs_risk = round(raw_crs * (0.35 if raw_crs > 30 else 0.70), 1)
            active_icans_risk = round(raw_icans * 0.50, 1)
            crs_mitigation_status = "Controlled / Low Severe CRS (Step-up dosing protocol active)"
            crs_badge = "emerald"
        else:
            active_crs_risk = round(raw_crs * 1.5, 1)
            active_icans_risk = round(raw_icans * 1.3, 1)
            crs_mitigation_status = "HIGH RISK OF GRADE >=3 CRS: Step-up dosing omitted!"
            crs_badge = "rose"

        # 4. Format Pharmacokinetics
        is_bite = "BiTE" in c["format"]
        is_immtac = "ImmTAC" in c["format"]
        if is_bite:
            pk_profile = "Rapid renal clearance due to lack of Fc domain (~54 kDa, t1/2 ~2h). Requires continuous 28-day ambulatory IV infusion via programmable pump to maintain steady-state serum levels."
        elif is_immtac:
            pk_profile = "Engineered soluble TCR format with 50 pM affinity for peptide-HLA complex (t1/2 ~7.5h). Weekly IV infusion with close initial hemodynamic monitoring."
        else:
            pk_profile = f"Full IgG-like DuoBody architecture containing engineered Fc domain (t1/2 ~{round(c['half_life_hours']/24.0, 1)} days). Allows convenient weekly or biweekly subcutaneous (SC) or IV administration."

        return {
            "construct_id": c["id"],
            "name": c["name"],
            "target_tumor": c["target_tumor"],
            "target_effector": c["target_effector"],
            "format": c["format"],
            "molecular_weight_kda": c["molecular_weight_kda"],
            "half_life_hours": c["half_life_hours"],
            "delivery_mode": c["delivery_mode"],
            "synapse_distance_nm": c["synapse_distance_nm"],
            "synapse_badge": synapse_badge,
            "tumor_antigen_expression": target_exp,
            "effector_to_target_ratio": et_ratio,
            "tumor_lysis_score": tumor_lysis_score,
            "biodistribution_score": biodistribution_score,
            "biodistribution_profile": biodistribution_profile,
            "step_up_dosing_adherent": step_up_dosing_adherent,
            "crs_risk_percent": active_crs_risk,
            "icans_risk_percent": active_icans_risk,
            "crs_mitigation_status": crs_mitigation_status,
            "crs_badge": crs_badge,
            "step_up_dosing_schedule": c["step_up_dosing"],
            "primary_toxicity": c["primary_toxicity"],
            "pk_profile": pk_profile,
            "clinical_trial": c["clinical_trial"],
            "nccn_guideline": c["nccn_guideline"]
        }


def analyze_bispecific_construct(
    construct_id: str,
    custom_target_expression: float = 75.0,
    effector_to_target_ratio: float = 5.0,
    step_up_dosing_adherent: bool = True
) -> Dict[str, Any]:
    """Helper convenience function for rapid bispecific antibody simulation."""
    engine = BispecificAntibodyEngine()
    return engine.evaluate_bispecific_mechanics(
        construct_id,
        custom_target_expression,
        effector_to_target_ratio,
        step_up_dosing_adherent
    )
