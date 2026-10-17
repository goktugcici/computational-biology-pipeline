"""
Cell Therapy, CAR-T & Adoptive Cell Engineering Studio Module
Comprehensive computational biology and cellular immunotherapy design engine for:
- 1st to 4th generation CAR-T constructs (CD28 vs. 4-1BB costimulation, armored TRUCKs, dnTGFBR2, IL-15).
- Canonical clinically approved & solid tumor targets (CD19, BCMA, Claudin 18.2, HER2, GPC3).
- ASTCT consensus grading and management for Cytokine Release Syndrome (CRS) and ICANS neurotoxicity.
- In-silico modeling of CAR-T expansion kinetics, exhaustion risk, and solid tumor immunosuppression penetration.
"""

from typing import Dict, Any, List, Optional


CANONICAL_CART_CONSTRUCTS: Dict[str, Dict[str, Any]] = {
    "CD19_B_ALL_KYMRIAH": {
        "id": "CD19_B_ALL_KYMRIAH",
        "name": "Tisagenlecleucel (Kymriah)",
        "target_antigen": "CD19",
        "indication": "Relapsed/Refractory B-Cell Acute Lymphoblastic Leukemia (B-ALL) & DLBCL",
        "binder_type": "Murine FMC63 Single-Chain Variable Fragment (scFv)",
        "hinge_transmembrane": "CD8alpha Hinge and Transmembrane Domain",
        "costimulatory_domain": "4-1BB (CD137)",
        "activation_domain": "CD3zeta (3 ITAMs)",
        "generation": "2nd Generation CAR",
        "armoring_features": "Standard 2nd Generation (High Mitochondrial Biogenesis, Long-term Persistence)",
        "clinical_efficacy": "ELIANA Trial: 82% Overall Remission Rate (CR/CRi) in pediatric/young adult B-ALL. Long-term functional CAR persistence >5 years.",
        "crs_profile": "ASTCT Grade 2-3 common (~77% all grade, 46% Grade >=3). Onset median day 3, duration ~8 days.",
        "icans_profile": "ICANS incidence ~40% (13% Grade >=3). Reversible with high-dose dexamethasone.",
        "nccn_guideline": "NCCN Category 1 for pediatric/young adult patients up to 25 years with B-ALL refractory or in 2nd or later relapse. Prioritize Tocilizumab (8 mg/kg IV) at onset of persistent Grade >=2 CRS.",
        "base_expansion_score": 86.0,
        "base_exhaustion_score": 24.0,
        "base_solid_penetrance": 35.0,
        "target_density_normal_tissues": "B-cell lineage only (Expected B-cell aplasia managed by IVIG replacement)"
    },
    "BCMA_MM_CARVYKTI": {
        "id": "BCMA_MM_CARVYKTI",
        "name": "Ciltacabtagene Autoleucel (Carvykti)",
        "target_antigen": "B-Cell Maturation Antigen (BCMA / TNFRSF17)",
        "indication": "Relapsed/Refractory Multiple Myeloma (Post-Proteasome Inhibitor, IMiD & Anti-CD38)",
        "binder_type": "Bi-Epitope Single-Domain Nanobody (Two heavy-chain-only VHH binders targeting distinct BCMA epitopes)",
        "hinge_transmembrane": "CD8alpha Hinge and Transmembrane",
        "costimulatory_domain": "4-1BB (CD137)",
        "activation_domain": "CD3zeta",
        "generation": "2nd Generation Bivalent Nanobody CAR",
        "armoring_features": "Avidity-driven dual BCMA epitope binding conferring ultra-high avidity even in low BCMA antigen density clones.",
        "clinical_efficacy": "CARTITUDE-4 & CARTITUDE-1: 97.9% Overall Response Rate (ORR), with 82.5% Stringent Complete Response (sCR) and sustained MRD negativity (10^-5).",
        "crs_profile": "Delayed CRS onset (median day 7, range 1-12 days). Grade >=3 CRS in <5% with modern early Tocilizumab intervention.",
        "icans_profile": "Neurotoxicity including delayed parkinsonian-like movement and neurocognitive symptoms (~9%). Requires vigilant neurological monitoring.",
        "nccn_guideline": "FDA Approved & NCCN Preferred for adult patients with relapsed/refractory multiple myeloma after >=1 prior lines of therapy including a PI and an IMiD. REMS safety certification mandated.",
        "base_expansion_score": 94.0,
        "base_exhaustion_score": 18.0,
        "base_solid_penetrance": 42.0,
        "target_density_normal_tissues": "Plasma cells and mature B-cells (Hypogammaglobulinemia expected; IVIG supportive therapy)"
    },
    "CD19_DLBCL_YESCARTA": {
        "id": "CD19_DLBCL_YESCARTA",
        "name": "Axicabtagene Ciloleucel (Yescarta)",
        "target_antigen": "CD19",
        "indication": "Relapsed/Refractory Large B-Cell Lymphoma (DLBCL, PMBCL, Transformed Follicular)",
        "binder_type": "Murine FMC63 scFv",
        "hinge_transmembrane": "CD28 Hinge and Transmembrane",
        "costimulatory_domain": "CD28",
        "activation_domain": "CD3zeta",
        "generation": "2nd Generation CD28-Costimulated CAR",
        "armoring_features": "Rapid effector differentiation, intense early glycolytic flux, high initial peak cytotoxic expansion (higher peak cytokines).",
        "clinical_efficacy": "ZUMA-1 & ZUMA-7: 83% ORR, 58% CR in refractory DLBCL. Statistically superior event-free survival (EFS HR 0.40) over salvage chemo-autologous stem cell transplantation.",
        "crs_profile": "Rapid onset CRS (median day 2). Higher initial peak IL-6 and IFN-gamma release.",
        "icans_profile": "ICANS incidence ~64% (28% Grade >=3). Protocolized early administration of dexamethasone or methylprednisolone.",
        "nccn_guideline": "Category 1 Preferred 2nd-line treatment for primary refractory or early relapsing (<12 months) DLBCL based on ZUMA-7 survival data.",
        "base_expansion_score": 92.0,
        "base_exhaustion_score": 42.0,
        "base_solid_penetrance": 38.0,
        "target_density_normal_tissues": "B-cell lineage exclusive"
    },
    "CLDN18_2_GASTRIC_CT041": {
        "id": "CLDN18_2_GASTRIC_CT041",
        "name": "Satricabtagene Autoleucel (CT041)",
        "target_antigen": "Claudin 18.2 (CLDN18.2)",
        "indication": "Advanced Gastric / Gastroesophageal Junction (GEJ) & Pancreatic Adenocarcinoma",
        "binder_type": "Humanized anti-CLDN18.2 scFv (Selective for exposed tumor tight junction conformers)",
        "hinge_transmembrane": "CD8alpha Hinge and Transmembrane",
        "costimulatory_domain": "4-1BB (CD137)",
        "activation_domain": "CD3zeta",
        "generation": "2nd Generation Solid Tumor CAR",
        "armoring_features": "Selective epitope exposure in malignant epithelial cells where tight junctions are lost, sparing normal gastric mucosal junctions.",
        "clinical_efficacy": "Phase 1b/2 Trial: 57.1% ORR and 75.0% Disease Control Rate (DCR) in heavily pre-treated Claudin 18.2-positive gastric cancer patients.",
        "crs_profile": "Predominantly Grade 1-2 CRS (~95% manageable). Low Grade >=3 CRS.",
        "icans_profile": "Low neurotoxicity risk (<3% ICANS). Well tolerated in solid gastrointestinal malignancies.",
        "nccn_guideline": "Breakthrough Therapy Designation (FDA / EMA): Evaluated for patients with CLDN18.2 expression >= 75% tumor cells with >= 2+ staining intensity. Compatible with Zolbetuximab-naïve and resistant cohorts.",
        "base_expansion_score": 78.0,
        "base_exhaustion_score": 38.0,
        "base_solid_penetrance": 68.0,
        "target_density_normal_tissues": "Gastric differentiated epithelial mucosa (Potential transient gastric mucosal erosion / nausea)"
    },
    "HER2_ARMORED_TRUCK_SOLID": {
        "id": "HER2_ARMORED_TRUCK_SOLID",
        "name": "Armored 4th Generation TRUCK (HER2-dnTGFBR2-IL15)",
        "target_antigen": "HER2 / ERBB2 (Extracellular Domain IV)",
        "indication": "HER2-Amplified Solid Tumors (Glioblastoma, Osteosarcoma & Metastatic Breast Cancer)",
        "binder_type": "Humanized Trastuzumab-Derived scFv (4D5)",
        "hinge_transmembrane": "CD8alpha Hinge and Transmembrane",
        "costimulatory_domain": "4-1BB (CD137)",
        "activation_domain": "CD3zeta",
        "generation": "4th Generation 'TRUCK' (T cells Redirected for Universal Cytokine Killing)",
        "armoring_features": "Dominant-negative TGF-beta Receptor Type II (dnTGFBR2) decoy + Constitutive Transgenic Membrane-Bound IL-15 expression to counter immunosuppressive stroma.",
        "clinical_efficacy": "Preclinical & Early Phase 1: Robust resistance against TGF-beta-mediated immunosuppression, sustained stem-like T memory (Tscm) phenotype, and elimination of heterogenous solid tumors.",
        "crs_profile": "Gradual moderate systemic inflammatory response with IL-15 autocrine loop. Requires IL-6 and ferritin surveillance.",
        "icans_profile": "Monitored for on-target off-tumor lung/cardiac HER2 cross-reactivity and CNS penetration in glioblastoma.",
        "nccn_guideline": "Investigational Phase 1/2 Clinical Trial Regimen: Solid tumor cell engineering paradigm addressing physical CAF desmoplasia and soluble TGF-beta barriers.",
        "base_expansion_score": 88.0,
        "base_exhaustion_score": 20.0,
        "base_solid_penetrance": 89.0,
        "target_density_normal_tissues": "Low baseline HER2 in cardiopulmonary/epithelial tissue (Requires safety hinge and fine affinity tuning)"
    },
    "GPC3_HCC_SYNTHETIC": {
        "id": "GPC3_HCC_SYNTHETIC",
        "name": "Armored 4th Generation GPC3-IL7-CCL19 CAR-T",
        "target_antigen": "Glypican-3 (GPC3)",
        "indication": "Advanced Unresectable Hepatocellular Carcinoma (HCC)",
        "binder_type": "Humanized GC33 scFv",
        "hinge_transmembrane": "CD8alpha Hinge and Transmembrane",
        "costimulatory_domain": "4-1BB (CD137)",
        "activation_domain": "CD3zeta",
        "generation": "4th Generation 7x19 Dual-Chemokine Armed CAR",
        "armoring_features": "Constitutive secretion of Interleukin-7 (IL-7) and Chemokine Ligand 19 (CCL19) to induce host dendritic cell and T-cell infiltration into liver parenchyma.",
        "clinical_efficacy": "Early Clinical Trial: Complete radiological remission and prolonged intra-tumoral persistence in GPC3+ hepatocellular carcinoma refractory to sorafenib/atezolizumab.",
        "crs_profile": "Mild to moderate CRS (Grade 1-2). Hepatic transaminase transient flare due to localized immune activation.",
        "icans_profile": "Negligible neurotoxicity risk.",
        "nccn_guideline": "Investigational Platform: Clinical trial enrollment recommended for GPC3-positive HCC progressing post-checkpoint inhibitor therapy.",
        "base_expansion_score": 82.0,
        "base_exhaustion_score": 26.0,
        "base_solid_penetrance": 84.0,
        "target_density_normal_tissues": "Fetal oncofetal protein; minimal to absent in normal adult human liver"
    }
}


class CellTherapyCartEngine:
    """
    Engine for CAR-T and cellular immunotherapy design, proliferation modeling,
    exhaustion risk quantification, and ASTCT toxicity grading.
    """

    def __init__(self):
        self.constructs = CANONICAL_CART_CONSTRUCTS

    def get_construct(self, construct_id: str) -> Dict[str, Any]:
        """Retrieve canonical CAR-T construct by ID."""
        return self.constructs.get(construct_id, self.constructs["CD19_B_ALL_KYMRIAH"])

    def list_constructs(self) -> List[Dict[str, str]]:
        """List all available canonical CAR-T and adoptive cell therapy models."""
        return [
            {
                "id": c["id"],
                "name": c["name"],
                "target_antigen": c["target_antigen"],
                "indication": c["indication"],
                "generation": c["generation"]
            }
            for c in self.constructs.values()
        ]

    def evaluate_cart_fitness_and_safety(
        self,
        construct_id: str,
        costimulatory_override: Optional[str] = None,
        tumor_antigen_density_percent: float = 85.0,
        tme_immunosuppression_percent: float = 25.0
    ) -> Dict[str, Any]:
        """
        Calculates fitness, proliferation index, exhaustion risk, ASTCT CRS/ICANS grades,
        and solid tumor penetration capacity.
        """
        construct = self.get_construct(construct_id)
        costim = costimulatory_override or construct["costimulatory_domain"]

        # Clamp parameters
        antigen_density = max(5.0, min(100.0, float(tumor_antigen_density_percent)))
        immunosuppression = max(0.0, min(100.0, float(tme_immunosuppression_percent)))

        # Costimulatory kinetics: CD28 drives higher peak, 4-1BB drives longevity and metabolic persistence
        if "CD28" in costim:
            costim_expansion_factor = 1.15
            costim_exhaustion_factor = 1.45
            costim_toxicity_weight = 1.25
        elif "4-1BB" in costim:
            costim_expansion_factor = 1.00
            costim_exhaustion_factor = 0.75
            costim_toxicity_weight = 0.90
        else:
            costim_expansion_factor = 0.85
            costim_exhaustion_factor = 1.10
            costim_toxicity_weight = 0.80

        # Calculate Proliferation & Expansion Score (0 - 100)
        raw_expansion = (construct["base_expansion_score"] * costim_expansion_factor) * (antigen_density / 85.0)
        expansion_score = round(max(10.0, min(100.0, raw_expansion)), 1)

        # Calculate Tonic Signaling & Exhaustion Risk Index (0 - 100)
        raw_exhaustion = (construct["base_exhaustion_score"] * costim_exhaustion_factor) + (immunosuppression * 0.45)
        exhaustion_score = round(max(5.0, min(100.0, raw_exhaustion)), 1)

        # Calculate Solid Tumor Penetrance Index (0 - 100)
        is_armored = "4th Generation" in construct["generation"] or "dnTGFBR2" in construct.get("armoring_features", "")
        if is_armored:
            # Armored constructs actively resist immunosuppression
            penetrance = construct["base_solid_penetrance"] - (immunosuppression * 0.20)
        else:
            penetrance = construct["base_solid_penetrance"] - (immunosuppression * 0.65)
        penetrance_index = round(max(10.0, min(100.0, penetrance)), 1)

        # ASTCT CRS Risk Stratification
        crs_index = (expansion_score * 0.50) + (costim_toxicity_weight * 25.0) + (antigen_density * 0.20)
        if crs_index >= 88.0:
            crs_grade = "Grade 3/4 (Severe - Hypotension refractory to fluids, Hypoxia requiring FiO2 >=40%)"
            crs_badge = "rose"
            crs_management = "Immediate ICU transfer. Tocilizumab 8 mg/kg IV (max 800 mg) + Dexamethasone 20 mg IV q6h. Repeat Tocilizumab in 8h if no improvement. Consider Anakinra (IL-1R antagonist) for refractory HLH/MAS."
        elif crs_index >= 65.0:
            crs_grade = "Grade 2 (Moderate - Hypotension responsive to fluids, Hypoxia requiring low-flow nasal cannula)"
            crs_badge = "amber"
            crs_management = "Administer Tocilizumab 8 mg/kg IV (max 800 mg). If persistent after 24-48 hours, add Dexamethasone 10 mg IV every 12 hours. Monitor vitals q2h."
        else:
            crs_grade = "Grade 1 (Mild - Fever >= 38.0°C without hypotension or hypoxia)"
            crs_badge = "teal"
            crs_management = "Supportive care: Acetaminophen, IV hydration, blood cultures, empiric broad-spectrum antibiotics for neutropenic fever. Tocilizumab reserved for progression."

        # ASTCT ICANS Neurotoxicity Risk Stratification
        icans_index = (expansion_score * 0.40) + (costim_toxicity_weight * 30.0)
        if icans_index >= 82.0:
            icans_grade = "Grade 3/4 (Severe - ICE Score 0-2, Coma, Motor Weakness, Cerebral Edema, or Status Epilepticus)"
            icans_management = "Emergency neurocritical care. High-dose Methylprednisolone 1,000 mg/day IV for 3 days followed by rapid taper. Levetiracetam 1,000 mg IV q12h seizure prophylaxis. Diagnostic urgent brain MRI."
        elif icans_index >= 60.0:
            icans_grade = "Grade 2 (Moderate - ICE Score 3-6, expressive dysphasia, tremors, lethargy)"
            icans_management = "Dexamethasone 10 mg IV every 6 hours. Frequent ICE (Immune Effector Cell Encephalopathy) assessments every 4 hours. Continuous aspiration precautions."
        else:
            icans_grade = "Grade 0/1 (Mild/None - ICE Score 7-10, mild confusion or somnolence)"
            icans_management = "Serial ICE scoring every 8-12 hours. Diagnostic handwriting assessments twice daily. Avoid sedating antihistamines or opioids."

        return {
            "construct_id": construct["id"],
            "name": construct["name"],
            "target_antigen": construct["target_antigen"],
            "indication": construct["indication"],
            "generation": construct["generation"],
            "costimulatory_domain": costim,
            "armoring_features": construct["armoring_features"],
            "antigen_density_percent": antigen_density,
            "immunosuppression_percent": immunosuppression,
            "expansion_score": expansion_score,
            "exhaustion_score": exhaustion_score,
            "solid_penetrance_index": penetrance_index,
            "crs_grade": crs_grade,
            "crs_badge": crs_badge,
            "crs_management": crs_management,
            "icans_grade": icans_grade,
            "icans_management": icans_management,
            "nccn_guideline": construct["nccn_guideline"],
            "clinical_efficacy": construct["clinical_efficacy"],
            "normal_tissue_safety": construct["target_density_normal_tissues"]
        }


def analyze_cart_construct(
    construct_id: str,
    costimulatory_override: Optional[str] = None,
    tumor_antigen_density_percent: float = 85.0,
    tme_immunosuppression_percent: float = 25.0
) -> Dict[str, Any]:
    """Helper convenience function for rapid CAR-T construct evaluation."""
    engine = CellTherapyCartEngine()
    return engine.evaluate_cart_fitness_and_safety(
        construct_id,
        costimulatory_override,
        tumor_antigen_density_percent,
        tme_immunosuppression_percent
    )
