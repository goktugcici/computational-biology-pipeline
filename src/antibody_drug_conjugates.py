"""
Antibody-Drug Conjugate (ADC) Engineering & Linker-Payload Studio Module
Comprehensive computational biology, pharmacology, and oncology engine for:
- Canonical clinically approved & late-phase ADCs (T-DXd Enhertu, SG Trodelvy, EV Padcev, T-DM1 Kadcyla, BV Adcetris, Dato-DXd).
- Drug-to-Antibody Ratio (DAR) optimization and aggregation kinetics.
- Cleavable (GGFG, CL2A, vc) vs. non-cleavable (SMCC thioether) linker mechanics.
- Bystander killing effect quantification (membrane permeability of free payload).
- Internalization kinetics, lysosomal trafficking, and target antigen density impact.
- Toxicity risk stratification (Interstitial Lung Disease - ILD/Pneumonitis, Peripheral Neuropathy, Neutropenia, Stomatitis).
- NCCN and ASCO clinical practice guidelines.
"""

from typing import Dict, Any, List, Optional


CANONICAL_ADC_CONSTRUCTS: Dict[str, Dict[str, Any]] = {
    "HER2_T_DXD_ENHERTU": {
        "id": "HER2_T_DXD_ENHERTU",
        "name": "Trastuzumab Deruxtecan (T-DXd / Enhertu)",
        "target_antigen": "HER2 (ERBB2 / Neu)",
        "indication": "HER2-Positive & HER2-Low Metastatic Breast Cancer, NSCLC (HER2 mutant), Gastric / GEJ Adenocarcinoma",
        "antibody_type": "Humanized IgG1 mAb (Trastuzumab backbone)",
        "linker_type": "Enzymatically Cleavable Maleimide Tetrapeptide (GGFG: Gly-Gly-Phe-Gly)",
        "linker_cleavage_mechanism": "Lysosomal Cathepsin B and L cleavage following receptor-mediated endocytosis",
        "payload": "DXd (Deruxtecan - Exatecan derivative, potent Topoisomerase I inhibitor)",
        "payload_potency": "10x higher potency than SN-38; sub-nanomolar IC50 (~1 nM)",
        "canonical_dar": 8.0,
        "membrane_permeable_payload": True,
        "bystander_effect": "High (DXd is lipophilic and crosses membrane to kill adjacent antigen-low/negative bystander tumor cells)",
        "plasma_half_life_days": 5.7,
        "clinical_efficacy": "DESTINY-Breast03: 75.8% 12-month PFS (HR=0.28 vs T-DM1). DESTINY-Breast04: First therapy to prolong PFS (9.9 vs 5.1 mos) and OS (23.4 vs 16.8 mos, HR=0.64) in HER2-low (IHC 1+ or IHC 2+/ISH-) breast cancer.",
        "ild_pneumonitis_risk_percent": 12.5,
        "neuropathy_risk_percent": 6.0,
        "neutropenia_risk_percent": 28.0,
        "nccn_guideline": "NCCN Category 1 Preferred 2nd-line for HER2+ metastatic breast cancer after taxane+trastuzumab+pertuzumab. Preferred option for metastatic HER2-low breast cancer post-chemotherapy. Mandatory baseline and serial pulmonary CT monitoring for ILD/Pneumonitis (hold immediately for Grade >=1).",
        "base_potency_score": 96.0,
        "base_bystander_score": 92.0,
        "base_stability_score": 88.0
    },
    "TROP2_SG_TRODELVY": {
        "id": "TROP2_SG_TRODELVY",
        "name": "Sacituzumab Govitecan (Trodelvy)",
        "target_antigen": "Trop-2 (Trophoblast cell-surface antigen 2 / TACSTD2)",
        "indication": "Metastatic Triple-Negative Breast Cancer (mTNBC) & HR+/HER2- Breast Cancer (pretreated), Urothelial Carcinoma",
        "antibody_type": "Humanized IgG1 kappa mAb (hRS7 backbone)",
        "linker_type": "Hydrolyzable pH-Sensitive CL2A Linker",
        "linker_cleavage_mechanism": "Dual cleavage: Intracellular acidic lysosome and extracellular mildly acidic tumor microenvironment",
        "payload": "SN-38 (Active metabolite of Irinotecan, Topoisomerase I inhibitor)",
        "payload_potency": "1,000-fold more active than irinotecan (IC50 ~1-10 nM)",
        "canonical_dar": 7.6,
        "membrane_permeable_payload": True,
        "bystander_effect": "High (Free SN-38 readily penetrates neighboring tumor cells in hypoxic Trop-2 heterogeneous tumors)",
        "plasma_half_life_days": 4.2,
        "clinical_efficacy": "ASCENT Phase 3 Trial (mTNBC): Prolonged median PFS (5.6 vs 1.7 months, HR=0.41) and median OS (12.1 vs 6.7 months, HR=0.48) vs single-agent chemotherapy. TROPiCS-02: Significant PFS and OS benefit in endocrine-resistant HR+/HER2- breast cancer.",
        "ild_pneumonitis_risk_percent": 1.5,
        "neuropathy_risk_percent": 4.0,
        "neutropenia_risk_percent": 58.0,
        "nccn_guideline": "NCCN Category 1 for metastatic TNBC who received >=2 prior systemic therapies (at least one in metastatic setting). FDA Black Box Warning for severe neutropenia (G-CSF support recommended) and severe diarrhea (screen for UGT1A1*28 poor metabolizers).",
        "base_potency_score": 90.0,
        "base_bystander_score": 94.0,
        "base_stability_score": 76.0
    },
    "NECTIN4_EV_PADCEV": {
        "id": "NECTIN4_EV_PADCEV",
        "name": "Enfortumab Vedotin (Padcev)",
        "target_antigen": "Nectin-4 (PVRL4 / Cell adhesion molecule)",
        "indication": "Locally Advanced or Metastatic Urothelial Carcinoma (Bladder & Urinary Tract Cancer)",
        "antibody_type": "Fully Human IgG1 kappa mAb (AGS-22C3 backbone)",
        "linker_type": "Protease-Cleavable Valine-Citrulline (vc) Dipeptide Linker",
        "linker_cleavage_mechanism": "Lysosomal Cathepsin B-mediated cleavage upon cellular internalization",
        "payload": "MMAE (Monomethyl Auristatin E - Synthetic Auristatin Microtubule Disruptor)",
        "payload_potency": "Antimitotic agent blocking tubulin polymerization (sub-nanomolar IC50 ~0.2 nM)",
        "canonical_dar": 3.8,
        "membrane_permeable_payload": True,
        "bystander_effect": "High (Neutral uncharged free MMAE freely diffuses across lipid bilayer membranes)",
        "plasma_half_life_days": 3.6,
        "clinical_efficacy": "EV-302 / KEYNOTE-A39 Trial: Combination of Enfortumab Vedotin + Pembrolizumab doubled median overall survival (31.5 vs 16.1 months, HR=0.47) and median PFS (12.5 vs 6.3 months) compared to platinum-based chemotherapy.",
        "ild_pneumonitis_risk_percent": 4.0,
        "neuropathy_risk_percent": 52.0,
        "neutropenia_risk_percent": 18.0,
        "nccn_guideline": "NCCN Preferred 1st-Line Regimen (Category 1) in combination with Pembrolizumab for patients with locally advanced or metastatic urothelial carcinoma regardless of cisplatin eligibility. Warning for severe cutaneous adverse reactions (SJS/TEN), peripheral neuropathy, and hyperglycemia.",
        "base_potency_score": 94.0,
        "base_bystander_score": 88.0,
        "base_stability_score": 84.0
    },
    "HER2_T_DM1_KADCYLA": {
        "id": "HER2_T_DM1_KADCYLA",
        "name": "Trastuzumab Emtansine (T-DM1 / Kadcyla)",
        "target_antigen": "HER2 (ERBB2 / Neu)",
        "indication": "HER2-Positive Early Breast Cancer with Residual Invasive Disease post-neoadjuvant; previously treated metastatic HER2+ breast cancer",
        "antibody_type": "Humanized IgG1 mAb (Trastuzumab)",
        "linker_type": "Non-Cleavable SMCC (Succinimidyl 4-(N-maleimidomethyl)cyclohexane-1-carboxylate) Thioether Linker",
        "linker_cleavage_mechanism": "Requires complete lysosomal proteolytic degradation of the mAb yielding Lys-SMCC-DM1",
        "payload": "DM1 (Maytansinoid derivative - Microtubule polymerization inhibitor)",
        "payload_potency": "Microtubule assembly arrest (IC50 ~0.1 nM)",
        "canonical_dar": 3.5,
        "membrane_permeable_payload": False,
        "bystander_effect": "Minimal / None (Lys-SMCC-DM1 metabolite retains a charged lysine residue, preventing passive cell membrane diffusion)",
        "plasma_half_life_days": 4.0,
        "clinical_efficacy": "KATHERINE Phase 3 Trial: Adjuvant T-DM1 in patients with residual invasive HER2+ disease reduced risk of recurrence or death by 50% (IDFS HR=0.50, p<0.001) compared to adjuvant trastuzumab alone.",
        "ild_pneumonitis_risk_percent": 1.2,
        "neuropathy_risk_percent": 12.0,
        "neutropenia_risk_percent": 14.0,
        "nccn_guideline": "NCCN Standard of Care (Category 1) for adjuvant treatment of HER2+ early breast cancer patients with residual invasive disease after completing neoadjuvant taxane-based chemotherapy plus anti-HER2 therapy. Monitor liver function tests (AST/ALT) and platelet counts (thrombocytopenia risk).",
        "base_potency_score": 85.0,
        "base_bystander_score": 15.0,
        "base_stability_score": 98.0
    },
    "CD30_BV_ADCETRIS": {
        "id": "CD30_BV_ADCETRIS",
        "name": "Brentuximab Vedotin (Adcetris)",
        "target_antigen": "CD30 (TNFRSF8 / Hodgkin Antigen)",
        "indication": "Classical Hodgkin Lymphoma (cHL), Systemic Anaplastic Large Cell Lymphoma (sALCL), CD30+ Cutaneous T-Cell Lymphoma",
        "antibody_type": "Chimeric IgG1 mAb (cAC10 backbone)",
        "linker_type": "Protease-Cleavable mc-vc-PABC (Maleimidocaproyl-valine-citrulline-p-aminobenzyloxycarbonyl)",
        "linker_cleavage_mechanism": "Lysosomal Cathepsin B selective cleavage with self-immolative PABC spacer release",
        "payload": "MMAE (Monomethyl Auristatin E - Microtubule Disruptor)",
        "payload_potency": "Sub-nanomolar tubulin disruption (IC50 ~0.2 nM)",
        "canonical_dar": 4.0,
        "membrane_permeable_payload": True,
        "bystander_effect": "Moderate-High (Diffuses into Reed-Sternberg microenvironment cells)",
        "plasma_half_life_days": 6.5,
        "clinical_efficacy": "ECHELON-1 Phase 3 Trial: A+AVD (Brentuximab + Doxorubicin, Vinblastine, Dacarbazine) demonstrated significant overall survival benefit (HR=0.59) over standard ABVD, eliminating bleomycin-induced pulmonary toxicity.",
        "ild_pneumonitis_risk_percent": 2.5,
        "neuropathy_risk_percent": 67.0,
        "neutropenia_risk_percent": 54.0,
        "nccn_guideline": "NCCN Category 1 preferred frontline therapy in combination with AVD for adult patients with previously untreated Stage III or IV classical Hodgkin lymphoma. Also Category 1 preferred for relapsed/refractory cHL prior to autologous stem cell transplant.",
        "base_potency_score": 92.0,
        "base_bystander_score": 78.0,
        "base_stability_score": 86.0
    },
    "TROP2_DATO_DXD": {
        "id": "TROP2_DATO_DXD",
        "name": "Datopotamab Deruxtecan (Dato-DXd)",
        "target_antigen": "Trop-2 (TACSTD2)",
        "indication": "Pretreated Advanced Non-Small Cell Lung Cancer (NSCLC) & HR+/HER2- Metastatic Breast Cancer",
        "antibody_type": "Humanized IgG1 anti-Trop-2 mAb",
        "linker_type": "Enzymatically Cleavable Tetrapeptide Linker (GGFG)",
        "linker_cleavage_mechanism": "Lysosomal Cathepsins B and L cleavage with stable plasma circulation",
        "payload": "DXd (Deruxtecan Topoisomerase I inhibitor)",
        "payload_potency": "High potency topoisomerase I inhibitor (IC50 ~1 nM)",
        "canonical_dar": 4.0,
        "membrane_permeable_payload": True,
        "bystander_effect": "High (Permeable DXd diffuses into antigen-negative tumor cells)",
        "plasma_half_life_days": 6.0,
        "clinical_efficacy": "TROPION-Lung01: Significant PFS benefit in non-squamous NSCLC (HR=0.63 in clinically relevant subgroups). TROPION-Breast01: Statistically significant PFS improvement (HR=0.63, median PFS 6.9 vs 4.9 months) in pretreated HR+/HER2- breast cancer.",
        "ild_pneumonitis_risk_percent": 4.5,
        "neuropathy_risk_percent": 3.0,
        "neutropenia_risk_percent": 11.0,
        "nccn_guideline": "Investigational / Emerging therapeutic option for endocrine-resistant metastatic HR+/HER2- breast cancer and pretreated EGFR/actionable genomic alteration non-squamous NSCLC. Lower DAR of 4.0 designed to widen therapeutic index, reduce hematologic/pulmonary toxicity, and manage stomatitis with steroid mouthwash.",
        "base_potency_score": 88.0,
        "base_bystander_score": 85.0,
        "base_stability_score": 94.0
    }
}


class AntibodyDrugConjugateEngine:
    """
    Bioinformatics and precision oncology engine for in-silico simulation of
    Antibody-Drug Conjugates (ADCs), linker-payload kinetics, DAR optimization,
    bystander antitumoral diffusion, and clinical toxicity risk stratification.
    """

    def __init__(self):
        self.constructs = CANONICAL_ADC_CONSTRUCTS

    def list_constructs(self) -> List[Dict[str, Any]]:
        """Returns all canonical ADC constructs with complete metadata."""
        return list(self.constructs.values())

    def get_construct(self, construct_id: str) -> Optional[Dict[str, Any]]:
        """Fetches a specific ADC construct by ID."""
        return self.constructs.get(construct_id)

    def evaluate_adc_kinetics_and_safety(
        self,
        construct_id: str,
        custom_dar: Optional[float] = None,
        target_antigen_expression_percent: float = 75.0,
        tumor_vascularization_percent: float = 60.0
    ) -> Dict[str, Any]:
        """
        Simulates pharmacokinetic & pharmacodynamic indices for an ADC construct:
        - Effective Payload Delivery Score (0-100)
        - Bystander Antitumor Index (0-100)
        - Plasma Stability & Clearance Index (0-100)
        - Interstitial Lung Disease (ILD) / Pneumonitis Risk
        - Peripheral Neuropathy Risk
        - Severe Neutropenia Risk
        - Overall Therapeutic Index & Clinical Guidance
        """
        construct = self.constructs.get(construct_id)
        if not construct:
            raise ValueError(f"Unknown ADC construct: '{construct_id}'. Valid IDs: {list(self.constructs.keys())}")

        dar = float(custom_dar) if custom_dar is not None else float(construct["canonical_dar"])
        antigen_density = max(5.0, min(100.0, float(target_antigen_expression_percent)))
        vascularization = max(10.0, min(100.0, float(tumor_vascularization_percent)))

        # 1. Effective Payload Delivery Score
        # Proportional to DAR * Antigen Density * Vascularization / 100
        # High DAR delivers more payload per antibody, but very high DAR (>8) risks faster hepatic clearance
        dar_efficiency = dar if dar <= 8.0 else (8.0 - (dar - 8.0) * 0.5)
        raw_delivery = (construct["base_potency_score"] * 0.50) + (dar_efficiency * 5.0) + (antigen_density * 0.20) + (vascularization * 0.10)
        payload_delivery_score = round(max(10.0, min(100.0, raw_delivery)), 1)

        # 2. Bystander Antitumor Index (0 - 100)
        # Driven by payload membrane permeability and cleavable vs non-cleavable linker
        if construct["membrane_permeable_payload"]:
            bystander_base = construct["base_bystander_score"]
            # Bystander effect scales with payload potency and diffusion
            bystander_calc = bystander_base + ((dar / 8.0) * 10.0) - ((100.0 - antigen_density) * 0.15)
            bystander_index = round(max(15.0, min(100.0, bystander_calc)), 1)
            bystander_summary = "High Bystander Killing (Permeable payload crosses cell membranes to destroy adjacent low-antigen clones, overcoming intratumoral heterogeneity)."
        else:
            bystander_index = construct["base_bystander_score"]
            bystander_summary = "Zero / Negligible Bystander Effect (Non-cleavable linker traps polar catabolite inside target cell; sparing antigen-negative neighboring cells)."

        # 3. Plasma Stability & Clearance Index (0 - 100)
        # Excessive DAR (>6) increases hydrophobicity and accelerated hepatic reticuloendothelial clearance
        dar_penalty = max(0.0, (dar - 4.0) * 2.5) if dar > 4.0 else 0.0
        stability_calc = construct["base_stability_score"] - dar_penalty
        plasma_stability_score = round(max(20.0, min(100.0, stability_calc)), 1)

        # 4. Toxicity Profile Modeling (Adjusted for DAR)
        dar_toxicity_multiplier = dar / float(construct["canonical_dar"])
        ild_risk = round(construct["ild_pneumonitis_risk_percent"] * dar_toxicity_multiplier, 1)
        neuropathy_risk = round(construct["neuropathy_risk_percent"] * dar_toxicity_multiplier, 1)
        neutropenia_risk = round(construct["neutropenia_risk_percent"] * dar_toxicity_multiplier, 1)

        # Determine Primary Toxicity Warning & Management
        if "DXd" in construct["payload"]:
            primary_toxicity = f"Interstitial Lung Disease (ILD) / Pneumonitis (Estimated Risk: %{ild_risk})"
            primary_toxicity_badge = "rose" if ild_risk >= 10.0 else "amber"
            toxicity_action = "Baseline high-resolution chest CT required. Discontinue permanently for Grade >=2 ILD/pneumonitis; initiate high-dose systemic corticosteroid therapy (>=1 mg/kg/day prednisolone equivalent) immediately."
        elif "MMAE" in construct["payload"]:
            primary_toxicity = f"Peripheral Sensory Neuropathy (Estimated Risk: %{neuropathy_risk})"
            primary_toxicity_badge = "rose" if neuropathy_risk >= 40.0 else "amber"
            toxicity_action = "Serial neurologic assessments before each cycle. For Grade 2 neuropathy, withhold dose until improvement to <=Grade 1, then resume at reduced dose. For Grade 3/4, discontinue permanently."
        elif "SN-38" in construct["payload"]:
            primary_toxicity = f"Severe Neutropenia (%{neutropenia_risk}) & Diarrhea"
            primary_toxicity_badge = "rose" if neutropenia_risk >= 50.0 else "amber"
            toxicity_action = "Prophylactic G-CSF support recommended. Test for UGT1A1*28 polymorphism. For severe diarrhea, initiate high-dose loperamide and oral hydration promptly."
        else:
            primary_toxicity = f"Hepatotoxicity & Thrombocytopenia (%{neutropenia_risk})"
            primary_toxicity_badge = "amber"
            toxicity_action = "Monitor serum transaminases (AST/ALT) and platelet counts prior to each dose infusion. Dose reduction required for recurrent Grade 3 thrombocytopenia."

        # Overall Therapeutic Index
        therapeutic_index = round((payload_delivery_score * 0.60) + (plasma_stability_score * 0.40) - (ild_risk * 0.30) - (neuropathy_risk * 0.15), 1)

        return {
            "construct_id": construct["id"],
            "name": construct["name"],
            "target_antigen": construct["target_antigen"],
            "indication": construct["indication"],
            "antibody_type": construct["antibody_type"],
            "linker_type": construct["linker_type"],
            "linker_cleavage_mechanism": construct["linker_cleavage_mechanism"],
            "payload": construct["payload"],
            "payload_potency": construct["payload_potency"],
            "dar": dar,
            "canonical_dar": construct["canonical_dar"],
            "antigen_expression_percent": antigen_density,
            "tumor_vascularization_percent": vascularization,
            "payload_delivery_score": payload_delivery_score,
            "bystander_index": bystander_index,
            "bystander_summary": bystander_summary,
            "plasma_stability_score": plasma_stability_score,
            "plasma_half_life_days": construct["plasma_half_life_days"],
            "therapeutic_index": therapeutic_index,
            "ild_risk_percent": ild_risk,
            "neuropathy_risk_percent": neuropathy_risk,
            "neutropenia_risk_percent": neutropenia_risk,
            "primary_toxicity": primary_toxicity,
            "primary_toxicity_badge": primary_toxicity_badge,
            "toxicity_action": toxicity_action,
            "clinical_efficacy": construct["clinical_efficacy"],
            "nccn_guideline": construct["nccn_guideline"]
        }


def analyze_adc_construct(
    construct_id: str,
    custom_dar: Optional[float] = None,
    target_antigen_expression_percent: float = 75.0,
    tumor_vascularization_percent: float = 60.0
) -> Dict[str, Any]:
    """Helper convenience function for rapid ADC construct evaluation."""
    engine = AntibodyDrugConjugateEngine()
    return engine.evaluate_adc_kinetics_and_safety(
        construct_id,
        custom_dar,
        target_antigen_expression_percent,
        tumor_vascularization_percent
    )
