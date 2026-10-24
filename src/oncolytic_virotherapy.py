"""
Oncolytic Virotherapy & Synthetic Immuno-Virology Studio Module (OV Engine)
Comprehensive computational virology, viral kinetic modeling, and immuno-oncology engine for:
- Canonical clinical and approved oncolytic viruses:
  1. Talimogene laherparepvec (T-VEC / Imlygic) - HSV-1 (ICP34.5-/ICP47-/GM-CSF+), OPTiM Phase 3, Melanoma.
  2. Teserpaturev / G47Δ (DELYTACT) - Triple-mutated HSV-1 (ICP34.5-/ICP47-/ICP6-), Glioblastoma (PMDA Approved).
  3. Pelareorep (Reolysin) - Wild-type dsRNA Reovirus targeting activated Ras, BRACELET-1 Phase 2, Breast/PDAC.
  4. CAN-2409 (aglatimagene besadenovec) - Ad5-HSV-tk suicide gene therapy + valacyclovir, PDAC/Prostate.
  5. Cretostimogene grenadavirec (CG0070) - E2F-1 Rb-regulated Ad5 + GM-CSF, BOND-003 Phase 3, Bladder NMIBC.
  6. VCN-01 - Stroma-degrading oncolytic adenovirus encoding human PH20 hyaluronidase, VIRAGE Phase 2b, PDAC.
- Dynamical infection & replication ODE kinetics:
  Uninfected tumor cells (U), infected tumor cells (I), and viral titer (V).
- Immunogenic Cell Death (ICD) biomarkers:
  Calreticulin (CRT) cell-surface translocation, extracellular ATP release,
  HMGB1 release, and Type-I Interferon / CXCL10 chemokine induction.
- Systemic Abscopal Effect Index:
  Local tumor lysis vs. systemic CD8+ T-cell priming against non-injected distant metastases.
- Neutralizing Antibody (NAb) Shielding & Pharmacodynamics:
  Pre-existing anti-vector immunity (anti-HSV-1, anti-Ad5) vs. intratumoral delivery and carrier cell protection.
- NCCN, FDA, EMA & PMDA clinical guidelines and trial benchmarks.
"""

import math
from typing import Dict, Any, List, Optional


CANONICAL_OV_CONSTRUCTS: Dict[str, Dict[str, Any]] = {
    "TVEC_IMLYGIC_MELANOMA": {
        "id": "TVEC_IMLYGIC_MELANOMA",
        "name": "Talimogene laherparepvec (T-VEC / Imlygic - Amgen)",
        "virus_family": "Herpesviridae (HSV-1, JS1 strain)",
        "genome_type": "Enveloped dsDNA (~152 kb)",
        "genetic_engineering": "ICP34.5 deletion (neurovirulence attenuation, selective replication in PKR-deficient tumor cells) + ICP47 deletion (restores MHC-I antigen presentation and early US11 expression) + Human GM-CSF (CSF2) insertion (2 copies)",
        "tumor_indication": "Unresectable Stage IIIB-IVM1a Cutaneous Melanoma",
        "clinical_trial": "OPTiM Phase 3 Trial (JCO 2015): Durable Response Rate (DRR) was significantly higher with T-VEC vs GM-CSF (16.3% vs 2.1%, OR=8.9, p<0.001); Overall Response Rate (ORR) 26.4%; median overall survival 23.3 vs 18.9 months; demonstrated local lysis and distant abscopal responses.",
        "primary_receptor": "Nectin-1 (CD111) & HVEM (TNFRSF14)",
        "burst_size": 1000,
        "lysis_cycle_hours": 24.0,
        "default_dose_pfu": 1.0e8,
        "delivery_route": "Intralesional ultrasound-guided direct injection (up to 4.0 mL)",
        "innate_icd_potency": 92.0,
        "gmcsf_expression": True,
        "nab_seroprevalence_percent": 65.0,  # Pre-existing anti-HSV-1 IgG
        "intratumoral_retention_factor": 0.88,
        "nccn_guideline": "NCCN Category 1 Preferred for recurrent, locally advanced, or in-transit unresectable cutaneous melanoma lesions; recommended as monotherapy or with Pembrolizumab (MASTERKEY-265).",
        "base_efficacy_score": 93.0
    },
    "DELYTACT_G47DELTA_GBM": {
        "id": "DELYTACT_G47DELTA_GBM",
        "name": "Teserpaturev / G47Δ (DELYTACT - Daiichi Sankyo)",
        "virus_family": "Herpesviridae (HSV-1, F strain)",
        "genome_type": "Enveloped dsDNA (~152 kb)",
        "genetic_engineering": "Triple deletion/insertion: ICP34.5 deletion + ICP47 deletion + LacZ gene insertion inactivating ICP6 (UL39 ribonucleotide reductase, strictly limiting replication to rapidly dividing glioblastoma cells)",
        "tumor_indication": "Residual or Recurrent Glioblastoma Multiforme (GBM)",
        "clinical_trial": "Phase 2 Clinical Trial (Nature Medicine 2022 - Tomoki Todo et al., Univ of Tokyo): 1-year overall survival rate was 84.2% (16/19 patients) compared to historical control of ~15% (p<0.0001); median OS from initial G47Δ inoculation was 20.2 months with biopsy-proven CD8+ T-cell infiltration into residual brain tumors.",
        "primary_receptor": "Nectin-1 (CD111) on high-grade glioma cells",
        "burst_size": 650,
        "lysis_cycle_hours": 28.0,
        "default_dose_pfu": 1.0e8,
        "delivery_route": "Stereotactic coordinate-guided intracranial injection into tumor bed",
        "innate_icd_potency": 90.0,
        "gmcsf_expression": False,
        "nab_seroprevalence_percent": 55.0,
        "intratumoral_retention_factor": 0.95,
        "nccn_guideline": "PMDA (Japan) Full Regulatory Approval (2021) for malignant glioma. World's first approved oncolytic virus for brain tumors.",
        "base_efficacy_score": 95.0
    },
    "PELAREOREP_REOLYSIN_RAS": {
        "id": "PELAREOREP_REOLYSIN_RAS",
        "name": "Pelareorep / Reolysin (Oncolytics Biotech)",
        "virus_family": "Reoviridae (Type 3 Dearing strain)",
        "genome_type": "Non-enveloped segmented dsRNA (10 segments, ~23.5 kb)",
        "genetic_engineering": "Non-engineered naturally oncolytic wild-type reovirus; selectively replicates in cells with hyperactivated Ras signaling pathway (KRAS, HRAS, NRAS, or upstream RTK/EGFR/HER2) which blocks PKR phosphorylation",
        "tumor_indication": "HR+/HER2- Metastatic Breast Cancer & Pancreatic Ductal Adenocarcinoma (PDAC)",
        "clinical_trial": "BRACELET-1 Phase 2 Trial (ASCO 2023): Pelareorep + Paclitaxel + Avelumab doubled progression-free survival (PFS HR=0.45) in metastatic HR+/HER2- breast cancer; AWARE-1 trial in early breast cancer confirmed profound intratumoral CD8+ T-cell expansion and neoepitope spreading.",
        "primary_receptor": "Junctional Adhesion Molecule-A (JAM-A) & Sialic Acid",
        "burst_size": 1800,
        "lysis_cycle_hours": 18.0,
        "default_dose_pfu": 3.0e10,
        "delivery_route": "Systemic Intravenous (IV) 60-min infusion on Days 1, 2, 8, 9 Q3W",
        "innate_icd_potency": 86.0,
        "gmcsf_expression": False,
        "nab_seroprevalence_percent": 85.0,  # Natural high population seropositivity
        "intratumoral_retention_factor": 0.45, # IV delivery requires high dosing
        "nccn_guideline": "FDA Fast Track Designation for metastatic breast cancer and pancreatic adenocarcinoma; active Phase 3 registrational development.",
        "base_efficacy_score": 88.0
    },
    "CAN2409_HSV_TK_ADENOVIRUS": {
        "id": "CAN2409_HSV_TK_ADENOVIRUS",
        "name": "CAN-2409 / aglatimagene besadenovec (Candel Therapeutics)",
        "virus_family": "Adenoviridae (Ad5 serotype)",
        "genome_type": "Non-enveloped replication-defective dsDNA (~36 kb)",
        "genetic_engineering": "E1/E3 deletion + Insertion of Herpes Simplex Virus Thymidine Kinase (HSV-tk) cDNA driven by RSV promoter; administered with oral prodrug Valacyclovir/Ganciclovir to generate cytotoxic phosphorylated nucleoside analogs",
        "tumor_indication": "Borderline Resectable Pancreatic Ductal Adenocarcinoma (PDAC) & Localized Prostate Cancer",
        "clinical_trial": "Phase 2 Randomized Trial (ASCO 2024): CAN-2409 + Valacyclovir + FOLFIRINOX chemotherapy achieved median Overall Survival of 28.8 months vs 12.5 months in control (HR=0.42, p=0.038) in borderline resectable PDAC; 81.4% R0 surgical resection rate.",
        "primary_receptor": "Coxsackievirus and Adenovirus Receptor (CAR)",
        "burst_size": 0,  # Suicide gene / metabolic bystander killing rather than replication burst
        "lysis_cycle_hours": 36.0,
        "default_dose_pfu": 1.0e11,
        "delivery_route": "Endoscopic Ultrasound (EUS) guided fine-needle injection directly into the pancreas",
        "innate_icd_potency": 89.0,
        "gmcsf_expression": False,
        "nab_seroprevalence_percent": 70.0,
        "intratumoral_retention_factor": 0.90,
        "nccn_guideline": "FDA Fast Track Designation & Orphan Drug Designation for pancreatic ductal adenocarcinoma and non-small cell lung cancer (NSCLC).",
        "base_efficacy_score": 91.0
    },
    "CRETOSTIMOGENE_CG0070_BLADDER": {
        "id": "CRETOSTIMOGENE_CG0070_BLADDER",
        "name": "Cretostimogene Grenadavirec / CG0070 (CG Oncology)",
        "virus_family": "Adenoviridae (Ad5 serotype)",
        "genome_type": "Non-enveloped conditionally replicative dsDNA (~36 kb)",
        "genetic_engineering": "E1A transcription under human E2F-1 promoter (selectively replicates in cells with defective Retinoblastoma [Rb] pathway, ubiquitous in bladder cancer) + Insertion of human GM-CSF cDNA",
        "tumor_indication": "BCG-Unresponsive High-Risk Non-Muscle Invasive Bladder Cancer (NMIBC) with Carcinoma In Situ (CIS)",
        "clinical_trial": "BOND-003 Phase 3 Trial (AUA / SUO 2024): Complete Response (CR) rate of 75.7% (56/74 evaluable patients); 74.4% of responders maintained durable CR at 12 months; eliminates need for radical cystectomy (bladder removal) with zero Grade >=3 treatment-related adverse events.",
        "primary_receptor": "CAR & Integrin alpha-v-beta-3/5",
        "burst_size": 1100,
        "lysis_cycle_hours": 22.0,
        "default_dose_pfu": 1.0e12,
        "delivery_route": "Intravesical instillation via Foley catheter into bladder with 60-min dwell time",
        "innate_icd_potency": 94.0,
        "gmcsf_expression": True,
        "nab_seroprevalence_percent": 68.0,
        "intratumoral_retention_factor": 0.85,
        "nccn_guideline": "FDA Breakthrough Therapy Designation & Fast Track Designation for BCG-unresponsive NMIBC; emerging frontline bladder-preserving therapeutic standard.",
        "base_efficacy_score": 96.0
    },
    "VCN01_HYALURONIDASE_ADENOVIRUS": {
        "id": "VCN01_HYALURONIDASE_ADENOVIRUS",
        "name": "VCN-01 / Stroma-Degrading Oncolytic Adenovirus (Theriva Biologics)",
        "virus_family": "Adenoviridae (Ad5 with RGDK capsid modification)",
        "genome_type": "Non-enveloped conditionally replicative dsDNA (~36 kb)",
        "genetic_engineering": "E1A delta-24 deletion (Rb-binding pRb-defective selectivity) + Late expression of human soluble PH20 hyaluronidase (sPH20) to enzymatically degrade dense desmoplastic tumor extracellular matrix + RGDK capsid fiber modification bypassing CAR dependence",
        "tumor_indication": "Metastatic Pancreatic Ductal Adenocarcinoma (PDAC) & Retinoblastoma",
        "clinical_trial": "VIRAGE Phase 2b Trial: VCN-01 breaks down dense intratumoral hyaluronic acid (HA), reducing interstitial fluid pressure from >40 mmHg to physiological levels (<10 mmHg), allowing Gemcitabine/nab-Paclitaxel to extravasate and driving CD8+ T cells into the necrotic tumor core.",
        "primary_receptor": "Heparan Sulfate Proteoglycans & Integrins (RGDK-retargeted)",
        "burst_size": 850,
        "lysis_cycle_hours": 26.0,
        "default_dose_pfu": 1.0e13,
        "delivery_route": "Systemic Intravenous (IV) Infusion or Endoscopic Ultrasound (EUS) Intratumoral",
        "innate_icd_potency": 91.0,
        "gmcsf_expression": False,
        "nab_seroprevalence_percent": 60.0,
        "intratumoral_retention_factor": 0.65,
        "nccn_guideline": "FDA Orphan Drug Designation for pancreatic cancer; prominent matrix-depleting oncolytic platform for dense stroma solid tumors.",
        "base_efficacy_score": 90.0
    }
}


def model_oncolytic_viral_kinetics(
    initial_tumor_cells: float = 1.0e8,
    initial_dose_pfu: float = 1.0e8,
    infectivity_rate: float = 1.0e-9,
    burst_size: int = 1000,
    lysis_rate_per_day: float = 1.0,  # ~24 hour lysis cycle
    viral_clearance_per_day: float = 0.5,
    days: int = 14
) -> Dict[str, Any]:
    """
    Simulates oncolytic viral replication and tumor lysis kinetics using coupled ODEs:
    dU/dt = r * U * (1 - (U + I)/K) - beta * U * V
    dI/dt = beta * U * V - delta * I
    dV/dt = b * delta * I - gamma * V - beta * U * V

    Where:
    - U: Uninfected tumor cells
    - I: Infected tumor cells
    - V: Free viral particles (PFU)
    - beta: Infectivity rate (PFU^-1 * day^-1)
    - delta: Lysis rate of infected cells (day^-1)
    - b: Burst size (virions released per lysed cell)
    - gamma: Immune / neutralizing clearance rate (day^-1)
    """
    # Carrying capacity K and tumor proliferation rate r
    K = initial_tumor_cells * 2.0
    r = 0.05  # slow basal tumor growth per day

    dt = 0.25  # 6-hour numerical integration steps
    steps = int(days / dt)

    U = float(initial_tumor_cells)
    I = 0.0
    V = float(initial_dose_pfu)

    u_history = [U]
    i_history = [I]
    v_history = [V]
    time_history = [0.0]

    peak_v = V
    peak_v_day = 0.0
    nadir_u = U
    nadir_u_day = 0.0

    for step in range(1, steps + 1):
        t = step * dt

        # Derivatives
        infection_flux = infectivity_rate * U * V
        tumor_growth = r * U * (1.0 - (U + I) / K)

        dU = (tumor_growth - infection_flux) * dt
        dI = (infection_flux - lysis_rate_per_day * I) * dt
        dV = (burst_size * lysis_rate_per_day * I - viral_clearance_per_day * V - infection_flux) * dt

        U = max(0.0, U + dU)
        I = max(0.0, I + dI)
        V = max(0.0, V + dV)

        if V > peak_v:
            peak_v = V
            peak_v_day = t

        if U < nadir_u:
            nadir_u = U
            nadir_u_day = t

        if step % 4 == 0:  # Record daily points
            time_history.append(round(t, 1))
            u_history.append(round(U, 0))
            i_history.append(round(I, 0))
            v_history.append(round(V, 0))

    tumor_lysis_percent = round(((initial_tumor_cells - nadir_u) / initial_tumor_cells) * 100.0, 1)
    tumor_lysis_percent = max(0.0, min(100.0, tumor_lysis_percent))

    viral_replication_fold = round(peak_v / max(1.0, initial_dose_pfu), 1)

    return {
        "initial_tumor_cells": initial_tumor_cells,
        "nadir_tumor_cells": nadir_u,
        "nadir_day": round(nadir_u_day, 1),
        "tumor_lysis_percent": tumor_lysis_percent,
        "initial_dose_pfu": initial_dose_pfu,
        "peak_viral_titer_pfu": peak_v,
        "peak_viral_day": round(peak_v_day, 1),
        "viral_replication_fold": viral_replication_fold,
        "days_simulated": days,
        "time_history": time_history,
        "uninfected_history": u_history,
        "infected_history": i_history,
        "viral_history": v_history
    }


def calculate_icd_biomarker_score(
    lysis_percent: float,
    virus_family: str,
    gmcsf_present: bool = True,
    interferon_defective: bool = True
) -> Dict[str, Any]:
    """
    Evaluates the molecular hallmarks of Immunogenic Cell Death (ICD):
    1. Calreticulin (CRT) cell surface translocation ('Eat-me' signal for dendritic cells).
    2. Extracellular ATP release via pannexin-1 ('Find-me' signal for antigen-presenting cells).
    3. High-Mobility Group Box 1 (HMGB1) nuclear-to-extracellular release (TLR4 activation).
    4. Type-I Interferon (IFN-alpha/beta) & CXCL9/CXCL10 chemokine gradient induction.
    """
    base_icd = (lysis_percent / 100.0) * 70.0

    # Herpesviridae and Reoviridae produce double-stranded RNA/DNA replication intermediates activating STING and TLR3/9
    if "Herpesviridae" in virus_family:
        dna_sensor_boost = 15.0
    elif "Reoviridae" in virus_family:
        dna_sensor_boost = 18.0  # dsRNA strongly triggers RIG-I / MDA5
    else:
        dna_sensor_boost = 12.0

    gmcsf_boost = 10.0 if gmcsf_present else 0.0

    composite_icd = min(99.0, max(15.0, base_icd + dna_sensor_boost + gmcsf_boost))
    composite_icd = round(composite_icd, 1)

    # Specific biomarker predictions
    crt_translocation_pct = round(min(85.0, composite_icd * 0.82), 1)
    atp_release_uM = round(min(12.5, max(0.5, composite_icd * 0.12)), 2)
    hmgb1_release_ng_ml = round(min(50.0, max(2.0, composite_icd * 0.48)), 1)
    ifn_beta_fold = round(min(65.0, max(1.0, composite_icd * 0.62)), 1)

    if composite_icd >= 80.0:
        microenvironment_status = "Inflamed / 'Hot' Microenvironment: High CRT, HMGB1 & CXCL10 induction"
    elif composite_icd >= 55.0:
        microenvironment_status = "Primed / Intermediate Microenvironment: Moderate immune influx"
    else:
        microenvironment_status = "Immune Desert / 'Cold' Microenvironment: Suboptimal ICD induction"

    return {
        "composite_icd_score": composite_icd,
        "calreticulin_translocation_percent": crt_translocation_pct,
        "extracellular_atp_uM": atp_release_uM,
        "hmgb1_release_ng_ml": hmgb1_release_ng_ml,
        "ifn_beta_induction_fold": ifn_beta_fold,
        "microenvironment_status": microenvironment_status
    }


def calculate_abscopal_systemic_index(
    local_lysis_score: float,
    icd_score: float,
    delivery_route: str,
    anti_pd1_combination: bool = True
) -> Dict[str, Any]:
    """
    Computes the Abscopal Effect Index (systemic anti-tumor immune response
    targeting distant un-injected metastatic lesions):
    - Intralesional injection creates in situ tumor vaccine.
    - Checkpoint blockade (anti-PD-1 / anti-CTLA-4) prevents systemic T-cell exhaustion.
    """
    route_lower = delivery_route.lower()
    route_factor = 1.0 if any(term in route_lower for term in ["intralesional", "intratumoral", "stereotactic", "eus", "intravesical"]) else 0.75

    raw_abscopal = (0.45 * local_lysis_score + 0.55 * icd_score) * route_factor

    if anti_pd1_combination:
        # Checkpoint blockade synergizes supra-additively
        abscopal_index = min(98.0, raw_abscopal * 1.25)
        synergy_label = "Supra-Additive Synergy with Anti-PD-1/PD-L1 Checkpoint Blockade"
    else:
        abscopal_index = min(75.0, raw_abscopal * 0.85)
        synergy_label = "Monotherapy (Susceptible to Adaptive Immune Resistance / PD-L1 Upregulation)"

    abscopal_index = round(abscopal_index, 1)

    # Distant metastasis regression probability (%)
    distant_control_prob = round(min(88.0, max(5.0, abscopal_index * 0.85)), 1)

    return {
        "abscopal_systemic_index": abscopal_index,
        "distant_metastasis_control_prob_percent": distant_control_prob,
        "checkpoint_synergy": synergy_label
    }


def evaluate_neutralizing_antibody_impact(
    virus_family: str,
    pre_existing_nab_percent: float,
    delivery_route: str,
    carrier_cell_shielding: bool = False
) -> Dict[str, Any]:
    """
    Models the pharmacokinetic impact of host Neutralizing Antibodies (NAbs):
    - High anti-HSV-1 or anti-Ad5 IgG titers in adult human populations (~60-80%).
    - Intralesional delivery bypasses systemic vascular neutralizing clearance.
    - Systemic IV delivery suffers rapid complement and antibody neutralization unless shielded.
    """
    delivery_lower = delivery_route.lower()
    is_local = any(term in delivery_lower for term in ["intralesional", "intratumoral", "stereotactic", "eus", "intravesical"])

    if is_local:
        effective_neutralization_penalty = pre_existing_nab_percent * 0.15  # 85% protection due to local compartmentalization
        delivery_advantage = "Favorable: Compartmentalized/Intralesional injection bypasses circulating IgG"
    else:
        if carrier_cell_shielding:
            effective_neutralization_penalty = pre_existing_nab_percent * 0.30  # Carrier cells (MSC/monocytes) shield virus
            delivery_advantage = "Shielded: Carrier Cell / Nanoparticle trojan horse evades vascular NAbs"
        else:
            effective_neutralization_penalty = pre_existing_nab_percent * 0.75  # High neutralization in vascular space
            delivery_advantage = "Challenged: Systemic IV route exposed to high titer circulating NAbs and liver Kupffer uptake"

    effective_payload_fraction = max(0.10, min(1.0, 1.0 - (effective_neutralization_penalty / 100.0)))
    effective_payload_percent = round(effective_payload_fraction * 100.0, 1)

    return {
        "pre_existing_nab_seroprevalence": pre_existing_nab_percent,
        "effective_payload_percent": effective_payload_percent,
        "neutralization_attenuation_percent": round(100.0 - effective_payload_percent, 1),
        "delivery_assessment": delivery_advantage
    }


def analyze_oncolytic_regimen(
    regimen_id: str,
    custom_dose_pfu: Optional[float] = None,
    custom_tumor_burden: Optional[float] = None,
    anti_pd1_combination: bool = True,
    carrier_shielding: bool = False
) -> Dict[str, Any]:
    """
    Main orchestration function analyzing complete Oncolytic Virotherapy mechanics.
    """
    construct = CANONICAL_OV_CONSTRUCTS.get(regimen_id, CANONICAL_OV_CONSTRUCTS["TVEC_IMLYGIC_MELANOMA"])

    dose = custom_dose_pfu if custom_dose_pfu is not None else construct["default_dose_pfu"]
    tumor_cells = custom_tumor_burden if custom_tumor_burden is not None else 1.0e8

    # 1. NAb impact
    nab_result = evaluate_neutralizing_antibody_impact(
        virus_family=construct["virus_family"],
        pre_existing_nab_percent=construct["nab_seroprevalence_percent"],
        delivery_route=construct["delivery_route"],
        carrier_cell_shielding=carrier_shielding
    )

    # Effective starting dose taking NAb into account
    effective_dose = dose * (nab_result["effective_payload_percent"] / 100.0)

    # 2. Kinetic viral replication simulation
    burst = construct["burst_size"]
    lysis_rate = 24.0 / max(8.0, construct["lysis_cycle_hours"])

    kinetics_result = model_oncolytic_viral_kinetics(
        initial_tumor_cells=tumor_cells,
        initial_dose_pfu=effective_dose,
        infectivity_rate=1.0e-9,
        burst_size=burst,
        lysis_rate_per_day=lysis_rate,
        viral_clearance_per_day=0.45,
        days=14
    )

    # 3. ICD Biomarkers
    icd_result = calculate_icd_biomarker_score(
        lysis_percent=kinetics_result["tumor_lysis_percent"],
        virus_family=construct["virus_family"],
        gmcsf_present=construct["gmcsf_expression"]
    )

    # 4. Abscopal Systemic Response
    abscopal_result = calculate_abscopal_systemic_index(
        local_lysis_score=kinetics_result["tumor_lysis_percent"],
        icd_score=icd_result["composite_icd_score"],
        delivery_route=construct["delivery_route"],
        anti_pd1_combination=anti_pd1_combination
    )

    # Composite Therapeutic Efficacy Score
    overall_efficacy = round(
        0.30 * kinetics_result["tumor_lysis_percent"] +
        0.25 * icd_result["composite_icd_score"] +
        0.25 * abscopal_result["abscopal_systemic_index"] +
        0.20 * construct["base_efficacy_score"],
        1
    )
    overall_efficacy = min(99.0, max(25.0, overall_efficacy))

    return {
        "construct": construct,
        "nab_pharmacokinetics": nab_result,
        "viral_kinetics": kinetics_result,
        "icd_biomarkers": icd_result,
        "abscopal_effect": abscopal_result,
        "overall_therapeutic_score": overall_efficacy
    }


class OncolyticVirotherapyEngine:
    """Class wrapper for module integration and testing."""

    def __init__(self):
        self.constructs = CANONICAL_OV_CONSTRUCTS

    def list_constructs(self) -> List[Dict[str, Any]]:
        return list(self.constructs.values())

    def get_construct(self, construct_id: str) -> Optional[Dict[str, Any]]:
        return self.constructs.get(construct_id)

    def evaluate_oncolytic_treatment(
        self,
        construct_id: str,
        custom_dose_pfu: Optional[float] = None,
        anti_pd1_combination: bool = True
    ) -> Dict[str, Any]:
        return analyze_oncolytic_regimen(
            regimen_id=construct_id,
            custom_dose_pfu=custom_dose_pfu,
            anti_pd1_combination=anti_pd1_combination
        )
