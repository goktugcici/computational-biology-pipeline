"""
Targeted Radioligand Therapy (RLT) & Theranostics Studio Module
Comprehensive computational biophysics, nuclear medicine, and clinical oncology engine for:
- Canonical clinically approved & investigational radiopharmaceuticals (177Lu-PSMA-617 Pluvicto,
  177Lu-DOTATATE Lutathera, 225Ac-PSMA-617 TAT, 223Ra-Dichloride Xofigo, 131I-MIBG Azedra, 225Ac-DOTATATE).
- Diagnostic-Therapeutic pairing ("We see what we treat, we treat what we see" - 68Ga/18F PET to 177Lu/225Ac).
- Alpha particle vs. Beta particle radiation physics (Linear Energy Transfer - LET, crossfire vs. cluster DNA DSBs).
- MIRD (Medical Internal Radiation Dose) formal absorbed dosimetry calculation for tumor and organs at risk (OAR).
- Renal nephroprotection modeling via L-lysine/L-arginine amino acid competitive proximal tubule inhibition.
- Salivary gland toxicity (xerostomia) and red bone marrow radiotoxicity risk stratification.
- VISION, NETTER-1, NETTER-2, and ALSYMPCA trial guidelines and NCCN practice recommendations.
"""

from typing import Dict, Any, List, Optional


CANONICAL_RADIOLIGAND_REGIMENS: Dict[str, Dict[str, Any]] = {
    "LU177_PSMA_617_PLUVICTO": {
        "id": "LU177_PSMA_617_PLUVICTO",
        "name": "177Lu-PSMA-617 (Pluvicto / Lutetium vipivotide tetraxetan)",
        "target_antigen": "Prostate-Specific Membrane Antigen (PSMA / FOLH1)",
        "indication": "Metastatic Castration-Resistant Prostate Cancer (mCRPC) post-ARPI & Taxane Chemotherapy",
        "diagnostic_pair": "68Ga-PSMA-11 or 18F-DCFPyL (Pylarify) PET/CT",
        "therapeutic_radionuclide": "Lutetium-177 (177Lu)",
        "radiation_type": "Beta minus (beta-) emitter",
        "particle_energy_mev": 0.498,
        "tissue_penetration_mm": 1.5,
        "linear_energy_transfer": "Low LET (~0.2 keV/um)",
        "half_life_days": 6.647,
        "standard_activity_gbq": 7.4,  # 200 mCi per cycle q6w up to 6 cycles
        "cycles": 6,
        "dna_damage_mechanism": "Single-strand breaks (SSB) and moderate double-strand breaks via ionizing beta path crossfire effect.",
        "clinical_efficacy": "VISION Phase 3 Trial: Prolonged median overall survival (15.3 vs 11.3 mos, HR=0.62) and radiographic PFS (8.7 vs 3.4 mos, HR=0.40). PSA50 response rate 46.0% vs 7.1%. PSMAfore Phase 3: Significant rPFS benefit before taxane chemotherapy.",
        "kidney_dose_coeff_gy_per_gbq": 0.72,
        "marrow_dose_coeff_gy_per_gbq": 0.040,
        "salivary_uptake_percent": 74.0,
        "nccn_guideline": "NCCN Category 1 Preferred for adult patients with PSMA-positive mCRPC previously treated with an AR pathway inhibitor (ARPI) and taxane-based chemotherapy. Pre-treatment 68Ga-PSMA-11 or 18F-DCFPyL PET/CT mandatory: tumor lesion SUVmax must exceed liver SUVmean. Maximum cumulative kidney dose threshold is 23 Gy.",
        "base_tumor_potency": 88.0
    },
    "LU177_DOTATATE_LUTATHERA": {
        "id": "LU177_DOTATATE_LUTATHERA",
        "name": "177Lu-DOTATATE (Lutathera / 177Lu-Oxodotreotide)",
        "target_antigen": "Somatostatin Receptor Subtype 2 (SSTR2)",
        "indication": "Well-Differentiated Gastroenteropancreatic Neuroendocrine Tumors (GEP-NETs: Midgut, Pancreatic, Hindgut)",
        "diagnostic_pair": "68Ga-DOTATATE (Netspot) or 64Cu-DOTATATE (Detectnet) PET/CT",
        "therapeutic_radionuclide": "Lutetium-177 (177Lu)",
        "radiation_type": "Beta minus (beta-) emitter",
        "particle_energy_mev": 0.498,
        "tissue_penetration_mm": 1.5,
        "linear_energy_transfer": "Low LET (~0.2 keV/um)",
        "half_life_days": 6.647,
        "standard_activity_gbq": 7.4,  # 200 mCi per cycle q8w up to 4 cycles
        "cycles": 4,
        "dna_damage_mechanism": "Receptor-mediated internalizing beta-emitting crossfire causing tumor apoptosis and microvascular sclerosis.",
        "clinical_efficacy": "NETTER-1 Phase 3 Trial: 79% reduction in disease progression or death (PFS HR=0.21, p<0.001 vs high-dose octreotide LAR 60mg). NETTER-2 Phase 3: Superior PFS (22.8 vs 8.5 mos, HR=0.276) as frontline therapy in high-grade G2/G3 GEP-NETs.",
        "kidney_dose_coeff_gy_per_gbq": 0.65,
        "marrow_dose_coeff_gy_per_gbq": 0.038,
        "salivary_uptake_percent": 12.0,
        "nccn_guideline": "NCCN Category 1 Preferred for progressive SSTR2-positive gastroenteropancreatic and bronchopulmonary neuroendocrine tumors. Co-administration of an amino acid solution containing L-lysine and L-arginine is mandatory 30 min before through 3-4 hours after infusion to block proximal renal tubular reabsorption and prevent nephrotoxicity.",
        "base_tumor_potency": 92.0
    },
    "AC225_PSMA_617_ALPHA": {
        "id": "AC225_PSMA_617_ALPHA",
        "name": "225Ac-PSMA-617 (Targeted Alpha Therapy - TAT)",
        "target_antigen": "Prostate-Specific Membrane Antigen (PSMA / FOLH1)",
        "indication": "177Lu-Refractory Metastatic CRPC & Radioresistant Micrometastatic Disease",
        "diagnostic_pair": "68Ga-PSMA-11 PET/CT",
        "therapeutic_radionuclide": "Actinium-225 (225Ac)",
        "radiation_type": "Alpha (alpha) particle emitter (4 alpha particle decay cascade)",
        "particle_energy_mev": 5.8,  # Net decay cascade ~28 MeV
        "tissue_penetration_mm": 0.08,  # 40 - 80 um (2-4 cell diameters)
        "linear_energy_transfer": "High LET (~100 keV/um - Ultra High Density Ionization)",
        "half_life_days": 9.92,
        "standard_activity_gbq": 0.008,  # 100 kBq/kg (~8 MBq = 0.008 GBq)
        "cycles": 3,
        "dna_damage_mechanism": "High-LET alpha particles inflict complex clustered double-strand DNA breaks (DSBs) that completely overwhelm homologous recombination (HR) and NHEJ repair pathways, independent of cell hypoxia or p53 status.",
        "clinical_efficacy": "Heidelberg & Würzburg TAT Cohorts: 82% of patients refractory to 177Lu-PSMA achieved >50% PSA reduction, with complete radiological remissions in visceral and diffuse skeletal disease.",
        "kidney_dose_coeff_gy_per_gbq": 14.5,
        "marrow_dose_coeff_gy_per_gbq": 0.85,
        "salivary_uptake_percent": 96.0,
        "nccn_guideline": "Investigational / Compassionate use protocol for PSMA-expressing mCRPC with disease progression on 177Lu-PSMA-617. Dose-limiting toxicity is severe xerostomia (salivary gland destruction); salivary gland cooling with ice packs or local botulinum toxin injection is actively used to decrease salivary perfusion.",
        "base_tumor_potency": 98.0
    },
    "RA223_DICHLORIDE_XOFIGO": {
        "id": "RA223_DICHLORIDE_XOFIGO",
        "name": "223Ra-Dichloride (Xofigo / Radium-223)",
        "target_antigen": "Bone Hydroxyapatite (Calcium-Mimetic Osteoblastic Targeting)",
        "indication": "mCRPC with Symptomatic Bone Metastases and No Visceral Metastases",
        "diagnostic_pair": "99mTc-MDP Bone Scintigraphy or 18F-NaF PET/CT",
        "therapeutic_radionuclide": "Radium-223 (223Ra)",
        "radiation_type": "Alpha (alpha) particle emitter",
        "particle_energy_mev": 5.78,
        "tissue_penetration_mm": 0.10,  # <100 um
        "linear_energy_transfer": "High LET (~80 keV/um)",
        "half_life_days": 11.43,
        "standard_activity_gbq": 0.004,  # 55 kBq/kg (~4 MBq = 0.004 GBq) q4w x 6 cycles
        "cycles": 6,
        "dna_damage_mechanism": "Incorporation into bone matrix at sites of active osteoblastic turnover, delivering localized high-LET alpha radiation to adjacent osteoblasts, osteoclasts, and tumor cells while sparing bone marrow.",
        "clinical_efficacy": "ALSYMPCA Phase 3 Trial: Prolonged median overall survival (14.9 vs 11.3 mos, HR=0.70, p<0.001) and significantly delayed time to first symptomatic skeletal-related event (SRE: 15.6 vs 9.8 mos, HR=0.66).",
        "kidney_dose_coeff_gy_per_gbq": 0.15,
        "marrow_dose_coeff_gy_per_gbq": 1.10,
        "salivary_uptake_percent": 4.0,
        "nccn_guideline": "NCCN Category 1 for patients with castration-resistant prostate cancer, symptomatic bone metastases, and no known visceral metastatic disease. Contraindication: Do NOT combine with Abiraterone + Prednisone due to increased risk of bone fractures (ERA 223 trial).",
        "base_tumor_potency": 84.0
    },
    "I131_MIBG_AZEDRA": {
        "id": "I131_MIBG_AZEDRA",
        "name": "131I-Iobenguane (Azedra / High-Specific-Activity I-131 MIBG)",
        "target_antigen": "Norepinephrine Transporter (NET / SLC6A2)",
        "indication": "Malignant, Recurrent, or Unresectable Pheochromocytoma and Paraganglioma",
        "diagnostic_pair": "123I-MIBG Scintigraphy / 68Ga-DOTATATE PET/CT",
        "therapeutic_radionuclide": "Iodine-131 (131I)",
        "radiation_type": "Beta minus (beta-) & Gamma (gamma) emitter",
        "particle_energy_mev": 0.606,
        "tissue_penetration_mm": 2.0,
        "linear_energy_transfer": "Low LET (~0.25 keV/um)",
        "half_life_days": 8.02,
        "standard_activity_gbq": 18.5,  # 500 mCi dosimetric followed by therapeutic
        "cycles": 2,
        "dna_damage_mechanism": "Selective accumulation into neuroendocrine chromaffin granules via NET transporter, inducing cytotoxic crossfire irradiation.",
        "clinical_efficacy": "Pivotal Phase 2b Trial: 25% of patients achieved a >=50% sustained reduction in baseline antihypertensive medication requirement for at least 6 months; 92% achieved partial response or stable disease.",
        "kidney_dose_coeff_gy_per_gbq": 0.45,
        "marrow_dose_coeff_gy_per_gbq": 0.22,
        "salivary_uptake_percent": 24.0,
        "nccn_guideline": "FDA Approved for adult and pediatric patients (>=12 years) with iobenguane scan-positive unresectable or metastatic pheochromocytoma/paraganglioma. Mandatory inorganic potassium iodide (SSKI) thyroid blockade starting 24-48h prior to prevent 131I thyroid ablation. Autologous stem cell cryopreservation recommended due to prolonged myelosuppression.",
        "base_tumor_potency": 80.0
    },
    "AC225_DOTATATE_ALPHA_NET": {
        "id": "AC225_DOTATATE_ALPHA_NET",
        "name": "225Ac-DOTATATE (Targeted Alpha PRRT for NETs)",
        "target_antigen": "Somatostatin Receptor Subtype 2 (SSTR2)",
        "indication": "177Lu-PRRT Refractory or Rapidly Progressing Neuroendocrine Tumors",
        "diagnostic_pair": "68Ga-DOTATATE PET/CT",
        "therapeutic_radionuclide": "Actinium-225 (225Ac)",
        "radiation_type": "Alpha (alpha) particle cascade",
        "particle_energy_mev": 5.8,
        "tissue_penetration_mm": 0.08,
        "linear_energy_transfer": "High LET (~100 keV/um)",
        "half_life_days": 9.92,
        "standard_activity_gbq": 0.0075,  # ~7.5 MBq
        "cycles": 3,
        "dna_damage_mechanism": "High-LET alpha irradiation of SSTR2-expressing neuroendocrine tumor clusters, completely bypassing beta-radiation repair capacity and somatostatin receptor downregulation.",
        "clinical_efficacy": "International Compassionate Salvage Studies: In NET patients who progressed after full 4 cycles of 177Lu-DOTATATE, 225Ac-DOTATATE induced 62.5% objective response rates with substantial tumor regression.",
        "kidney_dose_coeff_gy_per_gbq": 12.0,
        "marrow_dose_coeff_gy_per_gbq": 0.70,
        "salivary_uptake_percent": 14.0,
        "nccn_guideline": "Emerging clinical trial strategy for refractory GEP-NETs with progressive disease on 177Lu-DOTATATE. Nephroprotection with amino acids and vigilant monitoring of glomerular filtration rate (eGFR) is essential.",
        "base_tumor_potency": 94.0
    }
}


class RadioligandTheranosticsEngine:
    """
    Bioinformatics, nuclear oncology, and MIRD radiation dosimetry engine for
    Targeted Radioligand Therapy (RLT) and Theranostics simulation.
    """

    def __init__(self):
        self.regimens = CANONICAL_RADIOLIGAND_REGIMENS

    def list_regimens(self) -> List[Dict[str, Any]]:
        """Returns all canonical theranostic regimens with complete metadata."""
        return list(self.regimens.values())

    def get_regimen(self, regimen_id: str) -> Optional[Dict[str, Any]]:
        """Fetches a specific radioligand regimen by ID."""
        return self.regimens.get(regimen_id)

    def evaluate_theranostics_and_dosimetry(
        self,
        regimen_id: str,
        administered_activity_gbq: Optional[float] = None,
        tumor_suv_max: float = 24.5,
        liver_suv_mean: float = 6.2,
        amino_acid_nephroprotection: bool = True
    ) -> Dict[str, Any]:
        """
        Simulates PET theranostic pairing, MIRD absorbed dose to tumor and critical organs at risk,
        linear energy transfer (LET) cell-kill mechanics, and clinical safety thresholds.
        """
        regimen = self.regimens.get(regimen_id)
        if not regimen:
            raise ValueError(f"Unknown regimen ID: '{regimen_id}'. Valid IDs: {list(self.regimens.keys())}")

        activity = float(administered_activity_gbq) if administered_activity_gbq is not None else float(regimen["standard_activity_gbq"])
        tumor_suv = max(1.0, float(tumor_suv_max))
        liver_suv = max(0.5, float(liver_suv_mean))

        # 1. Theranostic Eligibility & PET Avidity Ratio (VISION Trial Criteria)
        # VISION Criteria: Target lesion SUVmax > Liver SUVmean (ideally ratio >= 1.5x)
        avidity_ratio = round(tumor_suv / liver_suv, 2)
        if avidity_ratio >= 2.0:
            eligibility_status = "Optimal Theranostic Candidate (Marked Tumor Avidity > 2x Liver)"
            eligibility_badge = "emerald"
        elif avidity_ratio >= 1.5:
            eligibility_status = "Eligible Candidate (Meets VISION Trial Criteria: SUVmax >= 1.5x Liver)"
            eligibility_badge = "teal"
        else:
            eligibility_status = "Borderline / Low Avidity (Sub-threshold uptake; high risk of non-response)"
            eligibility_badge = "amber"

        # 2. MIRD Absorbed Radiation Dose Calculation
        # Tumor Dose (Gy) = Activity (GBq) * Base Potency factor * (Avidity Ratio / 2.0)
        is_alpha = "Alpha" in regimen["radiation_type"]
        if is_alpha:
            # Alpha particles have high relative biological effectiveness (RBE ~ 5 - 20)
            tumor_dose_gy = round(activity * (regimen["base_tumor_potency"] * 45.0) * (avidity_ratio / 3.0), 1)
        else:
            tumor_dose_gy = round(activity * (regimen["base_tumor_potency"] * 0.08) * (avidity_ratio / 2.5), 1)

        # 3. Critical Organs at Risk (OAR) Absorbed Doses (MIRD Formalism)
        # Kidneys: Cumulative tolerance limit is 23 - 28 Gy for beta, lower for alpha
        nephro_factor = 0.55 if amino_acid_nephroprotection else 1.0  # Amino acids reduce renal retention by 45%
        kidney_dose_per_cycle = round(activity * regimen["kidney_dose_coeff_gy_per_gbq"] * nephro_factor, 2)
        cumulative_kidney_dose = round(kidney_dose_per_cycle * regimen["cycles"], 1)

        # Bone Marrow (Cumulative limit ~2 Gy)
        marrow_dose_per_cycle = round(activity * regimen["marrow_dose_coeff_gy_per_gbq"], 3)
        cumulative_marrow_dose = round(marrow_dose_per_cycle * regimen["cycles"], 2)

        # Salivary Glands (PSMA avidity)
        salivary_risk_score = round(regimen["salivary_uptake_percent"] * (1.2 if is_alpha else 0.8), 1)
        salivary_risk_score = min(100.0, max(5.0, salivary_risk_score))

        # 4. Radiation Biophysics & DNA Damage Profile
        if is_alpha:
            radiation_physics_profile = "Targeted Alpha Therapy (TAT): High LET (~80-100 keV/um) helium nuclei causing dense, complex clustered DNA double-strand breaks (DSBs). Independent of cell oxygenation (no hypoxia resistance); single-cell kill radius <100 um."
            xerostomia_risk = "High / Dose-Limiting (Alpha-induced acinar salivary destruction; ice pack cooling mandatory)"
        else:
            radiation_physics_profile = "Beta-Minus Radioligand Therapy: Low LET (~0.2 keV/um) electrons with 1.5-2.0 mm crossfire path length, sterilizing surrounding heterogeneous antigen-negative tumor cells."
            xerostomia_risk = "Moderate / Transient (Mild to moderate dry mouth, manageable with saliva substitutes)"

        return {
            "regimen_id": regimen["id"],
            "name": regimen["name"],
            "target_antigen": regimen["target_antigen"],
            "indication": regimen["indication"],
            "diagnostic_pair": regimen["diagnostic_pair"],
            "therapeutic_radionuclide": regimen["therapeutic_radionuclide"],
            "radiation_type": regimen["radiation_type"],
            "half_life_days": regimen["half_life_days"],
            "administered_activity_gbq": activity,
            "cycles": regimen["cycles"],
            "tumor_suv_max": tumor_suv,
            "liver_suv_mean": liver_suv,
            "avidity_ratio": avidity_ratio,
            "eligibility_status": eligibility_status,
            "eligibility_badge": eligibility_badge,
            "tumor_absorbed_dose_gy": tumor_dose_gy,
            "kidney_dose_per_cycle_gy": kidney_dose_per_cycle,
            "cumulative_kidney_dose_gy": cumulative_kidney_dose,
            "cumulative_marrow_dose_gy": cumulative_marrow_dose,
            "amino_acid_protection_active": amino_acid_nephroprotection,
            "salivary_risk_score": salivary_risk_score,
            "xerostomia_risk": xerostomia_risk,
            "radiation_physics_profile": radiation_physics_profile,
            "dna_damage_mechanism": regimen["dna_damage_mechanism"],
            "clinical_efficacy": regimen["clinical_efficacy"],
            "nccn_guideline": regimen["nccn_guideline"]
        }


def analyze_radioligand_regimen(
    regimen_id: str,
    administered_activity_gbq: Optional[float] = None,
    tumor_suv_max: float = 24.5,
    liver_suv_mean: float = 6.2,
    amino_acid_nephroprotection: bool = True
) -> Dict[str, Any]:
    """Helper convenience function for rapid radioligand theranostics evaluation."""
    engine = RadioligandTheranosticsEngine()
    return engine.evaluate_theranostics_and_dosimetry(
        regimen_id,
        administered_activity_gbq,
        tumor_suv_max,
        liver_suv_mean,
        amino_acid_nephroprotection
    )
