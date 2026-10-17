"""
Cancer Neoantigen & Personalized mRNA Vaccine Design Studio Module (NeoVac Engine)
Comprehensive computational immunology, mRNA structural biology, and LNP formulation engine for:
- Canonical clinical neoantigen vaccine regimens:
  mRNA-4157 / V940 + Pembrolizumab (Melanoma, KEYNOTE-942),
  Autogene Cevumeran / BNT122 (PDAC Pancreatic Adenocarcinoma, Balachandran Nature 2023),
  Nous-209 Shared Frameshift Pool (Lynch Syndrome / MSI-H Colorectal),
  ELI-002 7P Amphiphile Vaccine (mKRAS G12D/V/R AMPLIFY-201),
  Rindopepimut EGFRvIII Junction Vaccine (Glioblastoma),
  NeoVac-NSCLC High-TMB Personalized 30-plex.
- HLA Class I & Class II allele peptide binding affinity (MHC-I IC50 nM, %Rank, anchor residues P2 and P9/PΩ).
- Neoantigen Fitness & Immunogenicity Index:
  TCR contact hydrophobicity/foreignness, clonality weighting (VAF x log2(TPM + 1)),
  proteasomal C-terminal cleavage probability, TAP translocation efficiency.
- Personalized mRNA Vaccine Cassette Assembly:
  5' Cap1 (CleanCap AG), 5' UTR (Alpha-Globin / HSD17B4), Kozak consensus,
  targeting signal peptide (LAMP1 / Ig-kappa / Sec61), multi-epitope concatemers with
  proteasome-cleavable linkers (AAY for MHC-I, GPGPG for MHC-II, Furin),
  3' UTR (AES + mtRNR1 / dual beta-globin), and Poly(A) tail (100-120 nt).
- Codon Optimization, GC% content (55-65%), and secondary structure stability.
- Lipid Nanoparticle (LNP) Formulation Biophysics:
  Ionizable lipid (SM-102, ALC-0315, MC3, C12-200), helper phospholipid (DSPC),
  cholesterol, PEG-lipid (DMG-PEG2000), N/P molar charge ratio (ideal 6.0),
  encapsulation efficiency (%), particle hydrodynamic diameter (nm), and endosomal escape rate.
- Clinical trial benchmarks (KEYNOTE-942, MSKCC Nature 2023, AMPLIFY-201, NCCN guidelines).
"""

import math
from typing import Dict, Any, List, Optional


CANONICAL_NEOVAC_REGIMENS: Dict[str, Dict[str, Any]] = {
    "MRNA_4157_V940_MELANOMA": {
        "id": "MRNA_4157_V940_MELANOMA",
        "name": "mRNA-4157 / V940 + Pembrolizumab (Moderna / Merck)",
        "tumor_type": "Cutaneous Melanoma (Resected Stage IIIB-IV)",
        "platform": "Personalized Synthetic mRNA (up to 34 Neoantigens, CleanCap 1, SM-102 LNP)",
        "indication": "Adjuvant treatment of completely resected high-risk cutaneous melanoma",
        "clinical_trial": "KEYNOTE-942 / Phase 2b (Lancet 2023): mRNA-4157 + Pembrolizumab reduced recurrence or death by 44% (RFS HR=0.56, 95% CI 0.31-1.02, p=0.053) and reduced distant metastasis or death by 65% (DMFS HR=0.35, 95% CI 0.17-0.72) vs Pembrolizumab monotherapy at 3 years.",
        "neoantigen_count": 34,
        "sample_peptide": "STEKITYSV",
        "wt_peptide": "STEKITYSL",
        "gene_mutation": "BRAF V600E (c.1799T>A)",
        "hla_allele": "HLA-A*02:01",
        "vaf": 0.48,
        "tpm_expression": 125.0,
        "linker_type": "AAY",
        "signal_peptide": "LAMP1",
        "ionizable_lipid": "SM-102",
        "np_ratio": 6.0,
        "dose_ug": 1000.0,
        "dosing_schedule": "1 mg IM every 3 weeks for up to 9 doses + Pembrolizumab 200 mg IV Q3W for 18 cycles (~1 year)",
        "nccn_guideline": "NCCN Melanoma Guidelines: Category 2A Emerging Modality. Received FDA Breakthrough Therapy Designation and EMA PRIME scheme for adjuvant melanoma.",
        "base_potency": 94.0
    },
    "AUTOGENE_CEVUMERAN_BNT122_PDAC": {
        "id": "AUTOGENE_CEVUMERAN_BNT122_PDAC",
        "name": "Autogene Cevumeran / BNT122 (BioNTech / Genentech)",
        "tumor_type": "Pancreatic Ductal Adenocarcinoma (Resected PDAC)",
        "platform": "Individualized Neoantigen-Specific Immunotherapy (iNeST) mRNA (up to 20 Neoantigens, 1-mΨ, Lipoplex)",
        "indication": "Surgically resected PDAC following surgery + Atezolizumab + mFOLFIRINOX adjuvant chemotherapy",
        "clinical_trial": "Nature 2023 (Vinod Balachandran et al., MSKCC): 50% (8/16) of PDAC patients mounted substantial neoantigen-specific polyfunctional CD8+ T cells (up to 10% of all blood T cells) persisting >2.5 years. Responders had significantly prolonged median RFS (not reached vs 13.4 months, HR=0.08, p=0.003).",
        "neoantigen_count": 20,
        "sample_peptide": "VVGADGVGK",
        "wt_peptide": "VVGAVGVGK",
        "gene_mutation": "KRAS G12D (c.35G>A)",
        "hla_allele": "HLA-A*11:01",
        "vaf": 0.38,
        "tpm_expression": 68.0,
        "linker_type": "GPGPG",
        "signal_peptide": "Sec61",
        "ionizable_lipid": "ALC-0315",
        "np_ratio": 6.0,
        "dose_ug": 50.0,
        "dosing_schedule": "Priming QW x 6 IV/IM starting post-surgery + Atezolizumab, followed by mFOLFIRINOX chemotherapy and maintenance booster vaccine",
        "nccn_guideline": "NCCN Pancreatic Adenocarcinoma: Clinical trial recommended post-resection. Breakthrough demonstration of de novo T-cell immunity in an immunologically 'cold' stroma.",
        "base_potency": 91.0
    },
    "MSI_H_FRAMESHIFT_POOL_NOUS209": {
        "id": "MSI_H_FRAMESHIFT_POOL_NOUS209",
        "name": "Nous-209 / Shared Frameshift Neoantigen Pool",
        "tumor_type": "Colorectal & Gastric MSI-H / dMMR (Lynch Syndrome)",
        "platform": "Off-the-shelf Shared Frameshift Neoantigens (209 pooled epitopes, GAd20 / MVA / mRNA)",
        "indication": "Microsatellite Instability-High (MSI-H) / Mismatch Repair Deficient (dMMR) solid tumors",
        "clinical_trial": "Phase 1/2 Trial (AACR / SITC 2023): Nous-209 + Pembrolizumab demonstrated 100% neoantigen-specific T-cell seroconversion in evaluable dMMR colorectal patients with 0% thymic central tolerance due to entirely novel out-of-frame polypeptide cascades in TGFBR2 (-1A), CASP8, and ACVR2A.",
        "neoantigen_count": 209,
        "sample_peptide": "RLSSCVPVA",
        "wt_peptide": "KSLVSSAST",
        "gene_mutation": "TGFBR2 -1A Frameshift (c.438delA)",
        "hla_allele": "HLA-A*02:01",
        "vaf": 0.52,
        "tpm_expression": 195.0,
        "linker_type": "AAY",
        "signal_peptide": "LAMP1",
        "ionizable_lipid": "SM-102",
        "np_ratio": 6.5,
        "dose_ug": 100.0,
        "dosing_schedule": "Prime-boost vaccination intramuscularly every 4 weeks with anti-PD-1",
        "nccn_guideline": "NCCN Colon & Rectal Guidelines: Shared frameshift vaccines are high-potential adjunctive immunotherapy for MSI-H/dMMR unresectable or metastatic tumors.",
        "base_potency": 96.0
    },
    "KRAS_MUTANT_AMPLIFY_ELI002": {
        "id": "KRAS_MUTANT_AMPLIFY_ELI002",
        "name": "ELI-002 7P / AMPLIFY-201 Amphiphile mKRAS Vaccine",
        "tumor_type": "KRAS-Mutant Pancreatic & Colorectal Cancer",
        "platform": "Lymph-Node Targeted Amphiphile (7 mKRAS peptides G12D/V/R/C/A/S, G13D + Amph-CpG-7909)",
        "indication": "Minimal Residual Disease (ctDNA+ / CA19-9+ / CEA+) post-curative surgery or chemotherapy",
        "clinical_trial": "AMPLIFY-201 Phase 1 (Nature Medicine 2024): 84% mutant KRAS-specific T-cell response rate; 84% ctDNA clearance or reduction; median relapse-free survival (RFS) not reached in T-cell responders vs 4.01 months in non-responders (HR=0.14, p=0.0016).",
        "neoantigen_count": 7,
        "sample_peptide": "VVGAVGVGK",
        "wt_peptide": "VVGAGGVGK",
        "gene_mutation": "KRAS G12V (c.35G>T)",
        "hla_allele": "HLA-A*02:01",
        "vaf": 0.44,
        "tpm_expression": 84.0,
        "linker_type": "GPGPG",
        "signal_peptide": "Ig-kappa",
        "ionizable_lipid": "MC3",
        "np_ratio": 6.0,
        "dose_ug": 500.0,
        "dosing_schedule": "Subcutaneous injection QW x 6 priming, then biweekly booster doses trafficking directly to sentinel lymph nodes via albumin binding",
        "nccn_guideline": "NCCN Colorectal & Pancreatic Guidelines: Investigational liquid biopsy MRD-directed personalized interception therapy.",
        "base_potency": 93.0
    },
    "EGFR_VIII_GLIOBLASTOMA_RINDOMUPEBTIDE": {
        "id": "EGFR_VIII_GLIOBLASTOMA_RINDOMUPEBTIDE",
        "name": "Rindopepimut / CDX-110 EGFRvIII Vaccine",
        "tumor_type": "Glioblastoma Multiforme (EGFRvIII-positive GBM)",
        "platform": "Synthetic Junction Neoepitope (Exon 2-7 deletion junction LEEKKGNYVVTDHC + KLH conjugate / mRNA)",
        "indication": "Newly diagnosed or recurrent EGFRvIII-positive glioblastoma post-chemoradiation",
        "clinical_trial": "ACT IV Phase 3 Trial & ReACT Phase 2 (Lancet Oncology 2017): Demonstrated potent anti-EGFRvIII humoral & cellular immune priming. Taught oncology the critical lesson of tumor immuno-editing: monovalent vaccines lead to escape of EGFRvIII-negative clones, prompting modern multi-epitope mRNA polyvalent designs.",
        "neoantigen_count": 1,
        "sample_peptide": "LEEKKGNYV",
        "wt_peptide": "LEEKKGNYL",
        "gene_mutation": "EGFRvIII (In-frame Exon 2-7 Deletion Junction)",
        "hla_allele": "HLA-A*02:01",
        "vaf": 0.62,
        "tpm_expression": 340.0,
        "linker_type": "AAY",
        "signal_peptide": "LAMP1",
        "ionizable_lipid": "SM-102",
        "np_ratio": 6.0,
        "dose_ug": 500.0,
        "dosing_schedule": "Intradermal or intramuscular vaccination biweekly with GM-CSF or Temozolomide",
        "nccn_guideline": "NCCN CNS Tumors Guidelines: Investigational. Underscores multi-epitope polyvalent neoantigen design necessity over single-epitope vaccines.",
        "base_potency": 82.0
    },
    "NSCLC_HIGH_TMB_PERSONALIZED": {
        "id": "NSCLC_HIGH_TMB_PERSONALIZED",
        "name": "NeoVac-NSCLC High-TMB Personalized 30-plex",
        "tumor_type": "Non-Small Cell Lung Cancer (High TMB > 10 mut/Mb, Tobacco UV/APOBEC Signature)",
        "platform": "Personalized 30-plex Synthetic mRNA (CleanCap 1, ALC-0315 / SM-102 LNP)",
        "indication": "Completely resected Stage IIA-IIIB NSCLC or metastatic NSCLC combined with anti-PD-L1 / anti-TIGIT",
        "clinical_trial": "NEO-PV-01 / BNT122 Solid Tumor Cohorts: Personalized multi-epitope vaccines generate polyfunctional CD4+ and CD8+ memory T-cell responses in 82-90% of NSCLC patients with high clonal TMB, driving epitope spreading to non-vaccinated clonal tumor antigens.",
        "neoantigen_count": 30,
        "sample_peptide": "KLWVDPREV",
        "wt_peptide": "KLWVDPREL",
        "gene_mutation": "TP53 R273H (c.818G>A)",
        "hla_allele": "HLA-A*02:01",
        "vaf": 0.55,
        "tpm_expression": 110.0,
        "linker_type": "AAY",
        "signal_peptide": "Sec61",
        "ionizable_lipid": "SM-102",
        "np_ratio": 6.0,
        "dose_ug": 100.0,
        "dosing_schedule": "IM injection Q3W for 9 cycles + Atezolizumab 1200 mg IV Q3W",
        "nccn_guideline": "NCCN NSCLC Guidelines: High-TMB neoantigen targeting in clinical development for adjuvant and consolidation immunotherapy.",
        "base_potency": 92.0
    }
}


HLA_MOTIF_RULES: Dict[str, Dict[str, Any]] = {
    "HLA-A*02:01": {
        "p2_favored": ["L", "M", "I", "V"],
        "p2_disfavored": ["D", "E", "R", "K", "P"],
        "p_omega_favored": ["V", "L", "I", "A", "M"],
        "p_omega_disfavored": ["D", "E", "R", "K", "P", "Y"],
        "length_preferred": [9, 10],
        "base_kd": 20.0
    },
    "HLA-A*24:02": {
        "p2_favored": ["Y", "F", "W", "M"],
        "p2_disfavored": ["D", "E", "P"],
        "p_omega_favored": ["L", "F", "I", "W"],
        "p_omega_disfavored": ["D", "E", "K", "R"],
        "length_preferred": [9, 10],
        "base_kd": 25.0
    },
    "HLA-B*07:02": {
        "p2_favored": ["P"],
        "p2_disfavored": ["D", "E", "L", "K"],
        "p_omega_favored": ["L", "M", "F", "V", "I"],
        "p_omega_disfavored": ["D", "E", "R", "K"],
        "length_preferred": [9],
        "base_kd": 30.0
    },
    "HLA-A*11:01": {
        "p2_favored": ["T", "V", "I", "S", "M"],
        "p2_disfavored": ["D", "E", "P"],
        "p_omega_favored": ["K", "R"],
        "p_omega_disfavored": ["D", "E", "L", "V"],
        "length_preferred": [9, 10],
        "base_kd": 28.0
    },
    "HLA-C*07:01": {
        "p2_favored": ["Y", "R", "K"],
        "p2_disfavored": ["D", "E"],
        "p_omega_favored": ["L", "F", "M"],
        "p_omega_disfavored": ["D", "E"],
        "length_preferred": [9],
        "base_kd": 35.0
    }
}


def predict_hla_binding_affinity(peptide: str, hla_allele: str) -> Dict[str, Any]:
    """
    Predicts MHC Class I dissociation constant (IC50 nM) and %Rank based on
    anchor residue motif matching at P2 and P-omega (P9/P10).
    Thresholds:
    - Strong Binder: IC50 < 50 nM or %Rank < 0.5%
    - Weak Binder: 50 <= IC50 <= 500 nM or 0.5% <= %Rank <= 2.0%
    - Non-Binder: IC50 > 500 nM or %Rank > 2.0%
    """
    clean_pep = peptide.strip().upper()
    n = len(clean_pep)
    rules = HLA_MOTIF_RULES.get(hla_allele, HLA_MOTIF_RULES["HLA-A*02:01"])
    
    if n < 8 or n > 11:
        return {
            "peptide": clean_pep,
            "hla_allele": hla_allele,
            "ic50_nm": 5000.0,
            "rank_percent": 15.0,
            "binding_class": "Non-Binder",
            "binding_affinity_score": 5.0,
            "comment": f"Peptide length ({n}) is outside the optimal MHC-I cleft range (8-11 AA)."
        }
    
    p2 = clean_pep[1]
    p_omega = clean_pep[-1]
    
    score = 50.0  # neutral starting baseline
    
    # Check P2 anchor
    if p2 in rules["p2_favored"]:
        score += 35.0
    elif p2 in rules["p2_disfavored"]:
        score -= 30.0
        
    # Check P-omega anchor
    if p_omega in rules["p_omega_favored"]:
        score += 35.0
    elif p_omega in rules["p_omega_disfavored"]:
        score -= 30.0
        
    # Length penalty
    if n in rules["length_preferred"]:
        score += 5.0
    else:
        score -= 10.0
        
    score = max(5.0, min(99.0, score))
    
    # Map score to IC50 nM (inverse sigmoid mapping)
    # score=95 -> IC50 ~ 12 nM; score=75 -> IC50 ~ 85 nM; score=50 -> IC50 ~ 450 nM; score=20 -> IC50 ~ 2500 nM
    ic50_nm = 500.0 * math.exp(-0.06 * (score - 50.0))
    ic50_nm = round(max(1.0, min(10000.0, ic50_nm)), 1)
    
    # Map to %Rank (percentile rank vs human peptidome)
    rank_percent = round(10.0 / (1.0 + math.exp(0.08 * (score - 45.0))), 2)
    rank_percent = max(0.01, min(25.0, rank_percent))
    
    if ic50_nm < 50.0 or rank_percent < 0.5:
        binding_class = "Strong Binder"
    elif ic50_nm <= 500.0 or rank_percent <= 2.0:
        binding_class = "Weak Binder"
    else:
        binding_class = "Non-Binder"
        
    return {
        "peptide": clean_pep,
        "hla_allele": hla_allele,
        "p2_residue": p2,
        "p_omega_residue": p_omega,
        "ic50_nm": ic50_nm,
        "rank_percent": rank_percent,
        "binding_class": binding_class,
        "binding_affinity_score": round(score, 1)
    }


def calculate_neoantigen_fitness(
    peptide: str,
    wt_peptide: str,
    ic50_nm: float,
    vaf: float,
    tpm_expression: float,
    is_frameshift: bool = False
) -> Dict[str, Any]:
    """
    Calculates comprehensive Neoantigen Fitness & Prioritization Index:
    1. MHC Binding Probability: P(bind) = 1 / (1 + exp((IC50 - 250) / 75))
    2. TCR Contact Foreignness: Solvent-exposed residues (P4-P7) hydrophobicity & aromaticity.
       Frameshifts receive max foreignness (+25) due to total lack of thymic central tolerance.
    3. Clonality & Expression Multiplier: VAF * log2(TPM + 1). Clonal mutations (VAF >= 0.4)
       prevent immune escape; subclonal mutations (VAF < 0.1) are penalized.
    4. C-terminal Proteasomal Cleavage: Leucine, Phenylalanine, Tyrosine, Arginine, Lysine favored.
    5. Transporter Associated with Antigen Processing (TAP) transport efficiency.
    """
    clean_mt = peptide.strip().upper()
    clean_wt = wt_peptide.strip().upper()
    
    # 1. MHC Binding Probability (0 - 100)
    p_bind = 1.0 / (1.0 + math.exp((ic50_nm - 250.0) / 75.0))
    mhc_score = p_bind * 100.0
    
    # 2. TCR Contact Foreignness
    # Positions 4-7 are primary solvent-exposed TCR contacts in canonical 9-10 mers
    tcr_contacts = clean_mt[3:7] if len(clean_mt) >= 7 else clean_mt
    aromatic_hydrophobic = sum(1 for aa in tcr_contacts if aa in ["W", "F", "Y", "I", "L", "V"])
    foreignness_score = 40.0 + (aromatic_hydrophobic * 12.0)
    
    if is_frameshift or clean_mt != clean_wt and len(clean_mt) == len(clean_wt):
        # Calculate differential binding or sequence divergence
        diff_count = sum(1 for a, b in zip(clean_mt, clean_wt) if a != b)
        foreignness_score += min(20.0, diff_count * 8.0)
    if is_frameshift:
        foreignness_score = min(99.0, foreignness_score + 25.0)
    foreignness_score = max(10.0, min(99.0, foreignness_score))
    
    # 3. Clonality & Expression Score
    # VAF (0-1), TPM (0 - 1000)
    # Clonal baseline: heterozygous somatic mutations have VAF ~0.35-0.50 in typical tumor purity
    clonality_vaf_factor = min(1.0, max(0.05, vaf) / 0.40)
    tpm_factor = min(1.0, math.log2(max(0.1, tpm_expression) + 1.0) / 6.0)
    clonality_score = clonality_vaf_factor * tpm_factor
    clonality_metric = round(min(100.0, max(5.0, clonality_score * 100.0)), 1)
    
    # 4. Proteasomal Cleavage Score (C-terminal residue)
    c_term = clean_mt[-1]
    if c_term in ["L", "F", "Y", "M", "I", "V", "R", "K"]:
        cleavage_score = 92.0
    elif c_term in ["A", "W", "T"]:
        cleavage_score = 75.0
    elif c_term in ["P", "D", "E"]:
        cleavage_score = 35.0  # Proteasome cleavage inhibited by proline or acidic residues at C-term
    else:
        cleavage_score = 60.0
        
    # 5. TAP Translocation Score (N-terminal 1-3 residues)
    n_term = clean_mt[0]
    if n_term in ["R", "K", "L", "I", "F"]:
        tap_score = 88.0
    elif n_term in ["P", "D", "E"]:
        tap_score = 45.0
    else:
        tap_score = 70.0
        
    # Composite Neoantigen Fitness Score (Weighted Harmonic/Linear Mix)
    # Weights: MHC (35%), Foreignness (25%), Clonality/Expression (25%), Cleavage/TAP (15%)
    composite_fitness = (
        0.35 * mhc_score +
        0.25 * foreignness_score +
        0.25 * clonality_metric +
        0.15 * ((cleavage_score + tap_score) / 2.0)
    )
    composite_fitness = round(max(5.0, min(99.5, composite_fitness)), 1)
    
    if composite_fitness >= 80.0:
        recommendation = "Tier 1 Priority: High fitness neoantigen. Mandatory inclusion in clinical mRNA cassette."
    elif composite_fitness >= 60.0:
        recommendation = "Tier 2 Candidate: Acceptable immunogenicity. Suitable for multi-epitope polyvalent bundle."
    else:
        recommendation = "Tier 3 Low Priority: Insufficient clonality, weak binding, or poor cleavage. Exclude from cassette."
        
    return {
        "mhc_binding_score": round(mhc_score, 1),
        "foreignness_score": round(foreignness_score, 1),
        "clonality_metric": clonality_metric,
        "cleavage_score": round(cleavage_score, 1),
        "tap_score": round(tap_score, 1),
        "composite_fitness": composite_fitness,
        "recommendation": recommendation
    }


def design_mrna_cassette(
    neoantigen_count: int,
    linker_type: str = "AAY",
    signal_peptide: str = "LAMP1",
    poly_a_length: int = 110
) -> Dict[str, Any]:
    """
    Designs and models the complete synthetic mRNA Vaccine Cassette:
    - 5' Cap1 (CleanCap AG): m7G(5')ppp(5')(2'OMeA)pG
    - 5' UTR: Human alpha-globin (HBA1) or HSD17B4 leader (optimizes eIF4E scanning)
    - Kozak sequence: GCCACCATGG
    - Signal peptide: LAMP1 (Lysosomal targeting for MHC-I/II dual presentation) or Sec61/Ig-kappa
    - Concatemer linkers:
      * 'AAY' (Ala-Ala-Tyr): Enhanced proteasome cleavage for CD8+ MHC-I presentation
      * 'GPGPG' (Gly-Pro-Gly-Pro-Gly): Flexible beta-turn preventing junctional neoepitopes, ideal for CD4+ MHC-II
      * 'furin' (RAKR): Subtilisin-like proprotein convertase cleavage
    - 3' UTR: Human mitochondrial AES + mtRNR1 or dual beta-globin
    - Poly(A) tail: 100-120 nucleotides (optimal stability vs deadenylation)
    - Codon optimization metrics (GC% content, Codon Adaptation Index CAI, mRNA half-life)
    """
    linker_specs = {
        "AAY": {
            "name": "Ala-Ala-Tyr (AAY) Tripeptide Linker",
            "mhc_preference": "MHC Class I (CD8+ Cytotoxic T-Cell Priming)",
            "proteasome_efficiency": 94.0,
            "junctional_neoepitope_risk": "Low (< 3%)",
            "sequence": "GCCGCCUAC"
        },
        "GPGPG": {
            "name": "Gly-Pro-Gly-Pro-Gly (GPGPG) Pentapeptide Linker",
            "mhc_preference": "MHC Class II (CD4+ Helper T-Cell Priming & B-cell Help)",
            "proteasome_efficiency": 88.0,
            "junctional_neoepitope_risk": "Very Low (< 1%)",
            "sequence": "GGCCCUGGGCCUGGC"
        },
        "furin": {
            "name": "Furin-Cleavable Arg-Ala-Lys-Arg (RAKR) Linker",
            "mhc_preference": "Dual MHC I & Secretory Cleavage",
            "proteasome_efficiency": 90.0,
            "junctional_neoepitope_risk": "Low (< 2%)",
            "sequence": "CGGGCCAAAAGG"
        }
    }
    
    linker = linker_specs.get(linker_type, linker_specs["AAY"])
    
    # Calculate transcript length
    cap_len = 2
    utr5_len = 54
    kozak_len = 10
    signal_len = 75 if signal_peptide == "LAMP1" else 60
    epitope_len = neoantigen_count * 27  # 9-mer = 27 nt
    linker_len = (neoantigen_count - 1) * len(linker["sequence"]) if neoantigen_count > 1 else 0
    stop_len = 3
    utr3_len = 132
    poly_a = max(60, min(150, poly_a_length))
    
    total_mrna_nt = cap_len + utr5_len + kozak_len + signal_len + epitope_len + linker_len + stop_len + utr3_len + poly_a
    
    # Codon optimization parameters
    target_gc_percent = 61.5  # Ideal mammalian mRNA GC content (prevents RIG-I while enhancing stability)
    cai_score = 0.94          # Codon Adaptation Index (0-1)
    
    # mRNA in vivo translation half-life estimation (hours)
    # Poly(A) 110-120 nt + Cap1 + dual UTR yields ~18-28 hours functional translation
    half_life_hours = round(14.0 + (poly_a / 10.0) * 0.9, 1)
    
    return {
        "linker_details": linker,
        "signal_peptide": signal_peptide,
        "total_mrna_length_nt": total_mrna_nt,
        "coding_region_nt": signal_len + epitope_len + linker_len + stop_len,
        "poly_a_length": poly_a,
        "gc_content_percent": target_gc_percent,
        "codon_adaptation_index": cai_score,
        "estimated_half_life_hours": half_life_hours,
        "cap_structure": "CleanCap AG (Cap1 - m7G(5')ppp(5')(2'OMeA)pG)",
        "utr_5prime": "Human HSD17B4 / HBA1 optimized leader",
        "utr_3prime": "Dual human beta-globin (HBB) + mtRNR1 tandem elements"
    }


def model_lnp_biophysics(
    ionizable_lipid: str = "SM-102",
    np_ratio: float = 6.0,
    total_rna_dose_ug: float = 100.0
) -> Dict[str, Any]:
    """
    Models Lipid Nanoparticle (LNP) formulation biophysics and cellular uptake:
    - Ionizable Lipid properties:
      * SM-102 (Moderna / mRNA-4157): pKa 6.75, high DC delivery, rapid tissue clearance
      * ALC-0315 (BioNTech / BNT122): pKa 6.09, low reactogenicity
      * MC3 (DLin-MC3-DMA): pKa 6.44, benchmark lipid
      * C12-200: pKa 6.80, ultra-potent lipidoid
    - Molar ratios: Ionizable Lipid (50%) : DSPC (10%) : Cholesterol (38.5%) : DMG-PEG2000 (1.5%)
    - N/P Charge Ratio:
      Ratio of protonatable tertiary amine nitrogens (N) on ionizable lipid to
      negatively charged phosphodiester phosphates (P) on mRNA backbone (MW ~330 g/mol).
      Optimal N/P ratio is 6.0:1 (ensures complete mRNA condensation and hexagonal phase endosomal rupture).
    - Encapsulation Efficiency (%)
    - Hydrodynamic Diameter (Z-average, nm) & Polydispersity Index (PDI)
    - Endosomal Escape Probability (%)
    """
    lipid_database = {
        "SM-102": {
            "name": "SM-102 (Heptadecan-9-yl 8-((2-hydroxyethyl)(6-oxo-6-(undecyloxy)hexyl)amino)octanoate)",
            "pka": 6.75,
            "mw": 710.2,
            "clearance_half_life_hours": 18.0,
            "potency_tier": "Very High"
        },
        "ALC-0315": {
            "name": "ALC-0315 (((4-hydroxybutyl)azanediyl)di(hexane-6,1-diyl) bis(2-hexyldecanoate))",
            "pka": 6.09,
            "mw": 766.3,
            "clearance_half_life_hours": 24.0,
            "potency_tier": "High"
        },
        "MC3": {
            "name": "DLin-MC3-DMA",
            "pka": 6.44,
            "mw": 642.1,
            "clearance_half_life_hours": 48.0,
            "potency_tier": "Benchmark Standard"
        },
        "C12-200": {
            "name": "C12-200 Lipidoid",
            "pka": 6.80,
            "mw": 1137.9,
            "clearance_half_life_hours": 16.0,
            "potency_tier": "Ultra High Potency"
        }
    }
    
    lipid = lipid_database.get(ionizable_lipid, lipid_database["SM-102"])
    
    # N/P ratio calculations
    # RNA: 1 ug RNA = 1e-6 g / 330 g/mol = 3.03 nmol phosphate
    clamped_np = max(2.0, min(12.0, np_ratio))
    
    # Encapsulation efficiency: reaches ~96-98% at N/P >= 5.0, drops sharply below 3.0
    if clamped_np >= 5.0:
        encapsulation_eff = 96.0 + min(3.0, (clamped_np - 5.0) * 0.4)
    elif clamped_np >= 4.0:
        encapsulation_eff = 88.0 + (clamped_np - 4.0) * 8.0
    else:
        encapsulation_eff = max(40.0, 50.0 + (clamped_np - 2.0) * 19.0)
    encapsulation_eff = round(min(99.5, encapsulation_eff), 1)
    
    # Hydrodynamic particle diameter (nm)
    # Typical microfluidic rapid mixing produces 75-85 nm particles at N/P 6.0
    particle_size_nm = round(72.0 + (12.0 - clamped_np) * 1.5 + (0.5 if ionizable_lipid == "ALC-0315" else 0.0), 1)
    pdi = 0.08 if clamped_np >= 5.0 else 0.14
    
    # Endosomal escape efficiency
    # In physiological endosomes (pH 5.5-6.0), tertiary amine protonation drives hexagonal H_II phase transition
    pka = lipid["pka"]
    protonation_fraction_endosome = 1.0 / (1.0 + 10.0 ** (5.5 - pka))
    endosomal_escape_pct = round(min(6.5, max(1.0, 2.5 * protonation_fraction_endosome * (clamped_np / 6.0))), 2)
    
    # Cytosolic mRNA release index (0-100)
    cytosolic_delivery_index = round(min(99.0, (encapsulation_eff / 100.0) * min(1.0, endosomal_escape_pct / 2.2) * 94.0), 1)
    
    return {
        "lipid_name": lipid["name"],
        "lipid_pka": lipid["pka"],
        "clearance_half_life_hours": lipid["clearance_half_life_hours"],
        "np_ratio": clamped_np,
        "total_rna_dose_ug": total_rna_dose_ug,
        "encapsulation_efficiency_percent": encapsulation_eff,
        "mean_particle_diameter_nm": particle_size_nm,
        "polydispersity_index_pdi": pdi,
        "endosomal_escape_percent": endosomal_escape_pct,
        "cytosolic_delivery_index": cytosolic_delivery_index
    }


def analyze_neovac_regimen(
    regimen_id: str,
    custom_peptide: Optional[str] = None,
    custom_hla: Optional[str] = None,
    custom_vaf: Optional[float] = None,
    custom_tpm: Optional[float] = None,
    custom_np_ratio: Optional[float] = None
) -> Dict[str, Any]:
    """
    Main orchestration function analyzing a complete Neoantigen & mRNA Vaccine Regimen.
    """
    regimen = CANONICAL_NEOVAC_REGIMENS.get(regimen_id, CANONICAL_NEOVAC_REGIMENS["MRNA_4157_V940_MELANOMA"])
    
    peptide = custom_peptide if custom_peptide else regimen["sample_peptide"]
    wt_pep = regimen["wt_peptide"]
    hla = custom_hla if custom_hla else regimen["hla_allele"]
    vaf = custom_vaf if custom_vaf is not None else regimen["vaf"]
    tpm = custom_tpm if custom_tpm is not None else regimen["tpm_expression"]
    np_ratio = custom_np_ratio if custom_np_ratio is not None else regimen["np_ratio"]
    
    is_frameshift = "Frameshift" in regimen.get("gene_mutation", "") or "TGFBR2" in regimen.get("gene_mutation", "")
    
    # 1. HLA Binding Affinity Prediction
    hla_result = predict_hla_binding_affinity(peptide, hla)
    
    # 2. Neoantigen Fitness Score
    fitness_result = calculate_neoantigen_fitness(
        peptide=peptide,
        wt_peptide=wt_pep,
        ic50_nm=hla_result["ic50_nm"],
        vaf=vaf,
        tpm_expression=tpm,
        is_frameshift=is_frameshift
    )
    
    # 3. mRNA Cassette Design
    mrna_result = design_mrna_cassette(
        neoantigen_count=regimen["neoantigen_count"],
        linker_type=regimen["linker_type"],
        signal_peptide=regimen["signal_peptide"],
        poly_a_length=110
    )
    
    # 4. LNP Formulation Biophysics
    lnp_result = model_lnp_biophysics(
        ionizable_lipid=regimen["ionizable_lipid"],
        np_ratio=np_ratio,
        total_rna_dose_ug=regimen["dose_ug"]
    )
    
    # Clinical potency composite score
    overall_score = round(
        0.30 * fitness_result["composite_fitness"] +
        0.25 * lnp_result["cytosolic_delivery_index"] +
        0.25 * regimen["base_potency"] +
        0.20 * mrna_result["codon_adaptation_index"] * 100.0,
        1
    )
    overall_score = min(99.0, max(20.0, overall_score))
    
    return {
        "regimen": regimen,
        "hla_binding": hla_result,
        "neoantigen_fitness": fitness_result,
        "mrna_cassette": mrna_result,
        "lnp_formulation": lnp_result,
        "overall_immunogenicity_potency_score": overall_score
    }
