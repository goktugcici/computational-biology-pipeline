"""
Epigenetics, DNA Methylation & Epigenomic Biomarkers Studio Module
Comprehensive computational biology and clinical oncology engine for:
- MGMT promoter methylation & Temozolomide (TMZ) alkylator sensitivity in Glioblastoma.
- MLH1 promoter hypermethylation vs. germline Lynch syndrome triage in MSI-H Colorectal/Endometrial Cancer.
- Multi-Cancer Early Detection (MCED) cell-free DNA (cfDNA) methylation and Tissue-of-Origin (TOO) deconvolution.
- Histone methylation & chromatin remodeling vulnerabilities (SMARCB1/INI1 loss, EZH2 gain-of-function, Tazemetostat).
- IDH1/IDH2 mutation, 2-hydroxyglutarate (2-HG) accumulation, and Glioma CpG Island Methylator Phenotype (G-CIMP).
"""

from typing import Dict, Any, List


CANONICAL_EPIGENETIC_SCENARIOS: Dict[str, Dict[str, Any]] = {
    "GBM_MGMT_METHYLATED": {
        "id": "GBM_MGMT_METHYLATED",
        "title": "Glioblastoma: MGMT Promoter Hypermethylation (Alkylator Sensitization)",
        "cancer_type": "Glioblastoma Multiforme (WHO Grade 4)",
        "primary_gene": "MGMT (O-6-Methylguanine-DNA Methyltransferase)",
        "cpg_region": "Promoter CpG Island (Exon 1 / Intron 1 boundary, 16 CpGs)",
        "methylation_level_percent": 68.5,
        "assay_technique": "Bisulfite Pyrosequencing / ddPCR",
        "silencing_status": "Transcriptional Gene Silencing Confirmed (IHC 0% nuclear)",
        "clinical_significance": "Predictive of high sensitivity to Temozolomide (TMZ) alkylating chemotherapy and radiotherapy (Stupp regimen). Median overall survival prolonged to 21.7 months vs. 15.3 months in unmethylated cohorts.",
        "pathway_impact": "Loss of O-6-alkylguanine DNA repair capacity leaves Temozolomide-induced O6-methylguanine adducts unrepaired, generating lethal double-strand breaks during replication.",
        "nccn_esmo_guideline": "Category 1 evidence: Standard Stupp Protocol (RT 60 Gy + concurrent daily TMZ followed by 6 cycles of adjuvant maintenance TMZ). Consider clinical trial enrollment for TMZ + PARP or TMZ + CCNU combinations.",
        "associated_markers": {
            "IDH1_Status": "Wild-Type (p.R132 WT)",
            "TERT_Promoter": "Mutated (C228T)",
            "EGFR_Amplification": "Present (EGFRvIII negative)",
            "PTEN_Status": "Homozygous Deletion"
        },
        "target_therapies": [
            {"drug": "Temozolomide (TMZ)", "class": "Oral Alkylating Agent", "expected_response": "High Sensitivity / Prolonged OS", "tier": "Standard of Care (Tier 1A)"},
            {"drug": "Lomustine (CCNU)", "class": "Nitrosourea Alkylator", "expected_response": "Secondary Sensitivity in Recurrence", "tier": "Tier 2A"},
            {"drug": "Tumor Treating Fields (TTFields)", "class": "Antimitotic Biophysical Therapy", "expected_response": "Additive Synergistic OS Extension", "tier": "Category 1"}
        ],
        "cpg_methylation_profile": [
            {"gene": "MGMT", "beta_value": 0.69, "status": "Hypermethylated"},
            {"gene": "MLH1", "beta_value": 0.04, "status": "Unmethylated"},
            {"gene": "SEPT9", "beta_value": 0.06, "status": "Unmethylated"},
            {"gene": "GSTP1", "beta_value": 0.08, "status": "Unmethylated"},
            {"gene": "CDKN2A", "beta_value": 0.85, "status": "Hypermethylated / Deleted"},
            {"gene": "RASSF1A", "beta_value": 0.52, "status": "Hypermethylated"},
            {"gene": "BRCA1", "beta_value": 0.03, "status": "Unmethylated"},
            {"gene": "SHOX2", "beta_value": 0.07, "status": "Unmethylated"}
        ]
    },
    "CRC_MLH1_SPORADIC": {
        "id": "CRC_MLH1_SPORADIC",
        "title": "Colorectal Carcinoma: MLH1 Promoter Hypermethylation & Lynch Triage",
        "cancer_type": "Right-Sided Colorectal Adenocarcinoma (MSI-H / dMMR)",
        "primary_gene": "MLH1 (MutL Homolog 1)",
        "cpg_region": "Promoter CpG Island (Deng Region C/D, c.-248 to c.-178)",
        "methylation_level_percent": 84.0,
        "assay_technique": "Bisulfite PCR-Pyrosequencing / MS-MLPA",
        "silencing_status": "Epigenetic Silencing Confirmed (Absent MLH1/PMS2 on IHC)",
        "clinical_significance": "Crucial clinical triage: Confirms sporadic MSI-H/dMMR etiology driven by CIMP phenotype and rules out hereditary Lynch syndrome (HNPCC), sparing patient and family unnecessary germline testing. Confirms high response to anti-PD-1 immunotherapy.",
        "pathway_impact": "Epigenetic hypermethylation of MLH1 promoter prevents MutLalpha complex formation, producing profound DNA mismatch repair deficiency, hypermutation, and high neoantigen burden.",
        "nccn_esmo_guideline": "NCCN Colon Cancer v2.2024: MLH1 hypermethylation + BRAF V600E confirms sporadic cancer. Recommend first-line immune checkpoint blockade (Pembrolizumab or Nivolumab + Ipilimumab). Germline MLH1 sequencing not indicated.",
        "associated_markers": {
            "BRAF_Status": "Mutated (p.V600E)",
            "MMR_IHC": "Loss of MLH1 and PMS2; Intact MSH2 and MSH6",
            "MSI_Status": "MSI-High (Pentaplex panel 5/5 unstable loci)",
            "CIMP_Status": "CIMP-High (5/5 Weisenberger panel methylated)"
        },
        "target_therapies": [
            {"drug": "Pembrolizumab", "class": "Anti-PD-1 Monoclonal Antibody", "expected_response": "High (KEYNOTE-177 Paradigm: ORR ~44%, PFS HR 0.60)", "tier": "FDA First-Line (Tier 1A)"},
            {"drug": "Nivolumab + Ipilimumab", "class": "Dual Anti-PD-1 + Anti-CTLA-4", "expected_response": "High Durable Complete Responses (CheckMate-8HW)", "tier": "FDA First-Line (Tier 1A)"},
            {"drug": "5-FU Adjuvant Monotherapy", "class": "Fluoropyrimidine", "expected_response": "Resistance / Lack of Benefit in Stage II dMMR", "tier": "Not Recommended"}
        ],
        "cpg_methylation_profile": [
            {"gene": "MGMT", "beta_value": 0.12, "status": "Unmethylated"},
            {"gene": "MLH1", "beta_value": 0.84, "status": "Hypermethylated"},
            {"gene": "SEPT9", "beta_value": 0.76, "status": "Hypermethylated"},
            {"gene": "GSTP1", "beta_value": 0.15, "status": "Unmethylated"},
            {"gene": "CDKN2A", "beta_value": 0.78, "status": "Hypermethylated"},
            {"gene": "RASSF1A", "beta_value": 0.65, "status": "Hypermethylated"},
            {"gene": "BRCA1", "beta_value": 0.05, "status": "Unmethylated"},
            {"gene": "SHOX2", "beta_value": 0.09, "status": "Unmethylated"}
        ]
    },
    "MCED_LIQUID_BIOPSY_CRC": {
        "id": "MCED_LIQUID_BIOPSY_CRC",
        "title": "MCED Liquid Biopsy: cfDNA Methylation & Tissue-of-Origin Deconvolution",
        "cancer_type": "Multi-Cancer Early Detection (Circulating cfDNA Plasma)",
        "primary_gene": "SEPT9, BCAT1, IKZF1 (Colorectal Epigenetic Signature)",
        "cpg_region": "Promoter & Enhancer Methylation Densely Clustered CpG Regions",
        "methylation_level_percent": 76.2,
        "assay_technique": "Targeted Enzymatic Methyl-Seq (EM-seq) / cfDNA Methylation Array",
        "silencing_status": "Positive Circulating Epigenomic Signal Detected",
        "clinical_significance": "Demonstrates positive cancer signal detection with machine-learning tissue-of-origin (TOO) localization to the Lower GI / Colorectal tract (92.4% TOO confidence). Enables non-invasive early-stage interceptive detection.",
        "pathway_impact": "Tumor cell apoptosis releases densely hypermethylated promoter DNA fragments into circulating plasma, creating distinct fragmentation length and CpG methylation density patterns distinguishable from leukocyte background.",
        "nccn_esmo_guideline": "MCED liquid biopsy positive result: Prioritize diagnostic colonoscopy with complete bowel preparation. If negative, evaluate upper GI endoscopy and pelvic CT.",
        "associated_markers": {
            "cfDNA_Concentration": "14.8 ng/mL plasma",
            "Fragment_Size_Peak": "167 bp (Nucleosomal monosome)",
            "Too_Confidence": "Colorectal Tract (92.4%), Gastric (4.1%), Pancreatic (2.0%)",
            "Circulating_EpiScore": "78.4 / 100 (Threshold for Action: >= 50)"
        },
        "target_therapies": [
            {"drug": "Diagnostic Optical Colonoscopy", "class": "Endoscopic Interception", "expected_response": "High Diagnostic Yield for Early Adenoma/Adenocarcinoma", "tier": "Immediate Next Step"},
            {"drug": "cfDNA Minimal Residual Disease Monitoring", "class": "Longitudinal Epigenetic Tracking", "expected_response": "Post-Resection ctDNA Clearance Verification", "tier": "Surveillance (Tier 2)"}
        ],
        "cpg_methylation_profile": [
            {"gene": "MGMT", "beta_value": 0.09, "status": "Unmethylated"},
            {"gene": "MLH1", "beta_value": 0.18, "status": "Low Methylation"},
            {"gene": "SEPT9", "beta_value": 0.76, "status": "Hypermethylated"},
            {"gene": "GSTP1", "beta_value": 0.08, "status": "Unmethylated"},
            {"gene": "CDKN2A", "beta_value": 0.32, "status": "Intermediate"},
            {"gene": "RASSF1A", "beta_value": 0.44, "status": "Intermediate"},
            {"gene": "BRCA1", "beta_value": 0.02, "status": "Unmethylated"},
            {"gene": "SHOX2", "beta_value": 0.11, "status": "Unmethylated"}
        ]
    },
    "ES_SMARCB1_EZH2_ONCOGENIC": {
        "id": "ES_SMARCB1_EZH2_ONCOGENIC",
        "title": "Epithelioid Sarcoma: SMARCB1/INI1 Loss & EZH2 Epigenetic Vulnerability",
        "cancer_type": "Advanced / Metastatic Epithelioid Sarcoma",
        "primary_gene": "SMARCB1 (SWI/SNF Subunit INI1) & EZH2 (PRC2 Catalytic Subunit)",
        "cpg_region": "Polycomb Repressive Complex 2 (PRC2) H3K27me3 Epigenetic Repression",
        "methylation_level_percent": 91.5,
        "assay_technique": "IHC (INI1 negative) + NGS / Histone Mass Spectrometry",
        "silencing_status": "Complete Loss of SWI/SNF Chromatin Remodeling Complex",
        "clinical_significance": "Biallelic loss of SMARCB1 removes natural antagonism against EZH2, leading to unopposed hyperactive PRC2 histone methyltransferase activity, global H3K27 trimethylation, and transcriptional silencing of differentiation genes. Directly sensitive to Tazemetostat.",
        "pathway_impact": "SWI/SNF vs. PRC2 epigenetic imbalance locks mesenchymal progenitor cells into an aggressive, poorly differentiated, proliferative oncogenic state.",
        "nccn_esmo_guideline": "FDA Accelerated Approval: Tazemetostat (800 mg orally twice daily) indicated for adults and pediatric patients aged >= 16 years with metastatic or locally advanced epithelioid sarcoma not eligible for complete resection.",
        "associated_markers": {
            "SMARCB1_Status": "Biallelic Loss / Deletion (0% IHC nuclear staining)",
            "H3K27me3_Level": "Markedly Elevated (Oncogenic PRC2 Hyper-repression)",
            "EZH2_Mutation": "Wild-Type (Epigenetically hyperactivated via SWI/SNF loss)",
            "TP53_Status": "Wild-Type"
        },
        "target_therapies": [
            {"drug": "Tazemetostat (Tazverik)", "class": "Selective EZH2 Inhibitor", "expected_response": "Clinically Proven Disease Control / Tumor Shrinkage", "tier": "FDA Approved (Tier 1A)"},
            {"drug": "HDAC Inhibitor (Vorinostat / Panobinostat)", "class": "Histone Deacetylase Inhibitor", "expected_response": "Synergistic Chromatin De-repression in Preclinical Models", "tier": "Clinical Trial (Tier 2B)"},
            {"drug": "Doxorubicin-Based Chemotherapy", "class": "Anthracycline Chemotherapy", "expected_response": "Modest Response / High Toxicity in Advanced Disease", "tier": "Alternative (Tier 2A)"}
        ],
        "cpg_methylation_profile": [
            {"gene": "MGMT", "beta_value": 0.05, "status": "Unmethylated"},
            {"gene": "MLH1", "beta_value": 0.04, "status": "Unmethylated"},
            {"gene": "SEPT9", "beta_value": 0.10, "status": "Unmethylated"},
            {"gene": "GSTP1", "beta_value": 0.09, "status": "Unmethylated"},
            {"gene": "CDKN2A", "beta_value": 0.92, "status": "Epigenetically Silenced (PRC2)"},
            {"gene": "RASSF1A", "beta_value": 0.38, "status": "Intermediate"},
            {"gene": "BRCA1", "beta_value": 0.03, "status": "Unmethylated"},
            {"gene": "SHOX2", "beta_value": 0.08, "status": "Unmethylated"}
        ]
    },
    "LGG_IDH1_GCIMP": {
        "id": "LGG_IDH1_GCIMP",
        "title": "Low-Grade Glioma: IDH1 R132H & G-CIMP Hypermethylation Phenotype",
        "cancer_type": "Astrocytoma / Oligodendroglioma (WHO Grade 2/3)",
        "primary_gene": "IDH1 (Isocitrate Dehydrogenase 1) & TET2 / KDM Epigenetic Enzymes",
        "cpg_region": "Genome-Wide CpG Islands (Glioma CpG Island Methylator Phenotype: G-CIMP)",
        "methylation_level_percent": 82.5,
        "assay_technique": "IDH1 R132H IHC / Targeted Sequencing + Methylation Array",
        "silencing_status": "G-CIMP Positive (Extensive CpG Island Hypermethylation)",
        "clinical_significance": "Neomorphic IDH1 mutation produces (R)-2-hydroxyglutarate (2-HG), competitively inhibiting alpha-ketoglutarate-dependent TET DNA demethylases and Jumonji histone demethylases. Confers favorable prognosis and sensitivity to Vorasidenib/Ivosidenib.",
        "pathway_impact": "2-HG oncometabolite accumulation causes genome-wide DNA and histone hypermethylation, arresting neural stem cell differentiation and driving low-grade gliomagenesis.",
        "nccn_esmo_guideline": "INDIGO Phase 3 Trial / FDA Approval: Vorasidenib (dual oral IDH1/IDH2 inhibitor with high blood-brain barrier penetrance) demonstrates statistically significant and clinically meaningful improvement in PFS and delay in subsequent intervention.",
        "associated_markers": {
            "IDH1_Mutation": "p.R132H (c.395G>A)",
            "1p19q_Codeletion": "Retained (Astrocytoma) or Codeleted (Oligodendroglioma)",
            "ATRX_Status": "Loss of Nuclear Expression (Mutated)",
            "2HG_Level": "Markedly Elevated (>5-fold in tumor tissue / MRS)"
        },
        "target_therapies": [
            {"drug": "Vorasidenib", "class": "Dual Brain-Penetrant IDH1/IDH2 Inhibitor", "expected_response": "Marked PFS Extension (PFS 27.7 mo vs 11.1 mo, HR 0.39)", "tier": "FDA Approved / NCCN Preferred (Tier 1A)"},
            {"drug": "Ivosidenib (Tibsovo)", "class": "Selective IDH1 Inhibitor", "expected_response": "Established Clinical Activity in IDH1-Mutant Malignancies", "tier": "Tier 1B"},
            {"drug": "Temozolomide (TMZ)", "class": "Alkylating Chemotherapy", "expected_response": "Adjuvant Chemotherapy in High-Risk / Progression", "tier": "Standard Adjuvant"}
        ],
        "cpg_methylation_profile": [
            {"gene": "MGMT", "beta_value": 0.81, "status": "Hypermethylated (G-CIMP)"},
            {"gene": "MLH1", "beta_value": 0.22, "status": "Low Methylation"},
            {"gene": "SEPT9", "beta_value": 0.35, "status": "Intermediate"},
            {"gene": "GSTP1", "beta_value": 0.42, "status": "Intermediate"},
            {"gene": "CDKN2A", "beta_value": 0.74, "status": "Hypermethylated"},
            {"gene": "RASSF1A", "beta_value": 0.79, "status": "Hypermethylated"},
            {"gene": "BRCA1", "beta_value": 0.15, "status": "Low Methylation"},
            {"gene": "SHOX2", "beta_value": 0.31, "status": "Intermediate"}
        ]
    }
}


class EpigenomicsEngine:
    """
    Precision Oncology Epigenetics & DNA Methylation Analysis Engine.
    Evaluates promoter hypermethylation, epigenetic silencing, tissue-of-origin deconvolution,
    and targeted epigenetic drug sensitivities.
    """

    def __init__(self):
        self.scenarios = CANONICAL_EPIGENETIC_SCENARIOS

    def get_scenario(self, scenario_id: str) -> Dict[str, Any]:
        """Retrieve canonical epigenetic scenario by ID."""
        return self.scenarios.get(scenario_id, self.scenarios["GBM_MGMT_METHYLATED"])

    def list_scenarios(self) -> List[Dict[str, str]]:
        """List summary of all available canonical epigenetic profiles."""
        return [
            {
                "id": sc["id"],
                "title": sc["title"],
                "cancer_type": sc["cancer_type"],
                "primary_gene": sc["primary_gene"]
            }
            for sc in self.scenarios.values()
        ]

    def calculate_epigenetic_metrics(self, scenario_id: str, custom_methylation_override: float = None) -> Dict[str, Any]:
        """
        Calculates actionable epigenetic metrics including:
        - Methylation Silencing Index (0 - 100%)
        - Epigenetic Targetability Score (0 - 100)
        - Alkylator / Temozolomide Sensitivity Classification
        - Lynch Syndrome Exclusion / Sporadic MSI Triage Confidence
        - Tissue-of-Origin (TOO) Deconvolution Score
        - NCCN/ESMO Therapeutic Alignment
        """
        scenario = self.get_scenario(scenario_id)
        
        meth_level = (
            float(custom_methylation_override)
            if custom_methylation_override is not None
            else scenario["methylation_level_percent"]
        )
        meth_level = max(0.0, min(100.0, meth_level))

        # Silencing threshold is clinically defined around 10-25% depending on gene/assay
        silencing_index = meth_level
        is_silenced = silencing_index >= 25.0

        # Calculate Epigenetic Targetability Score based on drug-target coupling
        if scenario_id == "GBM_MGMT_METHYLATED":
            # MGMT methylation directly scales TMZ sensitivity
            targetability_score = min(100.0, 30.0 + (meth_level * 0.70))
            drug_prediction = "High Alkylator (Temozolomide) Sensitivity" if is_silenced else "Temozolomide Resistance (Consider Clinical Trial / Lomustine)"
            lynch_triage = "N/A (Central Nervous System Primary)"
            too_deconvolution = "CNS / Glial Lineage (Confidence: 98.6%)"
        elif scenario_id == "CRC_MLH1_SPORADIC":
            # High MLH1 methylation -> sporadic CRC -> high IO sensitivity
            targetability_score = min(100.0, 40.0 + (meth_level * 0.60))
            drug_prediction = "High Anti-PD-1 (Pembrolizumab / Nivolumab+Ipilimumab) Sensitivity"
            lynch_triage = "Sporadic MSI-H Confirmed (Lynch Syndrome Ruled Out; Germline Testing Spared)" if is_silenced else "Indeterminate (Germline MLH1 Testing Required)"
            too_deconvolution = "Colorectal Mucosa (Confidence: 94.2%)"
        elif scenario_id == "MCED_LIQUID_BIOPSY_CRC":
            targetability_score = min(100.0, 20.0 + (meth_level * 0.75))
            drug_prediction = "Immediate Diagnostic Colonoscopy & Endoscopic Interception Indicated"
            lynch_triage = "Screening Marker (Evaluate post-histologic confirmation)"
            too_deconvolution = "Colorectal Tract (92.4%), Gastric (4.1%), Pancreatic (2.0%)"
        elif scenario_id == "ES_SMARCB1_EZH2_ONCOGENIC":
            targetability_score = 95.0
            drug_prediction = "Direct Sensitivity to EZH2 Inhibitor Tazemetostat (FDA Approved)"
            lynch_triage = "N/A (Mesenchymal Soft Tissue Sarcoma)"
            too_deconvolution = "Soft Tissue / Epitheloid Sarcoma Lineage (Confidence: 91.5%)"
        elif scenario_id == "LGG_IDH1_GCIMP":
            targetability_score = 92.0
            drug_prediction = "High Sensitivity to Brain-Penetrant Dual IDH1/2 Inhibitor Vorasidenib"
            lynch_triage = "N/A (IDH1-driven Gliomagenesis)"
            too_deconvolution = "Diffuse Astrocytoma / Oligodendroglioma (Confidence: 96.8%)"
        else:
            targetability_score = 50.0
            drug_prediction = "Standard Precision Epigenetics Profiling"
            lynch_triage = "Undetermined"
            too_deconvolution = "Unclassified Tissue"

        # Epigenetic Stage / Risk Stratification
        if meth_level >= 60.0:
            epigenetic_class = "Hypermethylated (High Epigenetic Divergence)"
            color_badge = "emerald"
        elif meth_level >= 25.0:
            epigenetic_class = "Intermediate Methylation (Moderate Silencing)"
            color_badge = "amber"
        else:
            epigenetic_class = "Hypomethylated / Unmethylated (Normal Expression Intact)"
            color_badge = "blue"

        return {
            "scenario_id": scenario_id,
            "title": scenario["title"],
            "cancer_type": scenario["cancer_type"],
            "primary_gene": scenario["primary_gene"],
            "methylation_level_percent": round(meth_level, 1),
            "silencing_status": "Transcriptional Gene Silencing" if is_silenced else "Active Gene Expression / Unmethylated",
            "epigenetic_class": epigenetic_class,
            "color_badge": color_badge,
            "targetability_score": round(targetability_score, 1),
            "predicted_therapeutic_response": drug_prediction,
            "lynch_syndrome_triage": lynch_triage,
            "tissue_of_origin_deconvolution": too_deconvolution,
            "guideline": scenario["nccn_esmo_guideline"],
            "target_therapies": scenario["target_therapies"],
            "cpg_profile": scenario["cpg_methylation_profile"]
        }


def analyze_epigenetics(scenario_id: str, custom_methylation_override: float = None) -> Dict[str, Any]:
    """Convenience helper function for rapid API-level analysis."""
    engine = EpigenomicsEngine()
    return engine.calculate_epigenetic_metrics(scenario_id, custom_methylation_override)
