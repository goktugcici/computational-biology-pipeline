"""
Pan-Cancer Oncology Profiler Module
Analyzes high-impact oncogenic driver mutations (TP53, EGFR, KRAS, BRAF, BRCA1/2),
clinical actionability, and precision targeted therapeutics.
"""

import os
import json
from typing import Dict, List, Any

class OncologyProfiler:
    def __init__(self):
        # Curated precision oncology knowledgebase
        self.oncology_db = {
            "EGFR": {
                "full_name": "Epidermal Growth Factor Receptor",
                "cancer_types": ["Non-Small Cell Lung Cancer (NSCLC)", "Glioblastoma"],
                "hotspot_mutations": [
                    {"variant": "p.L858R", "exon": 21, "type": "Activating Missense", "response": "Sensitive to 1st/3rd gen TKIs"},
                    {"variant": "p.T790M", "exon": 20, "type": "Acquired Resistance", "response": "Sensitive to Osimertinib (3rd gen)"},
                    {"variant": "Exon 19 del", "exon": 19, "type": "In-frame Deletion", "response": "Highly sensitive to EGFR TKIs"}
                ],
                "fda_targeted_therapies": [
                    {"drug": "Osimertinib (Tagrisso)", "generation": "3rd Gen TKI", "line": "1st line metastatic NSCLC"},
                    {"drug": "Gefitinib (Iressa)", "generation": "1st Gen TKI", "line": "Standard TKI"},
                    {"drug": "Afatinib (Gilotrif)", "generation": "2nd Gen TKI", "line": "Irreversible ErbB Family Blocker"}
                ],
                "pdb_reference": "2J6M"
            },
            "KRAS": {
                "full_name": "KRAS Proto-Oncogene, GTPase",
                "cancer_types": ["Colorectal Cancer", "Pancreatic Ductal Adenocarcinoma", "NSCLC"],
                "hotspot_mutations": [
                    {"variant": "p.G12C", "exon": 2, "type": "Constitutive GTPase lock", "response": "Covalent switch-II pocket targetable"},
                    {"variant": "p.G12D", "exon": 2, "type": "Constitutive Activation", "response": "Emerging small molecule inhibitors"},
                    {"variant": "p.G12V", "exon": 2, "type": "GTPase Inactivation", "response": "MAPK/ERK pathway hyperactivation"}
                ],
                "fda_targeted_therapies": [
                    {"drug": "Sotorasib (Lumakras)", "target": "KRAS G12C", "line": "FDA Approved for G12C NSCLC"},
                    {"drug": "Adagrasib (Krazati)", "target": "KRAS G12C", "line": "FDA Approved for G12C NSCLC"}
                ],
                "pdb_reference": "6OIM"
            },
            "BRAF": {
                "full_name": "B-Raf Proto-Oncogene, Serine/Threonine Kinase",
                "cancer_types": ["Cutaneous Melanoma", "Colorectal Cancer", "Thyroid Carcinoma"],
                "hotspot_mutations": [
                    {"variant": "p.V600E", "exon": 15, "type": "Kinase domain phosphomimetic", "response": "Sensitive to BRAF+MEK inhibition"},
                    {"variant": "p.V600K", "exon": 15, "type": "Kinase Activation", "response": "Sensitive to Vemurafenib/Dabrafenib"}
                ],
                "fda_targeted_therapies": [
                    {"drug": "Dabrafenib (Tafinlar) + Trametinib (Mekinist)", "combination": "BRAF + MEK inhibitor", "line": "1st line BRAF-mutant melanoma"},
                    {"drug": "Vemurafenib (Zelboraf)", "combination": "Single agent BRAF inhibitor", "line": "Metastatic melanoma"}
                ],
                "pdb_reference": "4MNE"
            },
            "TP53": {
                "full_name": "Tumor Protein P53 (Guardian of the Genome)",
                "cancer_types": ["Pan-Cancer (>50% of all human malignancies)", "Li-Fraumeni Syndrome"],
                "hotspot_mutations": [
                    {"variant": "p.R175H", "domain": "DNA-binding core", "type": "Structural instability", "response": "Zinc coordination disruption"},
                    {"variant": "p.R248Q", "domain": "DNA-contact residue", "type": "Contact mutation", "response": "Loss of minor groove binding"},
                    {"variant": "p.R273H", "domain": "DNA-contact residue", "type": "Contact mutation", "response": "Loss of phosphate backbone contact"}
                ],
                "fda_targeted_therapies": [
                    {"drug": "WEE1/ATR Inhibitors (Adavosertib)", "status": "Clinical Trials", "mechanism": "Synthetic lethality in p53-deficient tumors"},
                    {"drug": "Arsenic Trioxide + Rescue chaperones", "status": "Investigational", "mechanism": "p53 structural refolding"}
                ],
                "pdb_reference": "1TSR"
            }
        }

    def list_supported_genes(self) -> List[str]:
        return list(self.oncology_db.keys())

    def profile_gene(self, gene: str) -> Dict[str, Any]:
        gene_upper = gene.upper()
        if gene_upper not in self.oncology_db:
            return {"error": f"Gene '{gene}' not found. Supported genes: {', '.join(self.list_supported_genes())}"}
        
        info = self.oncology_db[gene_upper]
        return {
            "gene": gene_upper,
            "full_name": info["full_name"],
            "cancer_types": info["cancer_types"],
            "hotspot_count": len(info["hotspot_mutations"]),
            "hotspots": info["hotspot_mutations"],
            "targeted_therapies": info["fda_targeted_therapies"],
            "recommended_pdb": info["pdb_reference"]
        }

    def pan_cancer_summary(self) -> Dict[str, Any]:
        total_hotspots = sum(len(g["hotspot_mutations"]) for g in self.oncology_db.values())
        total_therapies = sum(len(g["fda_targeted_therapies"]) for g in self.oncology_db.values())
        return {
            "genes_covered": self.list_supported_genes(),
            "total_target_genes": len(self.oncology_db),
            "total_hotspots": total_hotspots,
            "total_therapies_cataloged": total_therapies
        }

if __name__ == "__main__":
    profiler = OncologyProfiler()
    print("Pan-Cancer Overview:")
    print(json.dumps(profiler.pan_cancer_summary(), indent=2, ensure_ascii=False))
    print("\nEGFR Precision Profile:")
    print(json.dumps(profiler.profile_gene("EGFR"), indent=2, ensure_ascii=False))
