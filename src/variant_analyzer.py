"""
Genomic Variant Analyzer Module
Inspects genetic variants across ClinVar, dbSNP, and gnomAD datasets.
"""

import os
import json

class VariantAnalyzer:
    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir

    def load_clinvar_summary(self, filename: str = "summary_all.json"):
        path = os.path.join(self.data_dir, filename)
        if not os.path.exists(path):
            return []
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    def analyze_pathogenicity(self, filename: str = "summary_all.json"):
        variants = self.load_clinvar_summary(filename)
        pathogenic = []
        benign = []
        vus = []

        for v in variants:
            sig = v.get("clinical_significance", "").lower()
            item = {
                "variant_id": v.get("variant_id"),
                "title": v.get("title"),
                "clinical_significance": v.get("clinical_significance"),
                "review_status": v.get("review_status")
            }
            if "pathogenic" in sig:
                pathogenic.append(item)
            elif "benign" in sig:
                benign.append(item)
            else:
                vus.append(item)

        return {
            "total_variants": len(variants),
            "pathogenic_count": len(pathogenic),
            "benign_count": len(benign),
            "vus_other_count": len(vus),
            "pathogenic_samples": pathogenic[:3]
        }

if __name__ == "__main__":
    analyzer = VariantAnalyzer()
    print("Variant Pathogenicity Breakdown:")
    print(json.dumps(analyzer.analyze_pathogenicity(), indent=2, ensure_ascii=False))
