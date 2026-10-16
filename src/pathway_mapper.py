"""
Biological Pathway Mapper Module
Extracts cellular pathways and functional ontology terms from Reactome database records.
"""

import os
import json

class PathwayMapper:
    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir

    def get_gene_pathways(self, gene: str = "hbb"):
        path = os.path.join(self.data_dir, f"reactome_{gene.lower()}.json")
        if not os.path.exists(path):
            return []

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            pathways = data.get("pathways", [])
            results = []
            for p in pathways:
                results.append({
                    "stId": p.get("stId"),
                    "name": p.get("name"),
                    "species": p.get("species", {}).get("name") if isinstance(p.get("species"), dict) else "Homo sapiens"
                })
            return results

if __name__ == "__main__":
    pm = PathwayMapper()
    print("Top Pathways for HBB:")
    print(json.dumps(pm.get_gene_pathways("hbb")[:6], indent=2, ensure_ascii=False))
