"""
Literature Intelligence Module
Mines biomedical publications, clinical trials, and breakthrough articles from PubMed / PMC records.
"""

import os
import json

class LiteratureMiner:
    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir

    def get_abstract_summary(self, filename: str = "pubmed_abstract.json"):
        path = os.path.join(self.data_dir, filename)
        if not os.path.exists(path):
            return {}

        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            if isinstance(data, list) and len(data) > 0:
                item = data[0]
                return {
                    "pmid": item.get("pmid"),
                    "title": item.get("title"),
                    "journal": item.get("journal"),
                    "pubdate": item.get("pubdate"),
                    "authors": item.get("authors", [])[:4],
                    "abstract_snippet": (item.get("abstract") or "")[:450] + "..."
                }
        return {}

if __name__ == "__main__":
    lm = LiteratureMiner()
    print("Literature Insight:")
    print(json.dumps(lm.get_abstract_summary(), indent=2, ensure_ascii=False))
