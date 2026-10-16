"""
Live Fetcher & Dynamic Scientific API Engine
Queries live biological databases (RCSB PDB, NCBI ClinVar, PubChem, PubMed)
with an automatic persistent caching layer.
"""

import os
import json
import urllib.request
import urllib.parse
from typing import Dict, Any, Optional

class LiveFetcher:
    def __init__(self, cache_dir: str = "data/cache"):
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

    def _get_cache_path(self, namespace: str, key: str) -> str:
        safe_key = "".join([c if c.isalnum() or c in "-_" else "_" for c in key])
        return os.path.join(self.cache_dir, f"{namespace}_{safe_key}.json")

    def _read_cache(self, namespace: str, key: str) -> Optional[Dict[str, Any]]:
        path = self._get_cache_path(namespace, key)
        if os.path.exists(path):
            try:
                with open(path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                return None
        return None

    def _write_cache(self, namespace: str, key: str, data: Dict[str, Any]):
        path = self._get_cache_path(namespace, key)
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
        except Exception:
            pass

    def fetch_pdb_entry(self, pdb_id: str) -> Dict[str, Any]:
        """Queries RCSB PDB Data API for structural metadata."""
        pdb_id = pdb_id.upper()
        cached = self._read_cache("pdb", pdb_id)
        if cached:
            cached["_from_cache"] = True
            return cached

        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        req = urllib.request.Request(url, headers={"User-Agent": "ComputationalBioPipeline/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=12) as response:
                data = json.loads(response.read().decode("utf-8"))
                result = {
                    "pdb_id": pdb_id,
                    "title": data.get("struct", {}).get("title"),
                    "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                    "method": data.get("exptl", [{}])[0].get("method"),
                    "polymer_count": data.get("rcsb_entry_info", {}).get("polymer_entity_count_protein"),
                    "deposition_date": data.get("rcsb_accession_info", {}).get("deposit_date")
                }
                self._write_cache("pdb", pdb_id, result)
                result["_from_cache"] = False
                return result
        except Exception as e:
            return {"error": f"Failed to fetch PDB {pdb_id}: {str(e)}"}

    def fetch_pubchem_drug(self, drug_name: str) -> Dict[str, Any]:
        """Queries PubChem PUG-REST for chemical properties."""
        cached = self._read_cache("pubchem", drug_name.lower())
        if cached:
            cached["_from_cache"] = True
            return cached

        safe_name = urllib.parse.quote(drug_name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe_name}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,XLogP,TPSA/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": "ComputationalBioPipeline/1.0"})
        try:
            with urllib.request.urlopen(req, timeout=12) as response:
                data = json.loads(response.read().decode("utf-8"))
                props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                result = {
                    "name": drug_name.capitalize(),
                    "cid": props.get("CID"),
                    "formula": props.get("MolecularFormula"),
                    "molecular_weight": props.get("MolecularWeight"),
                    "smiles": props.get("CanonicalSMILES"),
                    "xlogp": props.get("XLogP"),
                    "tpsa": props.get("TPSA")
                }
                self._write_cache("pubchem", drug_name.lower(), result)
                result["_from_cache"] = False
                return result
        except Exception as e:
            return {"error": f"Failed to fetch PubChem compound {drug_name}: {str(e)}"}

    def fetch_pubmed_articles(self, term: str, max_results: int = 3) -> Dict[str, Any]:
        """Queries PubMed for articles matching a biomedical term."""
        cached = self._read_cache("pubmed", term.lower())
        if cached:
            cached["_from_cache"] = True
            return cached

        safe_term = urllib.parse.quote(term)
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={safe_term}&retmax={max_results}&retmode=json"
        try:
            with urllib.request.urlopen(search_url, timeout=12) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                id_list = data.get("esearchresult", {}).get("idlist", [])

            articles = []
            if id_list:
                ids_str = ",".join(id_list)
                summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={ids_str}&retmode=json"
                with urllib.request.urlopen(summary_url, timeout=12) as s_resp:
                    s_data = json.loads(s_resp.read().decode("utf-8"))
                    res_dict = s_data.get("result", {})
                    for pmid in id_list:
                        item = res_dict.get(pmid, {})
                        articles.append({
                            "pmid": pmid,
                            "title": item.get("title"),
                            "source": item.get("source"),
                            "pubdate": item.get("pubdate")
                        })

            result = {"query": term, "count": len(articles), "articles": articles}
            self._write_cache("pubmed", term.lower(), result)
            result["_from_cache"] = False
            return result
        except Exception as e:
            return {"error": f"Failed to search PubMed: {str(e)}"}

if __name__ == "__main__":
    fetcher = LiveFetcher()
    print("Testing Live PDB Fetch (1TSR - Tumor Suppressor p53):")
    print(json.dumps(fetcher.fetch_pdb_entry("1TSR"), indent=2, ensure_ascii=False))

    print("\nTesting Live PubChem Fetch (Imatinib - Gleevec):")
    print(json.dumps(fetcher.fetch_pubchem_drug("Imatinib"), indent=2, ensure_ascii=False))
