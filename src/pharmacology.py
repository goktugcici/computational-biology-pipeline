"""
Pharmacology and Cheminformatics Module
Processes chemical properties, SMILES representations, and drug mechanisms from PubChem & ChEMBL.
"""

import os
import json

class PharmacologyAnalyzer:
    def __init__(self, data_dir: str = "data"):
        self.data_dir = data_dir

    def get_compound_profile(self, name: str = "voxelotor"):
        props_path = os.path.join(self.data_dir, f"{name}_props.json")
        pharm_path = os.path.join(self.data_dir, f"{name}_pharm.json")

        profile = {"name": name.capitalize(), "properties": {}, "mechanism": "Not available"}

        if os.path.exists(props_path):
            with open(props_path, "r", encoding="utf-8") as f:
                p_data = json.load(f)
                props_list = p_data.get("PropertyTable", {}).get("Properties", [])
                if props_list:
                    profile["properties"] = props_list[0]

        if os.path.exists(pharm_path):
            with open(pharm_path, "r", encoding="utf-8") as f:
                ph_data = json.load(f)
                try:
                    sections = ph_data['Record']['Section'][0]['Section']
                    for s in sections:
                        if s.get('TOCHeading') == 'Mechanism of Action':
                            profile["mechanism"] = s['Information'][0]['Value']['StringWithMarkup'][0]['String']
                except Exception:
                    pass

        return profile

    def get_chembl_molecule_summary(self, filename: str = "chembl_metformin.json"):
        path = os.path.join(self.data_dir, filename)
        if not os.path.exists(path):
            return {}
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
            mols = data.get("molecules", [])
            if mols:
                m = mols[0]
                return {
                    "chembl_id": m.get("molecule_chembl_id"),
                    "name": m.get("pref_name"),
                    "type": m.get("molecule_type"),
                    "max_phase": m.get("max_phase"),
                    "formula": m.get("molecule_properties", {}).get("full_molformula"),
                    "mw": m.get("molecule_properties", {}).get("full_mwt")
                }
        return {}

if __name__ == "__main__":
    ph = PharmacologyAnalyzer()
    print("Voxelotor Profile:")
    print(json.dumps(ph.get_compound_profile("voxelotor"), indent=2, ensure_ascii=False))
    print("\nChEMBL Molecule Profile:")
    print(json.dumps(ph.get_chembl_molecule_summary(), indent=2, ensure_ascii=False))
