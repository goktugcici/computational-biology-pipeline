"""
Unit Tests for Computational Biology Pipeline
Tests structural, variant, pharmacology, pathway, and literature modules.
"""

import os
import unittest
import sys

# src dizinini path'e ekle
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from structural_biology import StructureAnalyzer
from variant_analyzer import VariantAnalyzer
from pharmacology import PharmacologyAnalyzer
from pathway_mapper import PathwayMapper
from literature_miner import LiteratureMiner

class TestComputationalBiologyPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

    def test_structural_analyzer(self):
        cif_path = os.path.join(self.data_dir, "pdb_00004hhb.cif")
        analyzer = StructureAnalyzer(cif_path)
        summary = analyzer.get_summary()
        self.assertIn("file", summary)
        self.assertEqual(summary["file"], "pdb_00004hhb.cif")
        self.assertGreater(summary["size_bytes"], 0)

    def test_variant_analyzer(self):
        analyzer = VariantAnalyzer(self.data_dir)
        res = analyzer.analyze_pathogenicity()
        self.assertGreaterEqual(res["total_variants"], 1)
        self.assertGreaterEqual(res["pathogenic_count"], 1)
        self.assertIsInstance(res["pathogenic_samples"], list)

    def test_pharmacology_analyzer(self):
        analyzer = PharmacologyAnalyzer(self.data_dir)
        vox = analyzer.get_compound_profile("voxelotor")
        self.assertEqual(vox["name"], "Voxelotor")
        props = vox.get("properties", {})
        self.assertIn("MolecularFormula", props)
        self.assertEqual(props.get("MolecularFormula"), "C19H19N3O3")

        # ChEMBL Metformin check
        chembl = analyzer.get_chembl_molecule_summary()
        self.assertEqual(chembl.get("chembl_id"), "CHEMBL1431")
        self.assertEqual(chembl.get("name"), "METFORMIN")

    def test_pathway_mapper(self):
        mapper = PathwayMapper(self.data_dir)
        pathways = mapper.get_gene_pathways("hbb")
        self.assertIsInstance(pathways, list)
        self.assertGreater(len(pathways), 0)
        st_ids = [p["stId"] for p in pathways]
        self.assertTrue(any("R-HSA" in sid for sid in st_ids))

    def test_literature_miner(self):
        miner = LiteratureMiner(self.data_dir)
        lit = miner.get_abstract_summary()
        self.assertEqual(lit.get("pmid"), "38661449")
        self.assertIn("Exagamglogene", lit.get("title", ""))

if __name__ == "__main__":
    unittest.main()
