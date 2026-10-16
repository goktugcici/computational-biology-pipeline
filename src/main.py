"""
Computational Biology & Science Pipeline - Main Runner
Executes an end-to-end multi-omics analysis spanning structural biology,
genomics, pharmacology, biological pathways, and literature intelligence.
"""

import sys
import json
import argparse
from structural_biology import StructureAnalyzer
from variant_analyzer import VariantAnalyzer
from pharmacology import PharmacologyAnalyzer
from pathway_mapper import PathwayMapper
from literature_miner import LiteratureMiner

def run_pipeline(data_dir: str = "data"):
    print("=" * 70)
    print("COMPUTATIONAL BIOLOGY & SCIENCE PIPELINE")
    print("Multi-Omics Structural, Genomic, Pharmacological & Literature Analysis")
    print("=" * 70)

    # 1. Structural Biology
    print("\n[1/5] 3D STRUCTURAL BIOLOGY (Hemoglobin 4HHB)")
    struct = StructureAnalyzer(f"{data_dir}/pdb_00004hhb.cif")
    print(f"  Loaded file: {struct.file_path}")
    print(f"  Detected ligands/co-factors: {struct.ligands}")

    # 2. Genomic Variant Analysis
    print("\n[2/5] GENOMIC VARIANT ANALYSIS (HBB rs334 & ClinVar Pathogenicity)")
    var_analyzer = VariantAnalyzer(data_dir)
    res = var_analyzer.analyze_pathogenicity()
    print(f"  Total variants screened: {res['total_variants']}")
    print(f"  Pathogenic: {res['pathogenic_count']} | Benign: {res['benign_count']} | VUS/Other: {res['vus_other_count']}")
    if res['pathogenic_samples']:
        top_p = res['pathogenic_samples'][0]
        print(f"  Key Pathogenic Mutation: {top_p['title']} [{top_p['clinical_significance']}]")

    # 3. Pharmacology & Cheminformatics
    print("\n[3/5] PHARMACOLOGY & DRUG DISCOVERY (Voxelotor & Metformin)")
    pharm = PharmacologyAnalyzer(data_dir)
    vox = pharm.get_compound_profile("voxelotor")
    props = vox.get("properties", {})
    print(f"  Drug: {vox['name']} (Formula: {props.get('MolecularFormula')}, MW: {props.get('MolecularWeight')} g/mol)")
    print(f"  Target Mechanism: {vox.get('mechanism')[:120]}...")
    met = pharm.get_chembl_molecule_summary()
    if met:
        print(f"  ChEMBL Molecule: {met.get('name')} ({met.get('chembl_id')}) - Type: {met.get('type')}")

    # 4. Biological Pathways
    print("\n[4/5] BIOLOGICAL PATHWAYS (Reactome HBB Mapping)")
    pm = PathwayMapper(data_dir)
    pathways = pm.get_gene_pathways("hbb")
    print(f"  Mapped Pathways: {len(pathways)} cellular events identified")
    for p in pathways[:3]:
        print(f"   * [{p['stId']}] {p['name']}")

    # 5. Scientific Literature Intelligence
    print("\n[5/5] SCIENTIFIC LITERATURE & CRISPR BREAKTHROUGH (PubMed)")
    lm = LiteratureMiner(data_dir)
    lit = lm.get_abstract_summary()
    print(f"  Landmark Study (PMID: {lit.get('pmid')}): {lit.get('title')}")
    print(f"  Journal: {lit.get('journal')}")

    print("\n" + "=" * 70)
    print("PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
    print("=" * 70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Computational Biology Pipeline")
    parser.add_argument("--data_dir", default="data", help="Path to data directory")
    args = parser.parse_args()
    run_pipeline(args.data_dir)
