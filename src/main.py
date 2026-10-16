"""
Computational Biology & Science Pipeline - Main Runner
Executes an end-to-end multi-omics analysis spanning structural biology,
genomics, pharmacology, biological pathways, and literature intelligence.
Includes both local curated dataset analysis, live dynamic API querying,
pan-cancer oncology profiling, and publication-ready HTML report generation.
"""

import sys
import json
import argparse
from structural_biology import StructureAnalyzer
from variant_analyzer import VariantAnalyzer
from pharmacology import PharmacologyAnalyzer
from pathway_mapper import PathwayMapper
from literature_miner import LiteratureMiner
from live_fetcher import LiveFetcher
from oncology_profiler import OncologyProfiler
from report_generator import ReportGenerator

def run_pipeline(data_dir: str = "data", live_pdb: str = None, live_drug: str = None, 
                 live_pubmed: str = None, pan_cancer: str = None, generate_report: bool = False):
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

    # OPTIONAL: Dynamic Live Fetcher Extension
    if live_pdb or live_drug or live_pubmed:
        print("\n" + "-" * 70)
        print("LIVE SCIENTIFIC API QUERIES (Dynamic Search Engine)")
        print("-" * 70)
        fetcher = LiveFetcher(f"{data_dir}/cache")

        if live_pdb:
            print(f"\n[LIVE PDB] Querying: {live_pdb}")
            pdb_info = fetcher.fetch_pdb_entry(live_pdb)
            from_cache = " (Cached)" if pdb_info.get("_from_cache") else " (Live API)"
            print(f"  Title: {pdb_info.get('title')}{from_cache}")
            print(f"  Resolution: {pdb_info.get('resolution')} Å | Method: {pdb_info.get('method')}")

        if live_drug:
            print(f"\n[LIVE PUBCHEM] Querying: {live_drug}")
            drug_info = fetcher.fetch_pubchem_drug(live_drug)
            from_cache = " (Cached)" if drug_info.get("_from_cache") else " (Live API)"
            print(f"  CID: {drug_info.get('cid')} | Formula: {drug_info.get('formula')}{from_cache}")
            print(f"  Molecular Weight: {drug_info.get('molecular_weight')} g/mol | XLogP: {drug_info.get('xlogp')}")

        if live_pubmed:
            print(f"\n[LIVE PUBMED] Searching: '{live_pubmed}'")
            pm_info = fetcher.fetch_pubmed_articles(live_pubmed, max_results=2)
            articles = pm_info.get("articles", [])
            for a in articles:
                print(f"   * [PMID: {a.get('pmid')}] {a.get('title')} ({a.get('source')})")

    # OPTIONAL: Pan-Cancer Oncology Profiler
    if pan_cancer:
        print("\n" + "-" * 70)
        print("PAN-CANCER ONCOLOGY & TARGETED THERAPEUTICS")
        print("-" * 70)
        profiler = OncologyProfiler()
        profile = profiler.profile_gene(pan_cancer)
        if "error" in profile:
            print(f"  Error: {profile['error']}")
        else:
            print(f"  Target Gene: {profile['gene']} ({profile['full_name']})")
            print(f"  Associated Malignancies: {', '.join(profile['cancer_types'])}")
            print(f"  Cataloged Hotspot Mutations ({profile['hotspot_count']}):")
            for h in profile['hotspots']:
                print(f"    * {h['variant']} [{h.get('type')}]: {h.get('response')}")
            print("  FDA-Approved Targeted Therapies:")
            for d in profile['targeted_therapies']:
                print(f"    -> {d['drug']} ({d.get('line', d.get('status'))})")

    # OPTIONAL: Publication HTML Report Generation
    if generate_report:
        print("\n" + "-" * 70)
        print("PUBLICATION REPORT GENERATION")
        print("-" * 70)
        rg = ReportGenerator(data_dir=data_dir, output_dir="reports")
        rep_file = rg.generate_html_report()
        print(f"  Executive Report Generated: {rep_file}")

    print("\n" + "=" * 70)
    print("PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
    print("=" * 70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Computational Biology Pipeline")
    parser.add_argument("--data_dir", default="data", help="Path to data directory")
    parser.add_argument("--live_pdb", default=None, help="Fetch live PDB metadata (e.g. 1TSR, 6LU7)")
    parser.add_argument("--live_drug", default=None, help="Fetch live PubChem drug properties (e.g. Imatinib, Aspirin)")
    parser.add_argument("--live_pubmed", default=None, help="Search live PubMed literature (e.g. 'CAR-T lymphoma')")
    parser.add_argument("--pan_cancer", default=None, help="Profile oncogene (e.g. EGFR, KRAS, BRAF, TP53)")
    parser.add_argument("--generate_report", action="store_true", help="Generate publication-grade HTML report")
    args = parser.parse_args()
    run_pipeline(args.data_dir, args.live_pdb, args.live_drug, args.live_pubmed, args.pan_cancer, args.generate_report)
