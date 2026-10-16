"""
Computational Biology & Science Pipeline - Main Runner (v2.0)
Executes an end-to-end multi-omics analysis spanning structural biology,
genomics, pharmacology, biological pathways, literature intelligence,
pairwise sequence alignments, in silico drug screening, and PPI interactomes.
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
from sequence_aligner import SequenceAligner, KNOWN_VARIANTS
from docking_scorer import DockingScorer
from interaction_network import InteractionNetwork
from report_generator import ReportGenerator

def run_pipeline(data_dir: str = "data", live_pdb: str = None, live_drug: str = None, 
                 live_pubmed: str = None, pan_cancer: str = None, align_variant: str = None,
                 screen_drugs_target: str = None, ppi_gene: str = None,
                 export_json_path: str = None, generate_report: bool = False):
    print("=" * 70)
    print("COMPUTATIONAL BIOLOGY & SCIENCE PIPELINE v2.0")
    print("Multi-Omics Structural, Genomic, Pharmacological & Network Analysis")
    print("=" * 70)

    pipeline_summary = {}

    # 1. Structural Biology
    print("\n[1/7] 3D STRUCTURAL BIOLOGY (Hemoglobin 4HHB)")
    struct = StructureAnalyzer(f"{data_dir}/pdb_00004hhb.cif")
    print(f"  Loaded file: {struct.file_path}")
    print(f"  Detected ligands/co-factors: {struct.ligands[:10]}...")
    pipeline_summary['structural_biology'] = {'file': struct.file_path, 'ligand_count': len(struct.ligands)}

    # 2. Genomic Variant Analysis
    print("\n[2/7] GENOMIC VARIANT ANALYSIS (HBB rs334 & ClinVar Pathogenicity)")
    var_analyzer = VariantAnalyzer(data_dir)
    res = var_analyzer.analyze_pathogenicity()
    print(f"  Total variants screened: {res['total_variants']}")
    print(f"  Pathogenic: {res['pathogenic_count']} | Benign: {res['benign_count']} | VUS/Other: {res['vus_other_count']}")
    if res['pathogenic_samples']:
        top_p = res['pathogenic_samples'][0]
        print(f"  Key Pathogenic Mutation: {top_p['title']} [{top_p['clinical_significance']}]")
    pipeline_summary['genomic_variants'] = res

    # 3. Pharmacology & Cheminformatics
    print("\n[3/7] PHARMACOLOGY & DRUG DISCOVERY (Voxelotor & Metformin)")
    pharm = PharmacologyAnalyzer(data_dir)
    vox = pharm.get_compound_profile("voxelotor")
    props = vox.get("properties", {})
    print(f"  Drug: {vox['name']} (Formula: {props.get('MolecularFormula')}, MW: {props.get('MolecularWeight')} g/mol)")
    print(f"  Target Mechanism: {vox.get('mechanism')[:110]}...")
    met = pharm.get_chembl_molecule_summary()
    if met:
        print(f"  ChEMBL Molecule: {met.get('name')} ({met.get('chembl_id')}) - Type: {met.get('type')}")
    pipeline_summary['pharmacology'] = {'voxelotor': props, 'chembl': met}

    # 4. Biological Pathways
    print("\n[4/7] BIOLOGICAL PATHWAYS (Reactome HBB Mapping)")
    pm = PathwayMapper(data_dir)
    pathways = pm.get_gene_pathways("hbb")
    print(f"  Mapped Pathways: {len(pathways)} cellular events identified")
    for p in pathways[:3]:
        print(f"   * [{p['stId']}] {p['name']}")
    pipeline_summary['pathways_count'] = len(pathways)

    # 5. Scientific Literature Intelligence
    print("\n[5/7] SCIENTIFIC LITERATURE & CRISPR BREAKTHROUGH (PubMed)")
    lm = LiteratureMiner(data_dir)
    lit = lm.get_abstract_summary()
    print(f"  Landmark Study (PMID: {lit.get('pmid')}): {lit.get('title')}")
    print(f"  Journal: {lit.get('journal')}")
    pipeline_summary['literature'] = {'pmid': lit.get('pmid'), 'title': lit.get('title')}

    # 6. Pairwise Sequence Alignment & Mutation Consequence
    print("\n[6/7] SEQUENCE ALIGNMENT & BLOSUM62 CONSERVATION")
    aligner = SequenceAligner()
    chosen_variant = align_variant if align_variant else 'HBB_SICKLE'
    if chosen_variant in KNOWN_VARIANTS:
        al_res = aligner.analyze_known_variant(chosen_variant)
        print(f"  Variant Target: {chosen_variant} ({al_res['metadata']['condition']})")
        print(f"  Identity: {al_res['identity_percent']}% | Score: {al_res['alignment_score']}")
        if al_res['mutations_detected']:
            mut = al_res['mutations_detected'][0]
            print(f"  Detected Mutation: {mut['description']} ({mut['type']}, BLOSUM62: {mut['blosum62_score']})")
        pipeline_summary['sequence_alignment'] = al_res
    else:
        print(f"  Unknown variant key '{chosen_variant}'. Available: {list(KNOWN_VARIANTS.keys())}")

    # 7. In Silico Drug Screening & Lipinski Scoring
    print("\n[7/7] IN SILICO DRUG SCREENING & LIPINSKI PROFILING")
    docking = DockingScorer()
    screen_target = screen_drugs_target if screen_drugs_target else (pan_cancer if pan_cancer else 'EGFR')
    ranked_drugs = docking.screen_library_against_target(screen_target)
    print(f"  Target Pocket: {ranked_drugs[0]['target_pocket']}")
    print(f"  Top In Silico Hit: #{ranked_drugs[0]['screen_rank']} {ranked_drugs[0]['drug_name']}")
    print(f"  Estimated dG: {ranked_drugs[0]['estimated_delta_g_kcal_mol']} kcal/mol (Fitness: {ranked_drugs[0]['contact_fitness_score']}%)")
    pipeline_summary['drug_screening'] = ranked_drugs[:5]

    # OPTIONAL: PPI Network Analysis
    if ppi_gene:
        print("\n" + "-" * 70)
        print("PROTEIN-PROTEIN INTERACTION (PPI) NETWORK TOPOLOGY")
        print("-" * 70)
        ppi = InteractionNetwork(ppi_gene)
        ppi_data = ppi.compute_network_metrics()
        print(f"  Network: {ppi_data['network_name']}")
        print(f"  Nodes: {ppi_data['total_nodes']} | Edges: {ppi_data['total_edges']} | Density: {ppi_data['network_density']}")
        print(f"  Primary Bottleneck Hub: {ppi_data['primary_hub_protein']} | Centrality Hubs: {', '.join(ppi_data['top_ranked_hubs'])}")
        pipeline_summary['ppi_network'] = ppi_data

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
            print(f"  Resolution: {pdb_info.get('resolution')} A | Method: {pdb_info.get('method')}")

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
        pipeline_summary['pan_cancer_profile'] = profile

    # OPTIONAL: Export JSON
    if export_json_path:
        with open(export_json_path, 'w', encoding='utf-8') as f:
            json.dump(pipeline_summary, f, indent=2, ensure_ascii=False)
        print(f"\n  JSON Export saved to: {export_json_path}")

    # OPTIONAL: Publication HTML Report Generation
    if generate_report:
        print("\n" + "-" * 70)
        print("PUBLICATION REPORT GENERATION")
        print("-" * 70)
        rg = ReportGenerator(data_dir=data_dir, output_dir="reports")
        rep_file = rg.generate_html_report(target_gene=screen_target)
        print(f"  Executive Report Generated: {rep_file}")

    print("\n" + "=" * 70)
    print("PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
    print("=" * 70)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Computational Biology Pipeline v2.0")
    parser.add_argument("--data_dir", default="data", help="Path to data directory")
    parser.add_argument("--live_pdb", default=None, help="Fetch live PDB metadata (e.g. 1TSR, 6LU7)")
    parser.add_argument("--live_drug", default=None, help="Fetch live PubChem drug properties (e.g. Imatinib, Aspirin)")
    parser.add_argument("--live_pubmed", default=None, help="Search live PubMed literature (e.g. 'CAR-T lymphoma')")
    parser.add_argument("--pan_cancer", default=None, help="Profile oncogene (e.g. EGFR, KRAS, BRAF, TP53)")
    parser.add_argument("--align", default=None, help="Pairwise sequence alignment (e.g. HBB_SICKLE, EGFR_L858R, KRAS_G12C)")
    parser.add_argument("--screen_drugs", default=None, help="Screen candidate drugs against target pocket (e.g. EGFR, KRAS, BRAF, HBB)")
    parser.add_argument("--ppi", default=None, help="Analyze PPI interactome network (e.g. EGFR, KRAS, TP53, HBB)")
    parser.add_argument("--export_json", default=None, help="Export complete pipeline results to JSON file")
    parser.add_argument("--generate_report", action="store_true", help="Generate publication-grade HTML report")
    args = parser.parse_args()
    
    run_pipeline(
        data_dir=args.data_dir,
        live_pdb=args.live_pdb,
        live_drug=args.live_drug,
        live_pubmed=args.live_pubmed,
        pan_cancer=args.pan_cancer,
        align_variant=args.align,
        screen_drugs_target=args.screen_drugs,
        ppi_gene=args.ppi,
        export_json_path=args.export_json,
        generate_report=args.generate_report
    )
