"""
Protein-Protein Interaction (PPI) Network and Interactome Analysis Engine.
Models protein interaction topologies, calculates degree centrality and hub bottlenecks,
and formats network graphs for interactive web visualization (STRING DB schema).
"""

import json
from pathlib import Path
from typing import Dict, List, Any, Optional

CURATED_PPI_NETWORKS: Dict[str, Dict[str, Any]] = {
    'EGFR': {
        'gene': 'EGFR',
        'name': 'Epidermal Growth Factor Receptor Interactome',
        'nodes': [
            {'id': 'EGFR', 'role': 'Hub Target', 'type': 'Receptor Kinase', 'score': 1.0},
            {'id': 'ERBB2', 'role': 'Heterodimerization Partner', 'type': 'Receptor Kinase', 'score': 0.99},
            {'id': 'GRB2', 'role': 'Adaptor Protein', 'type': 'SH2/SH3 Adaptor', 'score': 0.98},
            {'id': 'SOS1', 'role': 'GEF Activator', 'type': 'Exchange Factor', 'score': 0.95},
            {'id': 'SHC1', 'role': 'Phosphotyrosine Adaptor', 'type': 'Adaptor', 'score': 0.96},
            {'id': 'PIK3CA', 'role': 'Lipid Kinase', 'type': 'Survival Signaling', 'score': 0.94},
            {'id': 'STAT3', 'role': 'Transcription Factor', 'type': 'Transcriptional Regulation', 'score': 0.91},
            {'id': 'CBL', 'role': 'E3 Ubiquitin Ligase', 'type': 'Receptor Downregulation', 'score': 0.93}
        ],
        'edges': [
            {'source': 'EGFR', 'target': 'ERBB2', 'combined_score': 0.999, 'evidence': 'Biochemical Assay'},
            {'source': 'EGFR', 'target': 'GRB2', 'combined_score': 0.998, 'evidence': 'Phosphorylation Contact'},
            {'source': 'GRB2', 'target': 'SOS1', 'combined_score': 0.997, 'evidence': 'Complex Assembly'},
            {'source': 'EGFR', 'target': 'SHC1', 'combined_score': 0.995, 'evidence': 'Direct Phosphorylation'},
            {'source': 'SHC1', 'target': 'GRB2', 'combined_score': 0.991, 'evidence': 'Signal Transduction'},
            {'source': 'EGFR', 'target': 'PIK3CA', 'combined_score': 0.978, 'evidence': 'Downstream Cascade'},
            {'source': 'EGFR', 'target': 'STAT3', 'combined_score': 0.965, 'evidence': 'Transcriptional Activation'},
            {'source': 'EGFR', 'target': 'CBL', 'combined_score': 0.972, 'evidence': 'Ubiquitination'}
        ]
    },
    'KRAS': {
        'gene': 'KRAS',
        'name': 'KRAS GTPase Signalosome',
        'nodes': [
            {'id': 'KRAS', 'role': 'Hub Target', 'type': 'Small GTPase', 'score': 1.0},
            {'id': 'RAF1', 'role': 'Primary Kinase Effector', 'type': 'Ser/Thr Kinase', 'score': 0.99},
            {'id': 'BRAF', 'role': 'MAPK Cascade Kinase', 'type': 'Kinase', 'score': 0.98},
            {'id': 'PIK3CA', 'role': 'PI3K Pathway Effector', 'type': 'Lipid Kinase', 'score': 0.97},
            {'id': 'RALGDS', 'role': 'Ral GEF Effector', 'type': 'GEF', 'score': 0.92},
            {'id': 'SOS1', 'role': 'GDP/GTP Exchange Factor', 'type': 'Activator', 'score': 0.96},
            {'id': 'NF1', 'role': 'Neurofibromin GTPase Activating Protein', 'type': 'Tumor Suppressor', 'score': 0.95},
            {'id': 'MAP2K1', 'role': 'MEK1 Kinase', 'type': 'Downstream Kinase', 'score': 0.94}
        ],
        'edges': [
            {'source': 'KRAS', 'target': 'RAF1', 'combined_score': 0.999, 'evidence': 'Effector Binding'},
            {'source': 'KRAS', 'target': 'BRAF', 'combined_score': 0.995, 'evidence': 'Effector Binding'},
            {'source': 'KRAS', 'target': 'PIK3CA', 'combined_score': 0.991, 'evidence': 'Lipid Signaling'},
            {'source': 'KRAS', 'target': 'RALGDS', 'combined_score': 0.965, 'evidence': 'Secondary Effector'},
            {'source': 'SOS1', 'target': 'KRAS', 'combined_score': 0.998, 'evidence': 'Nucleotide Exchange'},
            {'source': 'NF1', 'target': 'KRAS', 'combined_score': 0.994, 'evidence': 'GTP Hydrolysis'},
            {'source': 'RAF1', 'target': 'MAP2K1', 'combined_score': 0.997, 'evidence': 'Phosphorylation Cascade'},
            {'source': 'BRAF', 'target': 'MAP2K1', 'combined_score': 0.998, 'evidence': 'Phosphorylation Cascade'}
        ]
    },
    'TP53': {
        'gene': 'TP53',
        'name': 'p53 Guardian of the Genome Network',
        'nodes': [
            {'id': 'TP53', 'role': 'Master Hub', 'type': 'Tumor Suppressor TF', 'score': 1.0},
            {'id': 'MDM2', 'role': 'Negative Feedback Regulator', 'type': 'E3 Ligase', 'score': 0.99},
            {'id': 'CDKN1A', 'role': 'p21 Cell Cycle Arrest', 'type': 'CDK Inhibitor', 'score': 0.99},
            {'id': 'ATM', 'role': 'DNA Damage Sensor', 'type': 'PI3K-related Kinase', 'score': 0.98},
            {'id': 'CHEK2', 'role': 'Checkpoint Kinase', 'type': 'Ser/Thr Kinase', 'score': 0.96},
            {'id': 'BAX', 'role': 'Pro-Apoptotic Factor', 'type': 'Bcl-2 Family', 'score': 0.97},
            {'id': 'EP300', 'role': 'Histone Acetyltransferase', 'type': 'Co-activator', 'score': 0.95}
        ],
        'edges': [
            {'source': 'TP53', 'target': 'MDM2', 'combined_score': 0.999, 'evidence': 'Auto-regulatory Feedback'},
            {'source': 'MDM2', 'target': 'TP53', 'combined_score': 0.999, 'evidence': 'Proteasomal Degradation'},
            {'source': 'TP53', 'target': 'CDKN1A', 'combined_score': 0.998, 'evidence': 'Direct Transcription'},
            {'source': 'ATM', 'target': 'TP53', 'combined_score': 0.995, 'evidence': 'Ser15 Phosphorylation'},
            {'source': 'ATM', 'target': 'CHEK2', 'combined_score': 0.994, 'evidence': 'Activation Phosphorylation'},
            {'source': 'CHEK2', 'target': 'TP53', 'combined_score': 0.989, 'evidence': 'Ser20 Phosphorylation'},
            {'source': 'TP53', 'target': 'BAX', 'combined_score': 0.992, 'evidence': 'Apoptotic Induction'},
            {'source': 'EP300', 'target': 'TP53', 'combined_score': 0.985, 'evidence': 'Lys382 Acetylation'}
        ]
    },
    'HBB': {
        'gene': 'HBB',
        'name': 'Hemoglobin Subunit Beta Erythroid Network',
        'nodes': [
            {'id': 'HBB', 'role': 'Target Subunit', 'type': 'Beta Globin', 'score': 1.0},
            {'id': 'HBA1', 'role': 'Tetramer Partner', 'type': 'Alpha Globin 1', 'score': 0.99},
            {'id': 'HBA2', 'role': 'Tetramer Partner', 'type': 'Alpha Globin 2', 'score': 0.99},
            {'id': 'AHSP', 'role': 'Alpha-Hemoglobin Stabilizing Protein', 'type': 'Chaperone', 'score': 0.97},
            {'id': 'BPGM', 'role': '2,3-DPG Regulator', 'type': 'Mutase/Phosphatase', 'score': 0.92},
            {'id': 'ALAS2', 'role': 'Heme Biosynthesis', 'type': 'Enzyme', 'score': 0.90}
        ],
        'edges': [
            {'source': 'HBB', 'target': 'HBA1', 'combined_score': 0.999, 'evidence': 'Tetramerization'},
            {'source': 'HBB', 'target': 'HBA2', 'combined_score': 0.999, 'evidence': 'Tetramerization'},
            {'source': 'HBA1', 'target': 'AHSP', 'combined_score': 0.996, 'evidence': 'Chaperone Stabilization'},
            {'source': 'HBA2', 'target': 'AHSP', 'combined_score': 0.996, 'evidence': 'Chaperone Stabilization'},
            {'source': 'BPGM', 'target': 'HBB', 'combined_score': 0.965, 'evidence': 'Allosteric Effector Modulation'},
            {'source': 'ALAS2', 'target': 'HBB', 'combined_score': 0.942, 'evidence': 'Heme Incorporation'}
        ]
    }
}


class InteractionNetwork:
    """
    Analyzes protein-protein interaction networks and calculates graph metrics
    such as degree centrality and hub bottlenecks.
    """

    def __init__(self, gene: str = 'EGFR'):
        self.gene = gene.upper()

    def get_network(self) -> Dict[str, Any]:
        """Retrieves PPI network data for the designated gene."""
        gene_key = self.gene if self.gene in CURATED_PPI_NETWORKS else 'EGFR'
        return CURATED_PPI_NETWORKS[gene_key]

    def compute_network_metrics(self) -> Dict[str, Any]:
        """
        Calculates node degrees, degree centrality, density,
        and identifies hub bottleneck proteins.
        """
        net = self.get_network()
        nodes = net['nodes']
        edges = net['edges']

        degree_map: Dict[str, int] = {n['id']: 0 for n in nodes}
        edge_weight_sum: Dict[str, float] = {n['id']: 0.0 for n in nodes}

        for edge in edges:
            src = edge['source']
            tgt = edge['target']
            weight = edge.get('combined_score', 0.5)
            if src in degree_map:
                degree_map[src] += 1
                edge_weight_sum[src] += weight
            if tgt in degree_map:
                degree_map[tgt] += 1
                edge_weight_sum[tgt] += weight

        total_nodes = len(nodes)
        max_possible_edges = (total_nodes * (total_nodes - 1)) / 2 if total_nodes > 1 else 1
        density = round(len(edges) / max_possible_edges, 3)

        # Degree centrality: degree / (N - 1)
        centrality = {}
        for nid, deg in degree_map.items():
            cent = round(deg / (total_nodes - 1), 3) if total_nodes > 1 else 1.0
            centrality[nid] = {
                'degree': deg,
                'centrality': cent,
                'weighted_connectivity': round(edge_weight_sum[nid], 3)
            }

        # Identify top hub nodes (highest degree)
        sorted_hubs = sorted(centrality.items(), key=lambda item: item[1]['degree'], reverse=True)
        primary_hub = sorted_hubs[0][0] if sorted_hubs else self.gene

        return {
            'gene': self.gene,
            'network_name': net['name'],
            'total_nodes': total_nodes,
            'total_edges': len(edges),
            'network_density': density,
            'primary_hub_protein': primary_hub,
            'node_centrality': centrality,
            'top_ranked_hubs': [item[0] for item in sorted_hubs[:3]],
            'raw_nodes': nodes,
            'raw_edges': edges
        }


if __name__ == '__main__':
    ppi = InteractionNetwork('EGFR')
    metrics = ppi.compute_network_metrics()
    print(f"Network: {metrics['network_name']}")
    print(f"Nodes: {metrics['total_nodes']} | Edges: {metrics['total_edges']} | Density: {metrics['network_density']}")
    print(f"Primary Hub: {metrics['primary_hub_protein']} | Top Hubs: {metrics['top_ranked_hubs']}")
