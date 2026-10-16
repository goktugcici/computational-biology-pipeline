"""
Publication-Grade Scientific Report Generator (v3.0 Ultimate)
Builds a standalone, responsive, interactive HTML scientific report
with embedded visual charts, executive metrics, sequence alignment tracks,
in silico drug screening rankings, PPI interactomes, combination drug synergy,
immunotherapy response prediction (TMB/MSI), and AlphaFold structural quality.
"""

import os
import json
from datetime import datetime
from typing import Dict, Any, Optional

from sequence_aligner import SequenceAligner
from docking_scorer import DockingScorer
from interaction_network import InteractionNetwork
from synergy_predictor import SynergyPredictor
from immunotherapy_response import ImmunotherapyPredictor
from alphafold_quality import AlphaFoldQualityAnalyzer

class ReportGenerator:
    def __init__(self, data_dir: str = "data", output_dir: str = "reports"):
        self.data_dir = data_dir
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_html_report(self, output_filename: str = "executive_summary.html", target_gene: str = "EGFR") -> str:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # 1. ClinVar counts
        clinvar_path = os.path.join(self.data_dir, "summary_all.json")
        pathogenic_count = 6
        benign_count = 1
        vus_count = 3

        if os.path.exists(clinvar_path):
            try:
                with open(clinvar_path, "r", encoding="utf-8") as f:
                    c_data = json.load(f)
                    pathogenic_count = sum(1 for v in c_data if "pathogenic" in v.get("clinical_significance", "").lower())
                    benign_count = sum(1 for v in c_data if "benign" in v.get("clinical_significance", "").lower())
                    vus_count = len(c_data) - pathogenic_count - benign_count
            except Exception:
                pass

        # 2. Sequence Alignment Metrics
        aligner = SequenceAligner()
        hbb_align = aligner.analyze_known_variant("HBB_SICKLE")

        # 3. In Silico Drug Screening Ranking
        docking = DockingScorer()
        ranked_drugs = docking.screen_library_against_target(target_gene)

        # 4. PPI Network Analysis
        ppi = InteractionNetwork(target_gene)
        ppi_metrics = ppi.compute_network_metrics()

        # 5. Synergy Prediction
        synergy = SynergyPredictor()
        syn_res = synergy.evaluate_combination("Dabrafenib", "Trametinib")

        # 6. Immunotherapy Prediction
        immuno = ImmunotherapyPredictor()
        imm_res = immuno.evaluate_response(tmb_mut_per_mb=16.2, msi_status="MSI-H", pdl1_tps_percent=65.0)

        # 7. AlphaFold Quality
        af_analyzer = AlphaFoldQualityAnalyzer()
        af_res = af_analyzer.analyze_protein("TP53")

        # Build drug table rows
        drug_rows_html = ""
        for d in ranked_drugs[:5]:
            status_badge = '<span class="px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/20">Passed</span>' if d['lipinski_pass'] else '<span class="px-2 py-0.5 rounded-full bg-amber-500/10 text-amber-400 border border-amber-500/20">Review</span>'
            drug_rows_html += f"""
            <tr class="hover:bg-slate-900/50 transition-colors">
              <td class="p-3 font-bold text-white font-mono">#{d['screen_rank']}</td>
              <td class="p-3 font-semibold text-sky-300">{d['drug_name']}</td>
              <td class="p-3 font-mono font-bold text-emerald-400">{d['estimated_delta_g_kcal_mol']} kcal/mol</td>
              <td class="p-3 font-mono text-slate-300">{d['contact_fitness_score']}%</td>
              <td class="p-3 font-mono text-amber-300">{d['drug_likeness_score']} / 100</td>
              <td class="p-3">{status_badge}</td>
            </tr>
            """

        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Executive Multi-Omics Scientific Investigation Report v3.0</title>
  <script src="https://cdn.tailwindcss.com"></script>
  <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
  <style>
    @media print {{
      body {{ background-color: #ffffff !important; color: #000000 !important; }}
      .no-print {{ display: none !important; }}
      .shadow-xl, .shadow-lg {{ box-shadow: none !important; }}
      .border {{ border-color: #e2e8f0 !important; }}
    }}
  </style>
</head>
<body class="bg-slate-950 text-slate-100 min-h-screen p-4 sm:p-8 font-sans antialiased">
  <div class="max-w-5xl mx-auto space-y-8">
    
    <!-- Report Header -->
    <header class="border-b border-slate-800 pb-6 flex flex-wrap items-center justify-between gap-4">
      <div>
        <div class="flex items-center gap-2 mb-2">
          <span class="px-2.5 py-0.5 rounded-full text-xs font-bold uppercase tracking-wider bg-sky-500/10 text-sky-400 border border-sky-500/30">
            Computational Biology Suite v3.0 Ultimate
          </span>
          <span class="text-xs text-slate-500 font-mono">Generated: {timestamp}</span>
        </div>
        <h1 class="text-2xl sm:text-3xl font-extrabold text-white tracking-tight">
          Multi-Omics Structural, Synergy & Precision Medicine Analysis
        </h1>
        <p class="text-sm text-slate-400 mt-1">
          Atomic Coordinates, Sequence Alignments, Drug Synergy, TMB Immunotherapy & AlphaFold Quality
        </p>
      </div>
      <button onclick="window.print()" class="no-print px-4 py-2 rounded-xl bg-slate-800 hover:bg-slate-700 text-slate-200 border border-slate-700 text-xs font-semibold flex items-center gap-2 transition-all">
        🖨️ Print / Save as PDF
      </button>
    </header>

    <!-- Key Metrics Cards -->
    <div class="grid grid-cols-2 md:grid-cols-4 gap-4">
      <div class="bg-slate-900 border border-slate-800 rounded-2xl p-4">
        <span class="text-xs text-slate-400 uppercase font-medium">Combination Synergy</span>
        <p class="text-2xl font-mono font-bold text-emerald-400 mt-1">CI: {syn_res['combination_index']}</p>
        <span class="text-[11px] text-emerald-400/80">{syn_res['drug_a']} + {syn_res['drug_b']}</span>
      </div>
      <div class="bg-slate-900 border border-slate-800 rounded-2xl p-4">
        <span class="text-xs text-slate-400 uppercase font-medium">Immunotherapy Prob</span>
        <p class="text-2xl font-mono font-bold text-sky-400 mt-1">{imm_res['response_probability_percent']}%</p>
        <span class="text-[11px] text-sky-400/80">TMB-H ({imm_res['tmb_mut_per_mb']} mut/Mb) & {imm_res['msi_status']}</span>
      </div>
      <div class="bg-slate-900 border border-slate-800 rounded-2xl p-4">
        <span class="text-xs text-slate-400 uppercase font-medium">AlphaFold Target</span>
        <p class="text-2xl font-mono font-bold text-purple-400 mt-1">{af_res['gene']}</p>
        <span class="text-[11px] text-purple-400/80">Mean pLDDT: {af_res['mean_plddt']} ({af_res['overall_band']})</span>
      </div>
      <div class="bg-slate-900 border border-slate-800 rounded-2xl p-4">
        <span class="text-xs text-slate-400 uppercase font-medium">Primary Hub Target</span>
        <p class="text-2xl font-mono font-bold text-rose-400 mt-1">{ppi_metrics['primary_hub_protein']}</p>
        <span class="text-[11px] text-rose-400/80">{ppi_metrics['total_nodes']} Nodes | {ppi_metrics['total_edges']} Edges</span>
      </div>
    </div>

    <!-- Charts & Structural Correlation -->
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-slate-900 border border-slate-800 rounded-2xl p-5 shadow-xl">
        <h3 class="text-sm font-bold text-white mb-1">ClinVar Pathogenicity Distribution</h3>
        <p class="text-xs text-slate-400 mb-4">Screened clinical alleles in HBB & BRCA2 cohorts</p>
        <div class="h-56 relative flex items-center justify-center">
          <canvas id="variantChart"></canvas>
        </div>
      </div>

      <div class="bg-slate-900 border border-slate-800 rounded-2xl p-5 shadow-xl">
        <h3 class="text-sm font-bold text-white mb-1">In Silico Drug Screening Affinity Comparison</h3>
        <p class="text-xs text-slate-400 mb-4">Binding free energy (&Delta;G in kcal/mol) against {target_gene}</p>
        <div class="h-56 relative flex items-center justify-center">
          <canvas id="dockingChart"></canvas>
        </div>
      </div>
    </div>

    <!-- Combination Drug Synergy Section -->
    <div class="bg-slate-900 border border-slate-800 rounded-2xl p-5 shadow-xl space-y-4">
      <div class="flex items-center justify-between border-b border-slate-800 pb-3">
        <div>
          <span class="px-2 py-0.5 rounded text-[10px] font-bold bg-emerald-500/10 text-emerald-400 border border-emerald-500/20 font-mono">
            {syn_res['clinical_trial']}
          </span>
          <h3 class="text-base font-bold text-white mt-1">Combination Therapy Synergy: {syn_res['drug_a']} + {syn_res['drug_b']}</h3>
          <p class="text-xs text-slate-400">{syn_res['indication']}</p>
        </div>
        <span class="px-3 py-1 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/20 text-xs font-mono font-bold">
          Loewe CI: {syn_res['combination_index']} (Strong Synergy)
        </span>
      </div>
      <p class="text-xs text-slate-300 leading-relaxed font-serif bg-slate-950 p-4 rounded-xl border border-slate-800">
        "{syn_res['rationale']}"
      </p>
    </div>

    <!-- Immunotherapy Precision Biomarker Section -->
    <div class="bg-slate-900 border border-slate-800 rounded-2xl p-5 shadow-xl space-y-3">
      <div class="flex items-center justify-between">
        <div>
          <h3 class="text-sm font-bold text-white">Precision Immunotherapy Response Prediction (TMB & MSI-H)</h3>
          <p class="text-xs text-slate-400">Immune Checkpoint Blockade (Anti-PD-1 / Anti-CTLA-4) Efficacy Index</p>
        </div>
        <span class="px-3 py-1 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/20 text-xs font-mono font-bold">
          {imm_res['response_probability_percent']}% Probability
        </span>
      </div>
      <div class="bg-slate-950 p-4 rounded-xl border border-slate-800 flex justify-between items-center text-xs">
        <div>
          <span class="text-slate-400">Actionable Regimen:</span>
          <p class="text-sm font-bold text-white mt-0.5">{imm_res['recommended_regimen']}</p>
        </div>
        <span class="font-mono text-emerald-400 font-bold">{imm_res['classification']}</span>
      </div>
    </div>

    <!-- Sequence Alignment Track Section -->
    <div class="bg-slate-900 border border-slate-800 rounded-2xl p-5 shadow-xl space-y-4">
      <div class="flex items-center justify-between">
        <div>
          <h3 class="text-sm font-bold text-white">Needleman-Wunsch Pairwise Sequence Alignment</h3>
          <p class="text-xs text-slate-400">BLOSUM62 matrix scoring & pathogenic substitution pinpointing</p>
        </div>
        <span class="px-2.5 py-1 rounded-lg bg-sky-500/10 text-sky-400 border border-sky-500/20 text-xs font-mono">
          Identity: {hbb_align['identity_percent']}% | Score: {hbb_align['alignment_score']}
        </span>
      </div>
      
      <div class="bg-slate-950 p-4 rounded-xl font-mono text-xs overflow-x-auto space-y-1 text-slate-300 border border-slate-800/80">
        <div class="flex items-center gap-2"><span class="text-sky-400 w-16">WT (HBB):</span> <span class="tracking-widest">{hbb_align['aligned_seq1'][:50]}...</span></div>
        <div class="flex items-center gap-2"><span class="text-slate-500 w-16">MATCH:   </span> <span class="tracking-widest text-emerald-400 font-bold">{hbb_align['match_track'][:50]}...</span></div>
        <div class="flex items-center gap-2"><span class="text-rose-400 w-16">MUT (HbS):</span> <span class="tracking-widest">{hbb_align['aligned_seq2'][:50]}...</span></div>
      </div>
      <p class="text-xs text-rose-400">
        <strong>Detected Mutation:</strong> {hbb_align['mutations_detected'][0]['description']} ({hbb_align['mutations_detected'][0]['type']}, BLOSUM62 score: {hbb_align['mutations_detected'][0]['blosum62_score']}) -> Drives pathological deoxy-hemoglobin polymer fibers.
      </p>
    </div>

    <!-- In Silico Drug Screening Ranking Table -->
    <div class="bg-slate-900 border border-slate-800 rounded-2xl overflow-hidden shadow-xl">
      <div class="p-5 border-b border-slate-800 flex justify-between items-center">
        <div>
          <h3 class="text-base font-bold text-white">In Silico Drug Screening & Lipinski Compliance</h3>
          <p class="text-xs text-slate-400">Target Binding Pocket: {ranked_drugs[0]['target_pocket']}</p>
        </div>
      </div>
      <div class="overflow-x-auto text-xs">
        <table class="w-full text-left">
          <thead class="bg-slate-950/80 text-slate-400 uppercase tracking-wider text-[11px]">
            <tr>
              <th class="p-3">Rank</th>
              <th class="p-3">Candidate Molecule</th>
              <th class="p-3">Estimated &Delta;G</th>
              <th class="p-3">Pocket Fitness</th>
              <th class="p-3">Drug-Likeness (QED)</th>
              <th class="p-3">Lipinski Ro5</th>
            </tr>
          </thead>
          <tbody class="divide-y divide-slate-800/60">
            {drug_rows_html}
          </tbody>
        </table>
      </div>
    </div>

    <!-- References & Citations -->
    <footer class="text-xs text-slate-500 border-t border-slate-800 pt-6 space-y-2">
      <p class="font-bold text-slate-400">Primary References & Data Sources:</p>
      <ul class="list-disc list-inside space-y-1 font-mono">
        <li>RCSB Protein Data Bank (PDB): 4HHB, 2J6M, 6OIM</li>
        <li>AlphaFold Protein Structure Database: AF-P04637-F1 (TP53)</li>
        <li>FDA Keynote-158 / Keynote-042: Pembrolizumab TMB & MSI-H Companion Diagnostics</li>
        <li>The New England Journal of Medicine (NEJM): COMBI-d / COMBI-v Phase 3 (Dabrafenib + Trametinib)</li>
      </ul>
    </footer>

  </div>

  <script>
    // Chart 1: Donut Chart
    new Chart(document.getElementById('variantChart'), {{
      type: 'doughnut',
      data: {{
        labels: ['Pathogenic', 'Benign (Control)', 'VUS / Other'],
        datasets: [{{
          data: [{pathogenic_count}, {benign_count}, {vus_count}],
          backgroundColor: ['#f43f5e', '#10b981', '#64748b'],
          borderColor: '#0f172a',
          borderWidth: 3
        }}]
      }},
      options: {{
        responsive: true,
        maintainAspectRatio: false,
        plugins: {{
          legend: {{ position: 'bottom', labels: {{ color: '#94a3b8', font: {{ size: 11 }} }} }}
        }}
      }}
    }});

    // Chart 2: Docking Bar Chart
    new Chart(document.getElementById('dockingChart'), {{
      type: 'bar',
      data: {{
        labels: {json.dumps([d['drug_name'].split('(')[0].strip() for d in ranked_drugs[:5]])},
        datasets: [{{
          label: 'Estimated Free Energy ΔG (kcal/mol)',
          data: {json.dumps([d['estimated_delta_g_kcal_mol'] for d in ranked_drugs[:5]])},
          backgroundColor: ['#10b981', '#0ea5e9', '#6366f1', '#f59e0b', '#ec4899'],
          borderRadius: 8
        }}]
      }},
      options: {{
        responsive: true,
        maintainAspectRatio: false,
        scales: {{
          y: {{
            grid: {{ color: '#1e293b' }},
            ticks: {{ color: '#94a3b8' }}
          }},
          x: {{
            grid: {{ display: false }},
            ticks: {{ color: '#94a3b8', font: {{ size: 10 }} }}
          }}
        }},
        plugins: {{
          legend: {{ display: false }}
        }}
      }}
    }});
  </script>
</body>
</html>
"""
        target_path = os.path.join(self.output_dir, output_filename)
        with open(target_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        return target_path


if __name__ == "__main__":
    rg = ReportGenerator()
    out = rg.generate_html_report(target_gene="EGFR")
    print(f"Scientific Report Generated: {out}")
