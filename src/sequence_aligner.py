"""
Pairwise Sequence Alignment and Mutation Consequence Predictor.
Implements Needleman-Wunsch (Global) alignment, BLOSUM62 scoring,
and amino acid mutation detection for oncogenic and genetic variants.
"""

from typing import Dict, List, Tuple, Any, Optional

# Standard BLOSUM62 substitution matrix subset for 20 amino acids
BLOSUM62_MATRIX: Dict[str, Dict[str, int]] = {
    'A': {'A':  4, 'R': -1, 'N': -2, 'D': -2, 'C':  0, 'Q': -1, 'E': -1, 'G':  0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S':  1, 'T':  0, 'W': -3, 'Y': -2, 'V':  0},
    'R': {'A': -1, 'R':  5, 'N':  0, 'D': -2, 'C': -3, 'Q':  1, 'E':  0, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R':  0, 'N':  6, 'D':  1, 'C': -3, 'Q':  0, 'E':  0, 'G':  0, 'H':  1, 'I': -3, 'L': -3, 'K':  0, 'M': -2, 'F': -3, 'P': -2, 'S':  1, 'T':  0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N':  1, 'D':  6, 'C': -3, 'Q':  0, 'E':  2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C':  9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R':  1, 'N':  0, 'D':  0, 'C': -3, 'Q':  5, 'E':  2, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  1, 'M':  0, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R':  0, 'N':  0, 'D':  2, 'C': -4, 'Q':  2, 'E':  5, 'G': -2, 'H':  0, 'I': -3, 'L': -3, 'K':  1, 'M': -2, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A':  0, 'R': -2, 'N':  0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G':  6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S':  0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R':  0, 'N':  1, 'D': -1, 'C': -3, 'Q':  0, 'E':  0, 'G': -2, 'H':  8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y':  2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I':  4, 'L':  2, 'K': -3, 'M':  1, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V':  3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I':  2, 'L':  4, 'K': -2, 'M':  2, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V':  1},
    'K': {'A': -1, 'R':  2, 'N':  0, 'D': -1, 'C': -3, 'Q':  1, 'E':  1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K':  5, 'M': -1, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q':  0, 'E': -2, 'G': -3, 'H': -2, 'I':  1, 'L':  2, 'K': -1, 'M':  5, 'F':  0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V':  1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I':  0, 'L':  0, 'K': -3, 'M':  0, 'F':  6, 'P': -4, 'S': -2, 'T': -2, 'W':  1, 'Y':  3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P':  7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A':  1, 'R': -1, 'N':  1, 'D':  0, 'C': -1, 'Q':  0, 'E':  0, 'G':  0, 'H': -1, 'I': -2, 'L': -2, 'K':  0, 'M': -1, 'F': -2, 'P': -1, 'S':  4, 'T':  1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A':  0, 'R': -1, 'N':  0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S':  1, 'T':  5, 'W': -2, 'Y': -2, 'V':  0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F':  1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y':  2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H':  2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F':  3, 'P': -3, 'S': -2, 'T': -2, 'W':  2, 'Y':  7, 'V': -1},
    'V': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I':  3, 'L':  1, 'K': -2, 'M':  1, 'F': -1, 'P': -2, 'S': -2, 'T':  0, 'W': -3, 'Y': -1, 'V':  4}
}

# Standard reference sequence catalog (wild-type & oncogenic mutants)
KNOWN_VARIANTS: Dict[str, Dict[str, Any]] = {
    'HBB_SICKLE': {
        'gene': 'HBB',
        'condition': 'Sickle Cell Anemia',
        'mutation': 'p.Glu7Val (HbS)',
        'wild_type': 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH',
        'mutant':    'VHLTPVEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
    },
    'EGFR_L858R': {
        'gene': 'EGFR',
        'condition': 'Non-Small Cell Lung Cancer (TKI Sensitivity)',
        'mutation': 'p.L858R (Exon 21)',
        'wild_type': 'KITDFGRAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQG',
        'mutant':    'KITDFGRAKLRGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQG'
    },
    'KRAS_G12C': {
        'gene': 'KRAS',
        'condition': 'Pancreatic / NSCLC Carcinoma',
        'mutation': 'p.G12C (Exon 2 Codon 12)',
        'wild_type': 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK',
        'mutant':    'MTEYKLVVVGACGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'
    },
    'BRAF_V600E': {
        'gene': 'BRAF',
        'condition': 'Metastatic Melanoma / Colorectal Cancer',
        'mutation': 'p.V600E (Kinase Domain Activation)',
        'wild_type': 'LDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYG',
        'mutant':    'LDFGLATVEISRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYG'
    }
}


class SequenceAligner:
    """
    Performs dynamic programming Needleman-Wunsch global alignment
    between wild-type and mutated biological sequences.
    """

    def __init__(self, gap_penalty: int = -4):
        self.gap_penalty = gap_penalty

    def _score_pair(self, a: str, b: str) -> int:
        """Returns substitution score from BLOSUM62 or default identity/mismatch."""
        a_up = a.upper()
        b_up = b.upper()
        if a_up in BLOSUM62_MATRIX and b_up in BLOSUM62_MATRIX[a_up]:
            return BLOSUM62_MATRIX[a_up][b_up]
        return 4 if a_up == b_up else -2

    def align(self, seq1: str, seq2: str) -> Dict[str, Any]:
        """
        Computes optimal global alignment using Needleman-Wunsch algorithm.
        Returns alignment strings, score, percent identity, similarity, and detected mutations.
        """
        n, m = len(seq1), len(seq2)
        score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

        for i in range(n + 1):
            score_matrix[i][0] = i * self.gap_penalty
        for j in range(m + 1):
            score_matrix[0][j] = j * self.gap_penalty

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                match = score_matrix[i - 1][j - 1] + self._score_pair(seq1[i - 1], seq2[j - 1])
                delete = score_matrix[i - 1][j] + self.gap_penalty
                insert = score_matrix[i][j - 1] + self.gap_penalty
                score_matrix[i][j] = max(match, delete, insert)

        # Traceback
        align1, align2, match_track = [], [], []
        i, j = n, m
        mutations = []

        while i > 0 or j > 0:
            current_score = score_matrix[i][j]
            if i > 0 and j > 0 and current_score == score_matrix[i - 1][j - 1] + self._score_pair(seq1[i - 1], seq2[j - 1]):
                char1 = seq1[i - 1]
                char2 = seq2[j - 1]
                align1.append(char1)
                align2.append(char2)
                if char1 == char2:
                    match_track.append('|')
                else:
                    pair_score = self._score_pair(char1, char2)
                    match_track.append(':' if pair_score > 0 else '.')
                    mutations.append({
                        'position': i,
                        'ref': char1,
                        'alt': char2,
                        'type': 'Substitution',
                        'blosum62_score': pair_score,
                        'description': f"p.{char1}{i}{char2}"
                    })
                i -= 1
                j -= 1
            elif i > 0 and current_score == score_matrix[i - 1][j] + self.gap_penalty:
                char1 = seq1[i - 1]
                align1.append(char1)
                align2.append('-')
                match_track.append(' ')
                mutations.append({
                    'position': i,
                    'ref': char1,
                    'alt': '-',
                    'type': 'Deletion',
                    'blosum62_score': self.gap_penalty,
                    'description': f"p.{char1}{i}del"
                })
                i -= 1
            else:
                char2 = seq2[j - 1]
                align1.append('-')
                align2.append(char2)
                match_track.append(' ')
                mutations.append({
                    'position': i + 1,
                    'ref': '-',
                    'alt': char2,
                    'type': 'Insertion',
                    'blosum62_score': self.gap_penalty,
                    'description': f"p.ins{i+1}{char2}"
                })
                j -= 1

        aligned_seq1 = "".join(reversed(align1))
        aligned_seq2 = "".join(reversed(align2))
        alignment_match = "".join(reversed(match_track))
        mutations.reverse()

        total_positions = len(aligned_seq1)
        identical_count = alignment_match.count('|')
        similar_count = identical_count + alignment_match.count(':')
        identity_pct = round((identical_count / total_positions) * 100, 2) if total_positions > 0 else 0.0
        similarity_pct = round((similar_count / total_positions) * 100, 2) if total_positions > 0 else 0.0

        return {
            'aligned_seq1': aligned_seq1,
            'aligned_seq2': aligned_seq2,
            'match_track': alignment_match,
            'alignment_score': score_matrix[n][m],
            'total_positions': total_positions,
            'identical_residues': identical_count,
            'identity_percent': identity_pct,
            'similarity_percent': similarity_pct,
            'mutations_detected': mutations
        }

    def analyze_known_variant(self, variant_key: str) -> Dict[str, Any]:
        """Runs alignment on a known canonical variant key (e.g. 'HBB_SICKLE', 'EGFR_L858R')."""
        if variant_key not in KNOWN_VARIANTS:
            raise KeyError(f"Variant '{variant_key}' not in catalog. Choose from: {list(KNOWN_VARIANTS.keys())}")
        
        entry = KNOWN_VARIANTS[variant_key]
        result = self.align(entry['wild_type'], entry['mutant'])
        result['metadata'] = {
            'variant_key': variant_key,
            'gene': entry['gene'],
            'condition': entry['condition'],
            'expected_mutation': entry['mutation']
        }
        return result


if __name__ == '__main__':
    aligner = SequenceAligner()
    res = aligner.analyze_known_variant('HBB_SICKLE')
    print(f"Alignment Score: {res['alignment_score']}")
    print(f"Identity: {res['identity_percent']}% | Similarity: {res['similarity_percent']}%")
    print(f"Mutations Detected: {res['mutations_detected']}")
