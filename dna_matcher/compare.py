from typing import List, Tuple, Dict
from .algorithms import calculate_similarity_and_alignment

def compare_species_matrix(species_sequences: Dict[str, str], marker: str, max_len: int = 1000):
    """Compare all pairs from species_sequences and return a similarity matrix."""
    species = list(species_sequences.keys())
    matrix = {name: {} for name in species}
    details = {}
    
    # Initialize trackers for the best alignments required by the UI
    best_global = {"pair": ("", ""), "score": -1.0}
    best_local = {"pair": ("", ""), "score": -1.0}

    for i, a in enumerate(species):
        for b in species[i:]:
            seq_a = species_sequences[a]
            seq_b = species_sequences[b]
            
            result = calculate_similarity_and_alignment(seq_a, seq_b, max_len=max_len)
            score = result["hirschberg_pct"]
            
            matrix[a][b] = score
            matrix[b][a] = score
            details[(a, b)] = result
            
            if a != b:
                details[(b, a)] = result
                
                # Track best matches
                if result["hirschberg_pct"] > best_global["score"]:
                    best_global = {"pair": (a, b), "score": result["hirschberg_pct"]}
                if result["sw_pct"] > best_local["score"]:
                    best_local = {"pair": (a, b), "score": result["sw_pct"]}

    return {
        "species": species,
        "matrix": matrix,
        "details": details,
        "best_global": best_global,
        "best_local": best_local
    }

def build_neighbor_joining_tree(matrix_data):
    """Builds a simple neighbor-joining tree from a similarity matrix."""
    species = matrix_data["species"]
    matrix = matrix_data["matrix"]
    distance = {
        a: {b: 100.0 - matrix[a][b] for b in matrix[a]} for a in matrix
    }

    weights = {}
    for sp in species:
        weights[sp] = sum(distance[sp].values()) / max(1, len(distance[sp]) - 1)

    sorted_species = sorted(weights, key=weights.get)
    tree_lines = ["(Neighbor-Joining Approximation):"]
    for sp in sorted_species:
        tree_lines.append(f"- {sp} (avg distance: {weights[sp]:.2f})")

    return "\n".join(tree_lines)