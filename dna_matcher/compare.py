from typing import Dict
from .algorithms import calculate_similarity_and_alignment
from .db import get_comparison_result, save_comparison_result

def compare_species_matrix(species_sequences: Dict[str, str], marker: str, max_len: int = 1000):
    """Compare all pairs from species_sequences and return a similarity matrix."""
    species = list(species_sequences.keys())
    matrix = {name: {} for name in species}
    details = {}
    
    best_global = {"pair": None, "score": -1.0}
    best_local = {"pair": None, "score": -1.0}

    for i, a in enumerate(species):
        for b in species[i:]:
            if a == b:
                # Trims the sequence to max_len to ensure self-comparisons match cross-comparisons
                seq = species_sequences[a]
                compare_len = min(len(seq), max_len)
                trimmed_seq = seq[:compare_len]
                
                result = {
                    "hirschberg_pct": 100.0,
                    "sw_pct": 100.0,
                    "sw_score": 2 * compare_len, # Match score is 2 per base
                    "hirschberg_match": trimmed_seq,
                    "sw_match": trimmed_seq,
                    "compare_len": compare_len,
                }
            else:
                cached = get_comparison_result(a, b, marker, max_len)
                if cached is not None:
                    result = cached
                else:
                    result = calculate_similarity_and_alignment(
                        species_sequences[a], species_sequences[b], max_len=max_len
                    )
                    save_comparison_result(a, b, marker, max_len, result)

            score = result["hirschberg_pct"]
            matrix[a][b] = score
            matrix[b][a] = score
            details[(a, b)] = result
            details[(b, a)] = result

            if a != b:
                if score > best_global["score"]:
                    best_global = {"pair": (a, b), "score": score}
                if result["sw_pct"] > best_local["score"]:
                    best_local = {"pair": (a, b), "score": result["sw_pct"]}

    for sp in species:
        matrix[sp][sp] = 100.0

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
        a: {b: 100.0 - matrix[a][b] for b in matrix[a]}
        for a in matrix
    }
    
    average_distance = {
        sp: sum(distance[sp].values()) / max(1, len(distance[sp]) - 1)
        for sp in species
    }
    
    sorted_species = sorted(average_distance, key=average_distance.get)

    lines = ["Neighbor-Joining Approximation:"]
    for sp in sorted_species:
        lines.append(f"- {sp} (avg distance: {average_distance[sp]:.2f})")

    return "\n".join(lines)