# dna_matcher/core.py
import os
from collections import OrderedDict
from Bio import Entrez
from dotenv import load_dotenv

from .fetch import fetch_marker_sequence
from .compare import compare_species_matrix, build_neighbor_joining_tree
from .db import (
    init_db,
    get_marker_sequence as db_get_marker_sequence,
    save_marker_sequence,
    get_search_results,
    save_search_results,
)

load_dotenv()
Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")

init_db()

MARKER_CHAINS = {
    "animal":  ["COI", "16S", "12S", "cytb"],
    "plant":   ["rbcL", "matK", "ITS"],
    "fungus":  ["ITS", "LSU", "SSU"],
    "default": ["COI", "16S", "ITS", "rbcL", "12S"],
}

MARKER_LENGTH_RANGE = {
    # Keep sensible minimums but allow a larger upper bound so long comparisons
    # (up to 10,000 bp) are possible when the frontend slider requests them.
    "COI":  (500,  10000),
    "16S":  (400,  10000),
    "12S":  (300,  10000),
    "cytb": (800,  10000),
    "rbcL": (400,  10000),
    "matK": (800,  10000),
    "ITS":  (400,  10000),
    "LSU":  (800,  10000),
    "SSU":  (1600, 10000),
}

NCBI_TIMEOUT = float(os.environ.get("NCBI_TIMEOUT", 20.0))

def search_species_options(common_name: str, max_results: int = 6):
    def _fetch():
        # Prefer explicit common-name and all-names field searches before falling
        # back to plain or quoted text. This avoids generic term searches that can
        # return unrelated taxa when the input matches other text fields.
        query_variants = [
            f'"{common_name}"[Common Name]',
            f'"{common_name}"[All Names]',
            f'{common_name}[Common Name]',
            f'{common_name}[All Names]',
            f'"{common_name}"[Text Word]',
            f'{common_name}[Text Word]',
            f'"{common_name}"',
            common_name,
        ]

        records = None
        for q in query_variants:
            try:
                handle = Entrez.esearch(db="taxonomy", term=q, retmax=max_results)
                record = Entrez.read(handle)
                handle.close()
            except Exception:
                record = {"IdList": []}

            if record and record.get("IdList"):
                try:
                    handle = Entrez.efetch(db="taxonomy", id=record["IdList"], retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()
                except Exception:
                    records = None

            if records:
                break

        if not records:
            return []

        results = []
        for r in records:
            sci = r.get("ScientificName", "")
            common = r.get("GenbankCommonName") or r.get("CommonName") or sci
            results.append((common, sci))
        return results


def get_scientific_name(common_name: str):
    options = search_species_options(common_name, max_results=1)
    return options[0][1] if options else None

def get_marker_sequence(species_name: str, gene: str):
    existing = db_get_marker_sequence(species_name, gene)
    if existing:
        return existing

    seq = fetch_marker_sequence(species_name, gene)
    if seq:
        save_marker_sequence(species_name, gene, seq)
    return seq(species_name, gene, seq)

def _build_marker_order(types):
    seen = set()
    order = []
    for t in types:
        chain = MARKER_CHAINS.get(t, MARKER_CHAINS["default"])
        for gene in chain:
            if gene not in seen:
                seen.add(gene)
                order.append(gene)
    for gene in MARKER_CHAINS["default"]:
        if gene not in seen:
            seen.add(gene)
            order.append(gene)
    return order


def fetch_common_marker_sequences(species_names: tuple, types: tuple):
    markers = _build_marker_order(types)
    for gene in markers:
        sequences = OrderedDict()
        for species_name in species_names:
            seq = get_marker_sequence(species_name, gene)
            if not seq:
                print(f"[NCBI] {species_name} missing {gene}")
                break
            sequences[species_name] = seq

        if len(sequences) == len(species_names):
            print(f"[NCBI] All species matched on marker: {gene}")
            return {"marker": gene, "sequences": sequences}

    return None


def compare_species_matrix(species_sequences: dict, marker: str, max_len: int = 1000):
    species = list(species_sequences.keys())
    matrix = {name: {} for name in species}
    details = {}
    best_global = {"pair": None, "score": -1.0}
    best_local = {"pair": None, "score": -1.0}

    for i, a in enumerate(species):
        for b in species[i:]:
            if a == b:
                result = {
                    "hirschberg_pct": 100.0,
                    "sw_pct": 100.0,
                    "sw_score": 0,
                    "hirschberg_match": species_sequences[a],
                    "sw_match": species_sequences[a],
                    "compare_len": len(species_sequences[a]),
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
        "best_local": best_local,
    }


def build_neighbor_joining_tree(matrix_data):
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
