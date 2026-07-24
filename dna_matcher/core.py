# dna_matcher/core.py — package-level orchestration for multi-species DNA matching

import os
import time
from collections import OrderedDict

import streamlit as st
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

from .algorithms import calculate_similarity_and_alignment
from .db import (
    init_db,
    get_marker_sequence as db_get_marker_sequence,
    save_marker_sequence,
    get_search_results,
    save_search_results,
    get_comparison_result,
    save_comparison_result,
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


def _fetch_with_retry(fn, retries=3, base_delay=1.0):
    for attempt in range(retries):
        try:
            result = fn()
            if result is not None:
                return result
        except Exception as e:
            print(f"[Retry {attempt + 1}/{retries}] {e}")
        time.sleep(base_delay * (2 ** attempt))
    return None


def search_species_options(common_name: str, max_results: int = 6):
    def _fetch():
        # Try a set of query variants to improve common-name resolution (e.g. "tiger",
        # "royal bengal tiger", or field-restricted queries). Stop at the first
        # non-empty result set.
        query_variants = [
            common_name,
            f'"{common_name}"',
            f"{common_name}[All Names]",
            f"{common_name}[Common Name]",
            f"{common_name}[Text Word]",
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
            # record last attempted queries and any final IdList for UI debugging
            try:
                st.session_state.setdefault("tax_debug", {})
                st.session_state["tax_debug"][common_name] = {
                    "queries_tried": query_variants,
                    "last_idlist": record.get("IdList") if record else [],
                }
            except Exception:
                # In contexts without Streamlit session state, ignore silently
                pass
            return []

        results = []
        for r in records:
            sci = r.get("ScientificName", "")
            common = r.get("GenbankCommonName") or r.get("CommonName") or sci
            results.append((common, sci))
        return results

    cached = get_search_results(common_name)
    if cached:
        return cached[:max_results]

    results = _fetch_with_retry(_fetch) or []
    if results:
        save_search_results(common_name, results)
    return results


def get_scientific_name(common_name: str):
    options = search_species_options(common_name, max_results=1)
    return options[0][1] if options else None


def _fetch_single_marker(species_name: str, gene: str):
    # Allow a larger default upper bound (10k) so long sequences are not
    # rejected by the SLEN filter when the UI requests them.
    min_len, max_len = MARKER_LENGTH_RANGE.get(gene, (300, 10000))
    search_term = (
        f"{species_name}[Organism] AND {gene}[Gene] "
        f"AND {min_len}[SLEN]:{max_len}[SLEN]"
    )

    def _fetch():
        handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=5,
            sort="relevance",
            timeout=NCBI_TIMEOUT,
        )
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return None

        for seq_id in record["IdList"]:
            handle = Entrez.efetch(
                db="nucleotide",
                id=seq_id,
                rettype="fasta",
                retmode="text",
                timeout=NCBI_TIMEOUT,
            )
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()
            seq_str = str(seq_record.seq)
            if min_len <= len(seq_str) <= max_len:
                print(f"[NCBI] {species_name} {gene}: accepted {len(seq_str)}bp")
                return seq_str
            print(f"[NCBI] {species_name} {gene}: skipped {len(seq_str)}bp")

        return None

    return _fetch_with_retry(_fetch)


def get_marker_sequence(species_name: str, gene: str):
    existing = db_get_marker_sequence(species_name, gene)
    if existing:
        return existing

    seq = _fetch_single_marker(species_name, gene)
    if seq:
        save_marker_sequence(species_name, gene, seq)
    return seq


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
