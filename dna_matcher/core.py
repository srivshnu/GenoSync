# dna_matcher/core.py
import os
from pathlib import Path
from collections import OrderedDict
from Bio import Entrez
from dotenv import load_dotenv

from .fetch import fetch_marker_sequence, fetch_with_retry
from .db import (
    init_db,
    get_marker_sequence as db_get_marker_sequence,
    save_marker_sequence,
    get_search_results,
    save_search_results,
)

# Explicitly resolve the path to the .env file in the root directory
env_path = Path(__file__).resolve().parent.parent / ".env"
load_dotenv(dotenv_path=env_path)

Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")

if not Entrez.email:
    print("[DEBUG] WARNING: NCBI_EMAIL is empty. The .env file was not loaded or is missing the variable.")

init_db()

MARKER_CHAINS = {
    "animal":  ["COI", "16S", "12S", "cytb"],
    "plant":   ["rbcL", "matK", "ITS"],
    "fungus":  ["ITS", "LSU", "SSU"],
    "default": ["COI", "16S", "ITS", "rbcL", "12S"],
}

def search_species_options(common_name: str, max_results: int = 6):
    def _fetch():
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
            except Exception as e:
                print(f"[DEBUG] NCBI esearch error for '{q}': {e}") 
                record = {"IdList": []}

            if record and record.get("IdList"):
                try:
                    handle = Entrez.efetch(db="taxonomy", id=record["IdList"], retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()
                except Exception as e:
                    print(f"[DEBUG] NCBI efetch error: {e}")
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

    cached = get_search_results(common_name)
    if cached:
        return cached[:max_results]

    results = fetch_with_retry(_fetch) or []
    if results:
        save_search_results(common_name, results)
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
    failed_species=set()
    for gene in markers:
        sequences = OrderedDict()
        for species_name in species_names:
            seq = get_marker_sequence(species_name, gene)
            if not seq:
                print(f"[NCBI] {species_name} missing {gene}")
                failed_species.add(species_name)
                break
            sequences[species_name] = seq

        if len(sequences) == len(species_names):
            print(f"[NCBI] All species matched on marker: {gene}")
            return {"marker": gene, "sequences": sequences}

    return {"error":list(failed_species)}