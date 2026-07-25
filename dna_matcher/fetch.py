import os
import time
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

load_dotenv()
Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")
NCBI_TIMEOUT = float(os.environ.get("NCBI_TIMEOUT", 20.0))

MARKER_LENGTH_RANGE = {
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


def fetch_with_retry(fn, retries=3, base_delay=1.0):
    for attempt in range(retries):
        try:
            result = fn()
            if result is not None:
                return result
        except Exception as e:
            print(f"[Retry {attempt + 1}/{retries}] {e}")
        time.sleep(base_delay * (2 ** attempt))
    return None


def fetch_marker_sequence(species_name: str, gene: str):
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

    return fetch_with_retry(_fetch)