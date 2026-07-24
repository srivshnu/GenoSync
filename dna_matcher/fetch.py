import time
import os
from dotenv import load_dotenv
from Bio import Entrez, SeqIO

load_dotenv()
Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")

MARKER_LENGTH_RANGE = {
    "COI":  (500,  1600),
    "16S":  (400,  1800),
    "12S":  (300,  1100),
    "cytb": (800,  1200),
    "rbcL": (400,  600),
    "matK": (800,  900),
    "ITS":  (400,  800),
    "LSU":  (800,  3500),
    "SSU":  (1600, 2000),
}


def _fetch_with_retry(fn, retries=3, base_delay=1.0):
    for attempt in range(retries):
        try:
            result = fn()
            if result is not None:
                return result
        except Exception:
            pass
        time.sleep(base_delay * (2 ** attempt))
    return None


def fetch_marker_sequence(species_name: str, gene: str):
    min_len, max_len = MARKER_LENGTH_RANGE.get(gene, (300, 5000))
    search_term = (
        f"{species_name}[Organism] AND {gene}[Gene] "
        f"AND {min_len}[SLEN]:{max_len}[SLEN]"
    )

    def _fetch():
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return None

        for seq_id in record["IdList"]:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()
            seq_str = str(seq_record.seq)
            if min_len <= len(seq_str) <= max_len:
                return seq_str
        return None

    return _fetch_with_retry(_fetch)
