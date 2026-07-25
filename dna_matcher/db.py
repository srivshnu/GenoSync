import json
import os
import sqlite3
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DB_FILE = os.environ.get("DNA_MATCHER_DB", "dna_matcher_cache.sqlite3")
DB_PATH = Path(DB_FILE).resolve()
CACHE_EXPIRY_DAYS = int(os.environ.get("DNA_MATCHER_CACHE_DAYS", 30))


def _connect() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(DB_PATH), timeout=30, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    return conn


def init_db() -> None:
    with _connect() as conn:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS species_search_results (
                common_name TEXT NOT NULL,
                candidate_common TEXT NOT NULL,
                scientific_name TEXT NOT NULL,
                result_rank INTEGER NOT NULL,
                searched_at TEXT NOT NULL,
                PRIMARY KEY (common_name, result_rank)
            );

            CREATE TABLE IF NOT EXISTS marker_sequences (
                species_name TEXT NOT NULL,
                gene TEXT NOT NULL,
                sequence TEXT NOT NULL,
                sequence_length INTEGER NOT NULL,
                fetched_at TEXT NOT NULL,
                PRIMARY KEY (species_name, gene)
            );

            CREATE TABLE IF NOT EXISTS comparison_results (
                species_a TEXT NOT NULL,
                species_b TEXT NOT NULL,
                marker TEXT NOT NULL,
                max_len INTEGER NOT NULL,
                details TEXT NOT NULL,
                computed_at TEXT NOT NULL,
                PRIMARY KEY (species_a, species_b, marker, max_len)
            );
            """
        )


def _normalize_common_name(common_name: str) -> str:
    return common_name.strip().lower()


def _get_cutoff_time() -> str:
    return (datetime.utcnow() - timedelta(days=CACHE_EXPIRY_DAYS)).isoformat()


def get_search_results(common_name: str) -> List[Tuple[str, str]]:
    normalized = _normalize_common_name(common_name)
    with _connect() as conn:
        rows = conn.execute(
            "SELECT candidate_common, scientific_name FROM species_search_results "
            "WHERE common_name = ? AND searched_at >= ? ORDER BY result_rank",
            (normalized, _get_cutoff_time()),
        ).fetchall()
    return [(row["candidate_common"], row["scientific_name"]) for row in rows]


def save_search_results(common_name: str, results: List[Tuple[str, str]]) -> None:
    normalized = _normalize_common_name(common_name)
    now = datetime.utcnow().isoformat()
    with _connect() as conn:
        for rank, (candidate_common, scientific_name) in enumerate(results):
            conn.execute(
                "INSERT OR REPLACE INTO species_search_results "
                "(common_name, candidate_common, scientific_name, result_rank, searched_at) "
                "VALUES (?, ?, ?, ?, ?)",
                (normalized, candidate_common, scientific_name, rank, now),
            )


def get_marker_sequence(species_name: str, gene: str) -> Optional[str]:
    with _connect() as conn:
        row = conn.execute(
            "SELECT sequence FROM marker_sequences "
            "WHERE species_name = ? AND gene = ? AND fetched_at >= ?",
            (species_name, gene, _get_cutoff_time()),
        ).fetchone()
    return row["sequence"] if row else None


def save_marker_sequence(species_name: str, gene: str, sequence: str) -> None:
    now = datetime.utcnow().isoformat()
    with _connect() as conn:
        conn.execute(
            "INSERT OR REPLACE INTO marker_sequences "
            "(species_name, gene, sequence, sequence_length, fetched_at) "
            "VALUES (?, ?, ?, ?, ?)",
            (species_name, gene, sequence, len(sequence), now),
        )


def _ordered_pair(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def get_comparison_result(
    species_a: str,
    species_b: str,
    marker: str,
    max_len: int,
) -> Optional[Dict]:
    a, b = _ordered_pair(species_a, species_b)
    with _connect() as conn:
        row = conn.execute(
            "SELECT details FROM comparison_results "
            "WHERE species_a = ? AND species_b = ? AND marker = ? AND max_len = ? AND computed_at >= ?",
            (a, b, marker, max_len, _get_cutoff_time()),
        ).fetchone()
        if not row:
            return None
        return json.loads(row["details"])


def save_comparison_result(
    species_a: str,
    species_b: str,
    marker: str,
    max_len: int,
    details: Dict,
) -> None:
    a, b = _ordered_pair(species_a, species_b)
    now = datetime.utcnow().isoformat()
    with _connect() as conn:
        conn.execute(
            "INSERT OR REPLACE INTO comparison_results "
            "(species_a, species_b, marker, max_len, details, computed_at) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (a, b, marker, max_len, json.dumps(details), now),
        )