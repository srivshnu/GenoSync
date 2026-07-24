from .alg_hirschberg import hirschberg_reconstruct
from .alg_smith_waterman import smith_waterman


def calculate_similarity_and_alignment(seq1: str, seq2: str, max_len: int = 1000):
    if not seq1 or not seq2:
        return {
            "hirschberg_pct": 0.0,
            "hirschberg_match": "",
            "sw_pct": 0.0,
            "sw_match": "",
            "sw_score": 0,
            "compare_len": 0,
        }

    compare_len = min(len(seq1), len(seq2), max_len)
    s1 = seq1[:compare_len]
    s2 = seq2[:compare_len]

    matched_global = hirschberg_reconstruct(s1, s2)
    hirschberg_pct = (len(matched_global) / compare_len) * 100

    matched_local, sw_score, sw_pct = smith_waterman(s1, s2)

    return {
        "hirschberg_pct": hirschberg_pct,
        "hirschberg_match": matched_global,
        "sw_pct": sw_pct,
        "sw_match": matched_local,
        "sw_score": sw_score,
        "compare_len": compare_len,
    }
