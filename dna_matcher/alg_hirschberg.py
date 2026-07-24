# =============================================================================
# dna_matcher/alg_hirschberg.py — Hirschberg's Divide-and-Conquer LCS Reconstruction
# =============================================================================

def get_lcs_row(seq1, seq2):
    """
    Sliding-window LCS scorer. Returns the entire final row of DP scores
    instead of just the last value — Hirschberg needs every column to
    locate the optimal split point across seq2.
    """
    prev_row = [0] * (len(seq2) + 1)

    for c1 in seq1:
        curr_row = [0] * (len(seq2) + 1)
        for j, c2 in enumerate(seq2):
            if c1 == c2:
                curr_row[j + 1] = prev_row[j] + 1
            else:
                curr_row[j + 1] = max(curr_row[j], prev_row[j + 1])
        prev_row = curr_row

    return prev_row


def hirschberg_reconstruct(seq1, seq2):
    """
    Recursively reconstructs the LCS string using divide-and-conquer.
    """
    if len(seq1) == 0 or len(seq2) == 0:
        return ""
    if len(seq1) == 1:
        return seq1 if seq1 in seq2 else ""

    mid = len(seq1) // 2
    seq1_left = seq1[:mid]
    seq1_right = seq1[mid:]

    score_left = get_lcs_row(seq1_left, seq2)
    score_right = get_lcs_row(seq1_right[::-1], seq2[::-1])
    score_right.reverse()

    split_index = max(
        range(len(seq2) + 1),
        key=lambda j: score_left[j] + score_right[j]
    )

    left_result = hirschberg_reconstruct(seq1_left, seq2[:split_index])
    right_result = hirschberg_reconstruct(seq1_right, seq2[split_index:])

    return left_result + right_result
