# =============================================================================
# dna_matcher/alg_smith_waterman.py — Smith-Waterman Local Alignment Algorithm
# =============================================================================

def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    dp = [[0] * cols for _ in range(rows)]
    best_score = 0
    best_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            score_up = dp[i-1][j] + gap
            score_left = dp[i][j-1] + gap

            dp[i][j] = max(0, diagonal, score_up, score_left)

            if dp[i][j] > best_score:
                best_score = dp[i][j]
                best_pos = (i, j)

    aligned = []
    i, j = best_pos

    while i > 0 and j > 0 and dp[i][j] > 0:
        current = dp[i][j]
        diag_val = dp[i-1][j-1]
        up_val = dp[i-1][j]

        if seq1[i-1] == seq2[j-1] and current == diag_val + match:
            aligned.append(seq1[i-1])
            i -= 1
            j -= 1
        elif current == diag_val + mismatch:
            i -= 1
            j -= 1
        elif current == up_val + gap:
            i -= 1
        else:
            j -= 1

    aligned_str = "".join(reversed(aligned))
    max_possible = match * min(len(seq1), len(seq2))
    percentage = (best_score / max_possible * 100) if max_possible > 0 else 0.0

    return aligned_str, best_score, percentage
