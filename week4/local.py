t1 = "AGC"
q1 = "ag"

def local_alignment(t, q, match=3, mismatch=-3, gap=-2):
    n = len(t)
    m = len(q)
    t = t.upper()
    q = q.upper()

    # Initialize scoring matrix
    S = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if t[i - 1] == q[j - 1]:
                score_diag = S[i - 1][j - 1] + match
            else:
                score_diag = S[i - 1][j - 1] + mismatch
            score_up = S[i - 1][j] + gap
            score_left = S[i][j - 1] + gap
            S[i][j] = max(0, score_diag, score_up, score_left)  # Local alignment allows zero

            if S[i][j] > max_score:
                max_score = S[i][j]
                max_pos = (i, j)

    # Backtrack to find the optimal local alignment
    align_t = []
    align_q = []
    i, j = max_pos
    while S[i][j] != 0:
        if S[i][j] == S[i - 1][j - 1] + (match if t[i - 1] == q[j - 1] else mismatch):
            align_t.append(t[i - 1])
            align_q.append(q[j - 1])
            i -= 1
            j -= 1
        elif S[i][j] == S[i - 1][j] + gap:
            align_t.append(t[i - 1])
            align_q.append('-')
            i -= 1
        else:  # S[i][j] == S[i][j - 1] + gap
            align_t.append('-')
            align_q.append(q[j - 1])
            j -= 1

    align_t.reverse()
    align_q.reverse()

    return ''.join(align_t), ''.join(align_q), max_score

align_t, align_q, score = local_alignment(t1, q1)
print(f"Local Alignment Score: {score}")
print(align_t)
print(align_q)