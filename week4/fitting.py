import sys

if len(sys.argv) != 3:
    print("Usage: python3 script.py <target> <query>")
    sys.exit(1)

t1 = sys.argv[1]
q1 = sys.argv[2]

def fitting_alignment(t, q, match=3, mismatch=-3, gap=-2):
    n = len(t)
    m = len(q)
    t = t.upper()
    q = q.upper()

    # Initialize scoring matrix
    S = [[0] * (m + 1) for _ in range(n + 1)]

    #Initialize backtracing matrix
    B = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize first row and column
    for i in range(1, n + 1):
        S[i][0] = 0  # Fitting alignment allows free gaps at the start of t
    for j in range(1, m + 1):
        S[0][j] = S[0][j - 1] + gap

    # Fill in the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if t[i - 1] == q[j - 1]:
                score_diag = S[i - 1][j - 1] + match
            else:
                score_diag = S[i - 1][j - 1] + mismatch
            score_up = S[i - 1][j] + gap
            score_left = S[i][j - 1] + gap

            max_score = max(score_diag, score_up, score_left)

            S[i][j] = max_score
            if max_score == score_diag:
                B[i][j] = 0  # diagonal
            elif max_score == score_up:
                B[i][j] = 1  # up
            else:
                B[i][j] = 2  # left

    # Find the maximum score in the last row
    max_score = float('-inf')
    max_pos = (0, m)
    for i in range(n + 1):
        if S[i][m] > max_score:
            max_score = S[i][m]
            max_pos = (i, m)

    # Backtrack to find the optimal fitting alignment
    align_t = []
    align_q = []
    i, j = max_pos
    while j > 0:
        if B[i][j] == 0:  # diagonal
            align_t.append(t[i - 1])
            align_q.append(q[j - 1])
            i -= 1
            j -= 1
        elif B[i][j] == 1:  # up
            align_t.append(t[i - 1])
            align_q.append('-')
            i -= 1
        else:  # left
            align_t.append('-')
            align_q.append(q[j - 1])
            j -= 1
    # Add any remaining gaps in t
    while i > 0:
        align_t.append(t[i - 1])
        align_q.append('-')
        i -= 1
    
    align_t.reverse()
    align_q.reverse()

    return ''.join(align_t), ''.join(align_q), max_score

align_t, align_q, score = fitting_alignment(t1, q1)
# print(f"Fitting Alignment Score: {score}")
# print(align_t)
# print(align_q)