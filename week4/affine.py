import sys

if len(sys.argv) != 3:
    print("Usage: python3 script.py <target> <query>")
    sys.exit(1)

t1 = sys.argv[1]
q1 = sys.argv[2]

NEG_INF = -sys.maxsize

def affine_alignment(t, q, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
    n = len(t)
    m = len(q)
    t = t.upper()
    q = q.upper()

    S = [[0]*(m+1) for _ in range(n+1)]
    E = [[NEG_INF]*(m+1) for _ in range(n+1)]
    F = [[NEG_INF]*(m+1) for _ in range(n+1)]

    # backtrace matrices: 0=diag, 1=left, 2=up
    B_S = [[None]*(m+1) for _ in range(n+1)]
    B_E = [[None]*(m+1) for _ in range(n+1)]
    B_F = [[None]*(m+1) for _ in range(n+1)]

    # initialization
    S[0][0] = 0
    for i in range(1, n+1):
        S[i][0] = gap_open + (i-1)*gap_extend
        F[i][0] = gap_open + (i-1)*gap_extend
        E[i][0] = NEG_INF
        B_F[i][0] = 2  # came from above (F[i-1][0])

    for j in range(1, m+1):
        S[0][j] = gap_open + (j-1)*gap_extend
        E[0][j] = gap_open + (j-1)*gap_extend
        F[0][j] = NEG_INF
        B_E[0][j] = 1  # came from left (E[0][j-1])

    # fill
    for i in range(1, n+1):
        for j in range(1, m+1):
            # E[i][j]: gap in target
            open_E = S[i][j-1] + gap_open
            extend_E = E[i][j-1] + gap_extend
            if open_E > extend_E:
                E[i][j] = open_E
                B_E[i][j] = 0
            else:
                E[i][j] = extend_E
                B_E[i][j] = 1

            # F[i][j]: gap in query
            open_F = S[i-1][j] + gap_open
            extend_F = F[i-1][j] + gap_extend
            if open_F > extend_F:
                F[i][j] = open_F
                B_F[i][j] = 0
            else:
                F[i][j] = extend_F
                B_F[i][j] = 2

            # S[i][j]: match/mismatch
            if t[i-1] == q[j-1]:
                diag = S[i-1][j-1] + match
            else:
                diag = S[i-1][j-1] + mismatch

            max_val = max(diag, E[i][j], F[i][j])
            S[i][j] = max_val
            if max_val == diag:
                B_S[i][j] = 0
            elif max_val == E[i][j]:
                B_S[i][j] = 1
            else:
                B_S[i][j] = 2

    # backtrace
    align_t = []
    align_q = []
    i, j = n, m
    state = 'S'

    while i > 0 or j > 0:
        if state == 'S':
            bt = B_S[i][j]
            if bt == 0:
                align_t.append(t[i-1])
                align_q.append(q[j-1])
                i -= 1
                j -= 1
            elif bt == 1:
                state = 'E'
            else:
                state = 'F'

        elif state == 'E':
            bt = B_E[i][j]
            align_t.append('-')
            align_q.append(q[j-1])
            j -= 1
            if bt == 0:
                state = 'S'

        elif state == 'F':
            bt = B_F[i][j]
            align_t.append(t[i-1])
            align_q.append('-')
            i -= 1
            if bt == 0:
                state = 'S'

    align_t.reverse()
    align_q.reverse()
    return ''.join(align_t), ''.join(align_q), S[n][m]


align_t, align_q, score = affine_alignment(t1, q1)
# print(f"Affine Alignment Score: {score}")
# print(align_t)
# print(align_q) 