The only code I had to change to get my python code to compile with codon was this: 
    B[i][j] = (score_diag, score_up, score_left).index(S[i][j]) # 0: diag, 1: up, 2: left
Which I used to cleanly populate the backtracing matrix. In codon, .index doesn't exist, so I had to split it up into a longer if/else chain.

Also, originally I was using float('-inf') in python, but switched to 
    NEG_INF = -sys.maxsize
since ints and floats don't play as nicely in python