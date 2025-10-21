The only code I had to change to get my python code to compile with codon was this: 
    B[i][j] = (score_diag, score_up, score_left).index(S[i][j]) # 0: diag, 1: up, 2: left
Which I used to cleanly populate the backtracing matrix. In codon, .index doesn't exist, so I had to split it up into a longer if/else chain.

Also, originally I was using float('-inf') in python, but switched to 
    NEG_INF = -sys.maxsize
since ints and floats don't play as nicely in python

I wasn't able to get the affine alignment of the human vs orangutan to run in GitHub CI because it takes so long and/or uses too much memory. It works in codon, though. It is omitted from the CI.

In total, I estimate this assignment took me 8 hours
