datasetA = ['AAAATCGATGCTTATGGACTGATTCATTCGTAAACTT', 'TCAAATTCACGCTTATGGTCTCATTTATTCTAGAGCG', 'TAATTGACGCGTATGGACTCATTTACTCGTAATTTT', 'GACTTGACGCTTTGGACTCGTTTATTCGAATGGCG', 'CAATGACGCTTATGACTCATTTCTTCGTAAAGC', 'ACAACTGACGTTATGGACTCATTATTCGTAGTG', 'AGACATTGACTCTTATGGACTCATTTAATCGTAGCATG', 'AACATTGAGGCTTATGGACTATTATTCGTAAAGGC', 'TAAATGGAGCTTATGGCCTCATTCATTCGTATTAA', 'CAGAAGTGAAGCTTATGGCCTCATTTATTCGTAGGCTA', 'CCCAATGAAGCTTATCGACTCATTTATTCGTAAAGAA', 'TAATTGACGCTATGAACTCATTTATGCGTATCTA', 'AAATTAGCTTATGGACACATTTCTTCGTATCTGG', 'AAAATGACGCTTTGGACTGATTTATTCGCAAGT', 'TAAAATTCACGCTTCTGGCTCATTTATTCTACTGTC']


import numpy as np

def global_alignment(A, B, alpha=2):    #20206, 21127 --> a=2
    table = [[0] * (len(B)+1) for _ in range(len(A)+1)]
    gap_penalty = - alpha
    match_score = +1
    mismatch_penalty = - alpha/2

   # Αρχικοποίηση
    for i in range(1, len(A) + 1):
        table[i][0] = i * gap_penalty
    for j in range(1, len(B) + 1):
        table[0][j] = j * gap_penalty

    # Υπολογισμός
    for i in range(1, len(A) + 1):
        for j in range(1, len(B) + 1):
            match = table[i-1][j-1] + (1 if A[i-1] == B[j-1] else mismatch_penalty)
            delete = table[i-1][j] + gap_penalty
            insert = table[i][j-1] + gap_penalty
            table[i][j] = max(match, delete, insert)
    for i in table:
        print(i)



     # Ανάκτηση της βέλτιστης στοίχισης
    align1, align2 = [], []
    i, j = len(A), len(B)
    #align1.append(A[i-1])
    #align2.append(B[j-1])

    while i > 0 and j > 0:
        score_current = table[i][j]
        score_diagonal = table[i-1][j-1]
        score_up = table[i][j-1]
        score_left = table[i-1][j]

        max_score = max(score_diagonal, score_up, score_left)
        if max_score == score_diagonal:
            i -= 1
            j -= 1
            align1.append(A[i-1])
            align2.append(B[j-1])
        elif max_score == score_up:
            j -= 1
            align1.append("-")
            align2.append(B[j-1])
        elif max_score == score_left:
            i -= 1
            align1.append(A[i-1])
            align2.append("-")

        
   

    align1.reverse()
    align2.reverse()

    return ''.join(align1), ''.join(align2)

   



# Παράδειγμα χρήσης
seq1 = "GCT"
seq2 = "AGTAC"
alpha = 2

alignment1, alignment2 = global_alignment(seq1, seq2, alpha)
print("Alignment 1:", alignment1)
print("Alignment 2:", alignment2)