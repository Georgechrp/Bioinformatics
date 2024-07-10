datasetA = ['AAAATCGATGCTTATGGACTGATTCATTCGTAAACTT', 'TCAAATTCACGCTTATGGTCTCATTTATTCTAGAGCG', 'TAATTGACGCGTATGGACTCATTTACTCGTAATTTT', 'GACTTGACGCTTTGGACTCGTTTATTCGAATGGCG', 'CAATGACGCTTATGACTCATTTCTTCGTAAAGC', 'ACAACTGACGTTATGGACTCATTATTCGTAGTG', 'AGACATTGACTCTTATGGACTCATTTAATCGTAGCATG', 'AACATTGAGGCTTATGGACTATTATTCGTAAAGGC', 'TAAATGGAGCTTATGGCCTCATTCATTCGTATTAA', 'CAGAAGTGAAGCTTATGGCCTCATTTATTCGTAGGCTA', 'CCCAATGAAGCTTATCGACTCATTTATTCGTAAAGAA', 'TAATTGACGCTATGAACTCATTTATGCGTATCTA', 'AAATTAGCTTATGGACACATTTCTTCGTATCTGG', 'AAAATGACGCTTTGGACTGATTTATTCGCAAGT', 'TAAAATTCACGCTTCTGGCTCATTTATTCTACTGTC']


import numpy as np

def global_alignment(A, B, alpha=2):
    len_A = len(A)
    len_B = len(B)
    
    # Initialize the scoring and traceback matrices
    scoring_matrix = np.zeros((len_A + 1, len_B + 1))
    traceback_matrix = np.zeros((len_A + 1, len_B + 1), dtype='object')
    
    # Initialize the first row and column
    for i in range(1, len_A + 1):
        scoring_matrix[i][0] = scoring_matrix[i-1][0] - alpha
        traceback_matrix[i][0] = 'up'
    
    for j in range(1, len_B + 1):
        scoring_matrix[0][j] = scoring_matrix[0][j-1] - alpha
        traceback_matrix[0][j] = 'left'
    
    # Fill the scoring and traceback matrices
    for i in range(1, len_A + 1):
        for j in range(1, len_B + 1):
            match = scoring_matrix[i-1][j-1] + (1 if A[i-1] == B[j-1] else -alpha / 2)
            delete = scoring_matrix[i-1][j] - alpha
            insert = scoring_matrix[i][j-1] - alpha
            
            max_score = max(match, delete, insert)
            scoring_matrix[i][j] = max_score
            
            if max_score == match:
                traceback_matrix[i][j] = 'diag'
            elif max_score == delete:
                traceback_matrix[i][j] = 'up'
            else:
                traceback_matrix[i][j] = 'left'
    
    # Traceback to get the alignment
    align_A = []
    align_B = []
    i = len_A
    j = len_B
    
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 'diag':
            align_A.append(A[i-1])
            align_B.append(B[j-1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'up':
            align_A.append(A[i-1])
            align_B.append('-')
            i -= 1
        elif traceback_matrix[i][j] == 'left':
            align_A.append('-')
            align_B.append(B[j-1])
            j -= 1
    
    # Reverse the alignments as we built them backwards
    align_A = align_A[::-1]
    align_B = align_B[::-1]
    
    return ''.join(align_A), ''.join(align_B)

# Example usage
A = "AAAATCGATGCTTATGGACTGATTCATTCGTAAACTT"
B = "TCAAATTCACGCTTATGGTCTCATTTATTCTAGAGCG"
align_A, align_B = global_alignment(A, B, alpha=2)

print("Alignment A:", align_A)
print("Alignment B:", align_B)

