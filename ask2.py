datasetA = ['AAAATCGATGCTTATGGACTGATTCATTCGTAAACTT', 'TCAAATTCACGCTTATGGTCTCATTTATTCTAGAGCG', 'TAATTGACGCGTATGGACTCATTTACTCGTAATTTT', 'GACTTGACGCTTTGGACTCGTTTATTCGAATGGCG', 'CAATGACGCTTATGACTCATTTCTTCGTAAAGC', 'ACAACTGACGTTATGGACTCATTATTCGTAGTG', 'AGACATTGACTCTTATGGACTCATTTAATCGTAGCATG', 'AACATTGAGGCTTATGGACTATTATTCGTAAAGGC', 'TAAATGGAGCTTATGGCCTCATTCATTCGTATTAA', 'CAGAAGTGAAGCTTATGGCCTCATTTATTCGTAGGCTA', 'CCCAATGAAGCTTATCGACTCATTTATTCGTAAAGAA', 'TAATTGACGCTATGAACTCATTTATGCGTATCTA', 'AAATTAGCTTATGGACACATTTCTTCGTATCTGG', 'AAAATGACGCTTTGGACTGATTTATTCGCAAGT', 'TAAAATTCACGCTTCTGGCTCATTTATTCTACTGTC']


def global_alignment(seq1, seq2, alpha):
    n = len(seq1)
    m = len(seq2)
    
    # Αρχικοποίηση πίνακα βαθμολογίας
    score = [[0 for _ in range(m+1)] for _ in range(n+1)]
    #print(f"Score is {score}")
    # Αρχικοποίηση πρώτης στήλης και πρώτης γραμμής
    for i in range(1, n+1):
        score[i][0] = score[i-1][0] - alpha
    for j in range(1, m+1):
        score[0][j] = score[0][j-1] - alpha
    
    # Υπολογισμός πίνακα βαθμολογίας
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -alpha/2)
            delete = score[i-1][j] - alpha
            insert = score[i][j-1] - alpha
            score[i][j] = max(match, delete, insert)
    
    # Ανακατασκευή στοίχισης
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = n, m
    
    while i > 0 or j > 0:
        current_score = score[i][j]
        if i > 0 and j > 0 and current_score == score[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -alpha/2):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and current_score == score[i-1][j] - alpha:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    return aligned_seq1, aligned_seq2, score[n][m]

def multiple_alignment(datasets, alpha):
    # Ξεκινάμε με την πρώτη συμβολοσειρά
    aligned_sequences = [datasets[0]]
    
    # Συγκρίνουμε κάθε επόμενη συμβολοσειρά με το ευθυγραμμισμένο σύνολο
    for i in range(1, len(datasets)):
        current_sequence = datasets[i]
        
        # Ευθυγράμμιση κάθε συμβολοσειράς με τις ήδη ευθυγραμμισμένες
        new_aligned_sequences = []
        for aligned_seq in aligned_sequences:
            aligned1, aligned2, _ = global_alignment(aligned_seq, current_sequence, alpha)
            new_aligned_sequences.append(aligned1)
        
        # Προσθέτουμε την τρέχουσα ευθυγραμμισμένη συμβολοσειρά στο σύνολο
        new_aligned_sequences.append(aligned2)
        aligned_sequences = new_aligned_sequences
    
    return aligned_sequences


'''# Παράδειγμα χρήσης της συνάρτησης
seq1 = "AATTGA"
seq2 = "CGCTTAT"
alpha = 2  # 20206, 21127

aligned_seq1, aligned_seq2, score = global_alignment(seq1, seq2, alpha)
print(f"Aligned Sequences:\n{aligned_seq1}\n{aligned_seq2}")
print(f"Alignment Score: {score}")
'''


'''
# Παράδειγμα χρήσης της συνάρτησης με datasetA
datasetA = ["AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA", "TTCGGA", "GGATC", "ATTGA", "CGCGTA", "GACTT", "TATTG", "TGACG", "CCTGA", "GGCGTA", "TTGGA", "GGATT"]
alpha = 1  # Υποθέτουμε ότι τα ΑΜ καταλήγουν σε περιττό ψηφίο

aligned_sequences = multiple_alignment(datasetA, alpha)

# Εκτύπωση του αποτελέσματος της πολλαπλής στοίχισης
for seq in aligned_sequences:
    print(seq)
'''