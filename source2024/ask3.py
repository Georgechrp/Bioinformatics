
from collections import Counter

multiple_alignment = [
    "ACA---ATG",
    "TCAACTATC",
    "ACAC--AGC",
    "ACA---ATC",
    "A-C---ATC"
]

threshold = 70
num_rows = len(multiple_alignment)
num_cols = len(multiple_alignment[0]) 
alphabeta = ['A', 'C', 'T', 'G', '-']

def is_conserved_region(align):
    # Αφαιρούμε τα παύλες για να μετρήσουμε μόνο τα σύμβολα A, C, T, G
    clean_align = align.replace("-", "")
    symbol_count = Counter(clean_align)
    
    for symbol, count in symbol_count.items():
        if (count / len(align)) * 100 >= threshold: 
            return True
    return False

def take_the_column(n):
    col = ""
    for i in range(len(multiple_alignment)):
        col = col + multiple_alignment[i][n]
    return col

def create_match_states():
    print("start", end = " ")
    i_for_conserved_region = 0
    for i in range(num_cols):
        if is_conserved_region(take_the_column(i)):
            i_for_conserved_region= i_for_conserved_region + 1
            print(f"-{find_match_states(take_the_column(i-1), take_the_column(i))}->", end = " ")
            print(f"[P{i_for_conserved_region}]", end = " ")
    print("-1-> end", end = " ")

def create_Emmision_Prob_table():
    emmision_table = {symbol: [0 for _ in range(num_cols)] for symbol in alphabeta}
    for i in range(num_cols):
        col = take_the_column(i)
        col_no_gaps = col.replace("-", "")
        symbol_counts = Counter(col_no_gaps)
        total_symbols = len(col_no_gaps)
        
        for symbol in alphabeta:
            if symbol in symbol_counts:
                prob = symbol_counts[symbol] / total_symbols
            else:
                prob = 0
            emmision_table[symbol][i] = prob

        #print(f"Στήλη {i+1}: {symbol_counts}")
    print("\nΠίνακας Πιθανοτήτων Εκπομπής (Emmision Probability Table):")
    for symbol in alphabeta:
        print(f"{symbol}: {emmision_table[symbol]}")


def find_deletions(a, b):
    deletions = 0
    if(is_conserved_region(a) and is_conserved_region(b)):
        for i in range(len(a)):
            if a[i]!='-' and b[i]=="-":
                deletions = deletions +  1/len(a)
    return deletions

def find_insertions(a, b):
    insertions = 0
    if(is_conserved_region(a) and is_conserved_region(b)==False):
        for i in range(len(a)):
            if a[i]!='-' and b[i]!="-": 
                insertions = insertions +  1/len(a)
    return insertions

def find_match_states(a, b):
    return 1 - (find_insertions(a, b) + find_deletions(a, b))

def create_Transition_Prob_table():
    match_states = []
    insertion_states = []
    deletion_states = []
    state_count = 1

    for i in range(num_cols - 1):
        current_state = f"P{state_count}"
        next_state_match = f"P{state_count + 1}"
        next_state_insert = f"I{state_count + 1}"
        next_state_delete = f"D{state_count + 1}"

        deletions = find_deletions(take_the_column(i), take_the_column(i + 1))
        insertions = find_insertions(take_the_column(i), take_the_column(i + 1))
        match = find_match_states(take_the_column(i), take_the_column(i + 1))
        
        if is_conserved_region(take_the_column(i)):
            if is_conserved_region(take_the_column(i + 1)):
                match_states.append((current_state, next_state_match, match))
                state_count += 1
            else:
                if insertions > 0:
                    insertion_states.append((current_state, next_state_insert, insertions))
        else:
            if deletions > 0 or take_the_column(i).count('-') > 0:
                deletion_states.append((current_state, next_state_delete, deletions))
    
    print("\nΠίνακας Μεταβάσεων (Transition Probability Table):")
    if match_states:
        print("Match States:   ", [f"{s[0]} -> {s[1]}: {s[2]:.4f}" for s in match_states])
    if insertion_states:
        print("Insertions:     ", [f"{s[0]} -> {s[1]}: {s[2]:.4f}" for s in insertion_states])
    if deletion_states:
        print("Deletions:      ", [f"{s[0]} -> {s[1]}: {s[2]:.4f}" for s in deletion_states])


#for i in range(num_cols):
    #print(i, take_the_column(i), is_conserved_region(take_the_column(i)))

create_match_states()
#create_Emmision_Prob_table()

'''Transition_Table = [[0 for _ in range(10)] for _ in range(10)]
for row in Transition_Table:
    print(row)'''

#create_Transition_Prob_table()
