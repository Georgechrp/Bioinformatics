
from collections import Counter
multiple_alignment = [
    "ACA---ATG",
    "TCAACTATC",
    "ACAC--AGC",
    "ACA---ATC",
    "A-C---ATC"
]
num_rows = len(multiple_alignment)
num_cols = len(multiple_alignment[0]) 
alphabeta = ['A', 'C', 'T', 'G', '-']

def is_conserved_region(align):
    align = align.replace("-", "")
    #print(align, end=" ")
    symbol_count = {}
    for symbol in align:
        if symbol in symbol_count:
            symbol_count[symbol] += 1
        else:
            symbol_count[symbol] = 1
    for i in symbol_count:
        if(symbol_count[i]) >= 4: #threhold 75% 
            return True
    return False

def take_the_column(n):
    col = ""
    for i in range(len(multiple_alignment)):
        col = col + multiple_alignment[i][n]
    return col

def create_hmm_profile():
    print("start", end = " ")
    i_for_conserved_region = 0
    for i in range(num_cols):
        if is_conserved_region(take_the_column(i)):
            i_for_conserved_region= i_for_conserved_region + 1
            print("-->", end = " ")
            print(f"P", {i_for_conserved_region}, end = " ")
    
create_hmm_profile()


def create_Emmision_Prob_table():
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

        print(f"Στήλη {i+1}: {symbol_counts}")

emmision_table = {symbol: [0 for _ in range(num_cols)] for symbol in alphabeta}
create_Emmision_Prob_table()

print("\nΠίνακας Πιθανοτήτων Εκπομπής (Emmision Probability Table):")
for symbol in alphabeta:
    print(f"{symbol}: {emmision_table[symbol]}")