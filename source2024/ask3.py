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


print("start", end = " ")

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


i_for_conserved_region = 0
for i in range(num_cols):
    if is_conserved_region(take_the_column(i)):
        i_for_conserved_region= i_for_conserved_region + 1
        print("-->", end = " ")
        print(f"P", {i_for_conserved_region}, end = " ")
    



