from collections import Counter

# Sample multiple alignment data
multiple_alignment = [
    "ACAA--ATG",
    "TCAACTATC",
    "ACAC--AGC",
    "AGA---ATC",
    "ACCG--ATC"
]

multiple_alignment = [
    "ACA---ATG",
    "TCAACTATC",
    "ACAC--AGC",
    "ACA---ATC",
    "A-----ATC"
]

# Define the threshold and initialize variables
threshold = 70
num_rows = len(multiple_alignment)
num_cols = len(multiple_alignment[0])
alphabeta = ['A', 'C', 'T', 'G']
conserved_states = []
deletion_states = []
insertion_states = []
hmmstates = []
states = []

def check_region(align):
    """
    This function checks if a region in the alignment is conserved,
    has deletions, or insertions based on the provided threshold.
    """
    clean_align = align.replace("-", "")
    symbol_count = len(clean_align)
    conserved_regions = False
    deletions = False
    insertions = False

    # Calculate the percentage of letters (non-dash characters)
    if (symbol_count / len(align)) * 100 >= threshold:
        conserved_regions = True
        if symbol_count < len(align):
            deletions = True
    else:
        insertions = True

    return conserved_regions, deletions, insertions

def search_list(lst, element):
    i = 0
    while i < len(lst):
        if lst[i] == element:
            return True
        i += 1
    return False

#print (search_list([1,2,3] , 1))

def take_the_column(n):
    """
    Extracts the nth column from the multiple alignment.
    """
    col = ""
    for i in range(len(multiple_alignment)):
        col += multiple_alignment[i][n]
    return col

def hmm_profile(multiple_alignment):
    for i in range(num_cols):
        conserved, deletion, insertion = check_region(take_the_column(i))
        states.append(i)

        if (conserved): conserved_states.append(i)
        if (deletion): deletion_states.append(i)
        if (insertion): insertion_states.append(i)

    for i in range(num_cols):
        if not (search_list(deletion_states, states[i])) and search_list(conserved_states, states[i]) :
            hmmstates.append("match") 
        elif search_list(deletion_states, states[i]) and search_list(conserved_states, states[i]) :
            hmmstates.append("delete")
        elif search_list(insertion_states, states[i] ):
            hmmstates.append("insert")
    print(conserved_states, deletion_states, insertion_states)
    print(hmmstates)

hmm_profile(multiple_alignment)


#Emmision Propability Table for HMM Profile
def create_Emmision_Prob_table():
    """
    Creates and prints the Emission Probability Table.
    """
    # Initialize the emission table with zeros for conserved columns only
    emission_table = {symbol: [0 for _ in range(len(conserved_states))] for symbol in alphabeta}
    
    conserved_index = 0  # Tracks the index for conserved states
    for i in range(num_cols):
        if search_list(conserved_states, states[i]):
            col = take_the_column(i)
            col_no_gaps = col.replace("-", "")  # Remove gaps
            symbol_counts = Counter(col_no_gaps)
            total_symbols = len(col_no_gaps)
            
            # Calculate probability of each symbol in the conserved column
            for symbol in alphabeta:
                if symbol in symbol_counts:
                    prob = symbol_counts[symbol] / total_symbols
                else:
                    prob = 0
                # Assign probability to the corresponding conserved state index
                emission_table[symbol][conserved_index] = prob
            conserved_index += 1  # Move to the next conserved state

    # Print the Emission Probability Table
    print("\nEmission Probability Table:")
    for symbol in alphabeta:
        print(f"{symbol}: {emission_table[symbol]}")

# Call the function to print the table
create_Emmision_Prob_table()


def create_Transition_Prob_table():
    j = 0
    prob_table = [[0 for _ in range(len(conserved_states)-1)] for _ in range(7)]
    count = [0 for _ in range(num_rows)]

    for i in range(len(hmmstates) - 1):  # Loop through valid column indices
        if hmmstates[i] == "match" and hmmstates[i + 1] == "delete":
            count_letters = take_the_column(i+1)
            count_letters = count_letters.replace("-", "")
            count_letters = len(count_letters)/num_rows
            print(count_letters)
            prob_table[0][j] = round(count_letters, 2)
            prob_table[2][j] = round(1 - count_letters, 2)
            j += 1
        elif (hmmstates[i] == "match" and hmmstates[i + 1] == "insert") or (hmmstates[i] == "delete" and hmmstates[i + 1] == "insert"):
            counter_of_inserts = 0
            prob_counter = 0

            if hmmstates[i] == "delete" and hmmstates[i + 1] == "insert":
                prob_table[6][j] = prob_table[5][j-1]

            i += 1       
            
            while(hmmstates[i] == "insert"): 
                inserts = take_the_column(i)

                for k in range(len(inserts)):
                    if inserts[k] == "-": count[k] += 1 #αν έχει μiα γραμμή μόνο παύλες
                counter_of_inserts += 1
                i += 1
            
            for k in range(len(inserts)):
                if count[k] == 0: prob_counter += 1  # διόρθωση εδώ, χρησιμοποιούμε count αντί για inserts
            
            prob_table[1][j] = round((counter_of_inserts - prob_counter) / num_rows, 2)
            prob_table[0][j] = round(1 - prob_table[1][j], 2)

            if counter_of_inserts > 0:  # Προσθήκη ελέγχου για αποφυγή διαίρεσης με το μηδέν
                prob_table[3][j] = round(prob_counter / counter_of_inserts, 2)  # Αλλαγή εδώ για σωστή αναλογία
                prob_table[4][j] = round(1 - prob_table[3][j], 2)  # Υπολογισμός του συμπληρώματος της πιθανότητας
            else:
                prob_table[3][j] = 0
                prob_table[4][j] = 1


            j += 1
        elif hmmstates[i] == "match" and hmmstates[i + 1] == "match":
            prob_table[0][j] = 1
            j += 1
        elif hmmstates[i] == "delete" and hmmstates[i + 1] == "delete":
            count_letters = take_the_column(i+1)
            count_letters = count_letters.replace("-", "")
            count_letters = len(count_letters)/num_rows
            print(count_letters)
            prob_table[0][j] = round(count_letters, 2)
            prob_table[5][j] = round(1 - count_letters, 2)
            j += 1
        elif hmmstates[i] == "delete" and hmmstates[i + 1] == "match":
            prob_table[6][j] = prob_table[2][j-1]
            prob_table[0][j] = prob_table[0][j-1] 
            j += 1

    print (prob_table)
    print(count)
    
    for i in range(len(prob_table[0])):
        if prob_table[0][i] != 0: print("Transition: M" + str(i + 1) + " to M" + str(i + 2) + ": " + str(prob_table[0][i]))
        if prob_table[2][i] != 0: print("Transition: M" + str(i + 1) + " to D" + str(i + 1) + ": " + str(prob_table[2][i]))
        if prob_table[1][i] != 0: print("Transition: M" + str(i + 1) + " to I" + str(i + 1) + ": " + str(prob_table[1][i]))
        if prob_table[3][i] != 0: print("Transition: I" + str(i + 1) + " to I" + str(i + 1) + ": " + str(prob_table[3][i]))
        if prob_table[4][i] != 0: print("Transition: I" + str(i + 1) + " to M" + str(i + 2) + ": " + str(prob_table[4][i]))
        if prob_table[5][i] != 0: print("Transition: D" + str(i) + " to D" + str(i + 1) + ": " + str(prob_table[5][i]))
        if prob_table[6][i] != 0: print("Transition: D" + str(i) + " to M" + str(i + 2) + ": " + str(prob_table[6][i]))

        print("--------------------------------")

        
create_Transition_Prob_table()




