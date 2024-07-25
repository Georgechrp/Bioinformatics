import random
import sys

# Ensure the console uses utf-8 encoding
sys.stdout.reconfigure(encoding='utf-8')

def random_symbols(n):
    alphabet = ['A', 'C', 'G', 'T'] # Αλφάβητο
    return ''.join(random.choices(alphabet, k=n)) # Επιλέγουμε n τυχαία σύμβολα από το αλφάβητο

def synthesize_string():
    num_symbols = random.randint(1, 3) # Επιλέγουμε τυχαία 1 έως 3 σύμβολα
    symbols = random_symbols(num_symbols) # Παίρνουμε τα τυχαία σύμβολα
    return symbols

def main_code_to_generate_a_string():
    # Δοκιμή με τα δεδομένα patterns
    patterns = ["AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA"]

    string = synthesize_string() # erotima (i-a)

    print(f"We start with the random string (Version 0): {string}")
    version = 1
    for pattern in patterns: # erotima (i-b)
        num_symbols = random.randint(1, 2) # πόσα σύμβολα θα αντικατστήσουμε(το πολύ 2)
        print(f"    -Θα αντικατασταθουν {num_symbols} σύμβολα με το pattern = {pattern}")
        for i in range(1, num_symbols + 1):
            x = int(random.randint(1, len(pattern) - 1 ))
            print(f"        -Αντικατασταση του συμβολου στην θεση {x}")
            print("         -[pattern[:x-1] = " + pattern[:x-1] + "|| pattern[x:] = " + pattern[x:], end = ' ')
            choices = ['A', 'C', 'G', 'T', '']
            choices.remove(pattern[x])  #ετσι ειμαστε σιγουροι οτι δεν θα αντικαταστησει ενα γραμμα με τον εαυτο του 
            extend_string = pattern[:x-1] + random.choice(choices) + pattern[x:] # αντικατάσταση με ένα άλλο τυχαία επιλεγμένο σύμβολο είτε με κενη συμβολοσειρα(διαγραφή)
            print(" extend_string= " + extend_string + "]")
            pattern = extend_string
        string += extend_string
        print(f"New string (Version {version}) (+ {extend_string}): {string}")
        version += 1

    for i in range(1, len(patterns)): # erotima (i-c)
        string = string + random_symbols(random.randint(1, 2))
    return string

strings = []

for i in range(50):
    strings.append(main_code_to_generate_a_string())

with open("FullDataset.fasta", "w") as file:
    file.write("\n".join(strings))
    print(" - - Δημιουργήθηκε το αρχείο FullDataset - - ")

random.shuffle(strings) #ανακατευουμε την λίστα με τα 50 strings

datasetA = strings[:15]
datasetB = strings[15:]

print("DatasetA:", datasetA)
print("DatasetB:", datasetB)


with open("datasetA.fasta", "w", encoding="utf-8") as file:
    file.write("\n".join(datasetA))
    print(" - - Δημιουργήθηκε το αρχείο datasetA - - ")


with open("datasetB.fasta", "w", encoding="utf-8") as file:
    file.write("\n".join(datasetB))
    print(" - - Δημιουργήθηκε το αρχείο datasetB - - ")