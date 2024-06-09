import random

# erotima (i)
def random_symbols(n):
    alphabet = ['A', 'C', 'G', 'T'] # Αλφάβητο
    return ''.join(random.choices(alphabet, k=n)) # Επιλέγουμε n τυχαία σύμβολα από το αλφάβητο

def synthesize_string():
    num_symbols = random.randint(1, 3) # Επιλέγουμε τυχαία 1 έως 3 σύμβολα
    symbols = random_symbols(num_symbols) # Παίρνουμε τα τυχαία σύμβολα
    return symbols


if __name__ == "__main__":
    # Δοκιμή με τα δεδομένα patterns
    patterns = ["AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA"]
    string = synthesize_string() # erotima (i-a)
    print(f"We start with string (Version 0): {string}")

    new_string =''
    version = 1
    for pattern in patterns: # erotima (i-b)
        num_symbols = random.randint(1, 3) # πόσα σύμβολα θα αντικατστήσουμε
        for i in range(1, num_symbols + 1):
            x = int(random.randint(1, len(pattern) - 1 ))
            plus = pattern[:x-1] + random.choice(['A', 'C', 'G', 'T', '']) + pattern[x:] # αντικατάσταση με ένα άλλο τυχαία επιλεγμένο σύμβολο είτε με κενη συμβολοσειρα(διαγραφή)
        new_string += plus
        print(f"New string (Version {version}) (+ {plus}): {new_string}")
        version += 1
    for i in range(1, len(patterns)): # erotima (i-c)
        new_string = new_string + random_symbols(random.randint(1, 2))

