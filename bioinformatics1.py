import random
print(random.randint(1, 2))
pattern = "ACGATGS"
x=2

plus = pattern[:x-1] + random.choice(['A', 'C', 'G', 'T', '']) + pattern[x:] # αντικατάσταση με ένα άλλο τυχαία επιλεγμένο σύμβολο είτε με κενη συμβολοσειρα(διαγραφή)

print(plus)