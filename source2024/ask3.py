#αυτός ο κώδικας δεν παίζει να είναι σωστός
import numpy as np
from collections import defaultdict

# Διαβάζουμε τα αποτελέσματα της πολλαπλής στοίχισης από το αρχείο 'multiple_alignment_result.txt'
file_name = 'auxiliary2024/multiple_alignment_result.txt'
try:
    with open(file_name, 'r') as file:
        lines = file.readlines()
        msa_result = [line.strip() for line in lines]
except FileNotFoundError:
    print(f"The file {file_name} does not exist.")
    msa_result = []
except Exception as e:
    print(f"An error occurred: {e}")
    msa_result = []

# Διαβάζουμε τις ακολουθίες από το αρχείο 'datasetB'
file_name = 'auxiliary2024/datasetB.txt'
try:
    with open(file_name, 'r') as file:
        lines = file.readlines()
        datasetB = [line.strip() for line in lines]
except FileNotFoundError:
    print(f"The file {file_name} does not exist.")
    datasetB = []
except Exception as e:
    print(f"An error occurred: {e}")
    datasetB = []

# Δημιουργούμε το προφίλ HMM από το αποτέλεσμα της MSA
def build_hmm_profile(msa):
    n = len(msa[0])
    counts = defaultdict(lambda: defaultdict(int))
    for seq in msa:
        for i, char in enumerate(seq):
            counts[i][char] += 1

    states = ['M', 'I', 'D']
    emission_probs = {state: [{} for _ in range(n)] for state in states}
    transition_probs = {state: defaultdict(float) for state in states}
    
    for i in range(n):
        total = sum(counts[i].values())
        for char in 'ACGT-':
            emission_probs['M'][i][char] = counts[i][char] / total if char in counts[i] else 1e-10
            emission_probs['I'][i][char] = 1e-10
            emission_probs['D'][i][char] = 1e-10

    for state in states:
        for next_state in states:
            transition_probs[state][next_state] = 1e-10

    for i in range(1, n):
        transition_probs['M']['M'] += 1
        transition_probs['M']['I'] += 1
        transition_probs['M']['D'] += 1
        transition_probs['I']['I'] += 1
        transition_probs['D']['D'] += 1

    for state in states:
        total = sum(transition_probs[state].values())
        for next_state in states:
            transition_probs[state][next_state] /= total

    return emission_probs, transition_probs

emission_probs, transition_probs = build_hmm_profile(msa_result)
print(emission_probs)
print("-------------------------------------------------------------------")
print(transition_probs)

# Αλγόριθμος Viterbi για την ευθυγράμμιση των ακολουθιών
def viterbi(obs, emission_probs, transition_probs):
    n = len(obs)
    states = ['M', 'I', 'D']
    V = [{} for _ in range(n + 1)]
    path = {state: [] for state in states}

    for state in states:
        V[0][state] = 0
        path[state] = [state]

    for i in range(1, n + 1):
        V[i] = {}
        newpath = {}

        for state in states:
            if i - 1 < len(emission_probs[state]):
                (prob, prev_state) = max(
                    (V[i - 1][y0] + (np.log(emission_probs[y0][i - 1].get(obs[i - 1], 1e-10))) + np.log(transition_probs[y0][state]), y0)
                    for y0 in states
                )
                V[i][state] = prob
                newpath[state] = path[prev_state] + [state]

        path = newpath

    (prob, state) = max((V[n][state], state) for state in states)
    return prob, path[state]

# Υπολογίζουμε τα σκορ ευθυγράμμισης και τις διαδρομές για τις ακολουθίες στο datasetB
results = []
for sequence in msa_result:
    prob, path = viterbi(sequence, emission_probs, transition_probs)
    results.append((sequence, prob, path))

# Εκτυπώνουμε τα αποτελέσματα
for result in results:
    sequence, prob, path = result
    print(f"Sequence: {sequence}")
    print(f"Alignment Probability: {prob}")
    print(f"Alignment Path: {''.join(path)}\n")
