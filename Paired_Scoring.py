import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import math
import numpy as np
from rdkit.Chem import rdPartialCharges
from collections import defaultdict
from itertools import product
import pickle
import gc
import matplotlib.pyplot as plt
import seaborn as sns

# Read SMILES strings from file
smiles = []
with open(r'smiles.txt', 'r') as fp:
    for line in fp:
        x = line.strip()
        smiles.append(x)

# Function to calculate Gasteiger charges for each atom
def calculate_gasteiger_charges(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    rdPartialCharges.ComputeGasteigerCharges(molecule)
    gasteiger_charges = defaultdict(list)
    for atom in molecule.GetAtoms():
        atom_symbol = atom.GetSymbol()
        atom_charge = atom.GetProp("_GasteigerCharge")
        gasteiger_charges[atom_symbol].append(float(atom_charge))
    return gasteiger_charges

# Collect all charges by atom type
all_charges = defaultdict(list)
for smile in smiles:
    charges = calculate_gasteiger_charges(smile)
    for atom, charge in charges.items():
        all_charges[atom].extend(charge)

# Remove NaN and inf values
atom_dict = {}
for atom, charges in all_charges.items():
    charges = [x for x in charges if not math.isnan(x) and not math.isinf(x)]
    atom_dict[atom] = charges
all_gast_charges = [charge for charges in atom_dict.values() for charge in charges]

plt.figure(figsize=(8, 6))
plt.hist(all_gast_charges, bins=50, color='skyblue', edgecolor='black', density=True)
plt.xlabel("Gasteiger Charge")
plt.ylabel("Density")
plt.title("Gasteiger Charge Distribution (Histogram Only)")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 6))
sns.kdeplot(all_gast_charges, fill=True, color="skyblue", linewidth=1.5)
plt.xlabel("Gasteiger Charge")
plt.ylabel("Density")
plt.title("Gasteiger Charge Density Plot")
sns.despine()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 6))
sns.histplot(all_gast_charges, bins=50, kde=True, stat="density", color="skyblue", edgecolor="black")
plt.xlabel("Gasteiger Charge")
plt.ylabel("Density")
plt.title("Gasteiger Charge Distribution (Histogram + KDE)")
sns.despine()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

# Build list of atom type pairs
elements = atom_dict.keys()
atom_pairs = list(product(elements, repeat=2))

# Compute pairwise charge differences using chunking to avoid memory errors
def compute_chunked_combinations(charges, chunk_size=1000):
    diffs = []
    n = len(charges)
    for i in range(0, n, chunk_size):
        chunk_i = charges[i:i+chunk_size]
        for j in range(i+1, n, chunk_size):
            chunk_j = charges[j:j+chunk_size]
            diffs.extend([abs(a - b) for a in chunk_i for b in chunk_j])
        chunk_j = chunk_i
        for x in range(len(chunk_j)):
            for y in range(x+1, len(chunk_j)):
                diffs.append(abs(chunk_j[x] - chunk_j[y]))
    return diffs

charge_diffs = {}
for atom1, atom2 in atom_pairs:
    if atom1 in atom_dict and atom2 in atom_dict:
        if atom1 == atom2:
            diffs = compute_chunked_combinations(atom_dict[atom1])
        else:
            diffs = []
            list1 = atom_dict[atom1]
            list2 = atom_dict[atom2]
            len1, len2 = len(list1), len(list2)
            for i in range(0, len1, 1000):
                chunk1 = list1[i:i+1000]
                for j in range(0, len2, 1000):
                    chunk2 = list2[j:j+1000]
                    diffs.extend([abs(a - b) for a in chunk1 for b in chunk2])
        charge_diffs[(atom1, atom2)] = diffs
        gc.collect()

# Build dictionary of raw diffs for scoring
pair_dif = {}
for pair, diffs in charge_diffs.items():
    pair_dif[pair] = diffs

# Build scoring dictionary with log2 probabilities
score_pair = {}
for pair in pair_dif.keys():
    gast = pair_dif[pair]
    total_count = len(gast)
    if total_count == 0:
        continue  
    score = {}
    for i in range(0, 30):
        test = round(i * 0.1, 1)
        count = len([x for x in gast if x >= test])
        prob = count / total_count
        score[test] = np.log2(prob) if prob > 0 else float('-inf')
    score_pair[pair] = score

# Replace -inf with most negative non-inf value
def replace_neg_inf_with_max(dictionary):
    most_neg_value = min(
        value for inner_dict in dictionary.values()
        for value in inner_dict.values()
        if value != float('-inf')
    )
    for inner_dict in dictionary.values():
        for key, value in inner_dict.items():
            if value == float('-inf'):
                inner_dict[key] = most_neg_value
    return dictionary

# Apply -inf fix
result = replace_neg_inf_with_max(score_pair)

# Save to pickle file
with open('score_pair_doubles.pkl', 'wb') as f:
    pickle.dump(result, f)