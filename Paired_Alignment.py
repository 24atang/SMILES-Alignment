import numpy as np
import pubchempy
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import math
import pickle 
import pandas as pd
import time
import seaborn as sns
import matplotlib.pyplot as plt

#variable neccessary for list creation
CHAR = set(['(',')','{','}','[',']','=','@','.','#','+','-','/','\\''0','1','2','3','4','5','6','7','8','9'])
DOUB_ELEM = set(['Na','Mg','Al','Si','Cl','Ca','Cr','Mn', 'Fe','Co','Cu','Zn','Se','Mo','Cd','Sn','Br','As','Pb','Li'])
DOUBLE_TROUB = {'n':['M','Z','S'], 'o':['C','M']}
CHARs = set(['(',')','{','}','[',']','=','@','.','#','+','-','/','\\','0','1','2','3','4','5','6','7','8','9'])
elements = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", 
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
    "Tl", "Pb", "Bi", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", 
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", 
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", 
    "Ts", "Og"
]


#opens dictionary of dictionary scoring matrix
with open('score_pair_doubles.pkl', 'rb') as f:
    pair_dif_score = pickle.load(f)


#function that changes a string of SMILES to a list
def parse_smiles(smile, reverse=False):
    smiles = []
    for i in range(len(smile) - 1):
        elem = smile[i]
        doub = smile[i:i+2]
        
        if elem.isdigit() or elem in CHAR:
            smiles.append(elem)

        elif doub in DOUB_ELEM:
            smiles.append(doub)

        elif elem in ('a','l','g','i','r','e','u','d','s','b'):
            continue

        elif elem in ('n','o',) and smile[i-1] in DOUBLE_TROUB[elem]:
            continue
        
        else:
            smiles.append(elem)

    # If the last character is not processed, process it
    if smile[-1] not in DOUB_ELEM:
        smiles.append(smile[-1])
        
    smiles = [i for i in smiles if i not in CHAR]

    if reverse:
        return "".join(smiles[::-1])
    else:
        return smiles

#returns a list of gastieger charges in a smile, list includes the index of the SMILE
#where atoms are given a value, and characters are give '-'
def gast_smiles(smile):
    mol = Chem.MolFromSmiles(smile)
    AllChem.ComputeGasteigerCharges(mol)
    charges = [round(atom.GetDoubleProp('_GasteigerCharge'), 1) for atom in mol.GetAtoms()]
    
    sm1 = parse_smiles(smile)
    sm1_red = ['-' if i in CHARs else i for i in sm1]

    gast_smile = []
    count = 0
    for i in sm1_red:
        if i == '-':
            gast_smile.append(i)
        else:
            gast_smile.append(charges[count])
            count += 1

    return(gast_smile)

#scoring function that returns score. characters are either matched or mismatched,
#and atoms are given a score based of the scoring matrix. alpha and beta are the 
#atoms or characters, and gast1 and gast2 are gasteiger charges (or '-' if characters)
def score(alpha, beta , gast1, gast2):
    alpha, beta = alpha.capitalize() , beta.capitalize()
    
    if alpha in elements and beta in elements:
        dif = round(abs(gast1-gast2), 1)
        try:
            return(pair_dif_score[alpha, beta][dif])
        except:
            return(pair_dif_score[beta, alpha][dif])
    
    elif alpha.isdigit() and beta.isdigit():
        return(match_score)
    
    elif alpha == beta:
        return(match_score)
    
    else: 
        return(mis_match_score)
    
from collections import deque

def reconstruct_smiles(modified_smiles, original_smiles):
    modified_queue = deque(modified_smiles)
    original_queue = deque(original_smiles)

    result = []

    while modified_queue:
        char = modified_queue.popleft()

        # Directly append dashes
        if char == '-':
            result.append(char)
            continue

        # For atoms: match with original and append all characters till the next atom
        while original_queue:
            orig_char = original_queue.popleft()
            result.append(orig_char)
            if orig_char == char:
                break

    # Append any remaining characters from the original queue
    result.extend(original_queue)

    return ''.join(result)

#alignment function based on Needelman-Wunsch
def align(seq1, seq2):
    
    seq1_orig = seq1
    seq2_orig = seq2
    
    gast1 , gast2 = gast_smiles(seq1) , gast_smiles(seq2)
    seq1, seq2 = parse_smiles(seq1) , parse_smiles(seq2)
    len1 , len2 = len(seq1) , len(seq2)
    
    M = np.zeros((len1+1, len2+1))
    X = np.zeros((len1+1, len2+1))
    Y = np.zeros((len1+1, len2+1))
    
    M[0, :] = -np.inf
    M[:, 0] = -np.inf
    M[0, 0] = 0

    X[:, 0] = gap_open_penalty + gap_extend_penalty * np.arange(len1+1)
    X[0, :] = -np.inf

    Y[0, :] = gap_open_penalty + gap_extend_penalty * np.arange(len2+1)
    Y[:, 0] = -np.inf
    
    for i in range(1, len1+1):
            for j in range(1, len2+1):
                 # Calculate match/mismatch score using the scoring matrix
                match = score(seq1[i-1], seq2[j-1], gast1[i-1] , gast2[j-1])
                M[i, j] = max(M[i-1, j-1] + match, X[i-1, j-1] + match, Y[i-1, j-1] + match)
                X[i, j] = max(M[i-1, j] + gap_open_penalty, X[i-1, j] + gap_extend_penalty, gap_open_penalty + gap_extend_penalty + Y[i-1,j])
                Y[i, j] = max(M[i, j-1] + gap_open_penalty, Y[i, j-1] + gap_extend_penalty)
                
                
    seq1_align, seq2_align = '', ''
    i, j = len1, len2
    while i > 0 or j > 0:
        if i > 0 and j > 0 and max(M[i, j], X[i, j], Y[i, j]) == M[i, j]:
            seq1_align += str(seq1[i-1])
            seq2_align += str(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and max(M[i, j], X[i, j], Y[i, j]) == X[i, j]:
            seq1_align += str(seq1[i-1])
            seq2_align += '-'
            i -= 1
        elif j > 0 and max(M[i, j], X[i, j], Y[i, j]) == Y[i, j]:
            seq1_align += '-'
            seq2_align += str(seq2[j-1])
            j -= 1
        else:
            break
            
    smiles_out1 = reconstruct_smiles(seq1_align[::-1], seq1_orig)
    smiles_out2 = reconstruct_smiles(seq2_align[::-1], seq2_orig)
    
    
    return(seq1_align[::-1] , seq2_align[::-1], M[len1,len2],smiles_out1,smiles_out2)

#Levenshtein similairty for comparing the similarity of the truth and observed
def levenshtein_similarity(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_similarity(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    distance = previous_row[-1]
    max_len = max(len(s1), len(s2))
    return 1 - distance / max_len


def string_similarity(truth, experimental):
    """
    Calculate the percentage similarity between two strings based on exact match. Extra characters in the 
    experimental string count against similarity.

    Parameters:
    - truth: The ground truth string.
    - experimental: The experimental string to compare against the truth.

    Returns:
    - Percentage similarity between the two strings.
    """
    max_len = max(len(truth), len(experimental))
    matches = sum(t == e for t, e in zip(truth, experimental))
    
    similarity_percentage = (matches / max_len)
    return similarity_percentage

def tanimoto_aligned_characters(align, truth):

    intersection = 0
    union = 0

    for c1, c2 in zip(align, truth):
        # Skip positions where both are gaps
        if c1 == '-' and c2 == '-':
            continue

        # Count matching positions
        if c1 == c2:
            intersection += 1

        # Count union for all non double-gap positions
        union += 1

    return intersection / union if union > 0 else 0.0


#canonized SMILES
canonized = ('[O-]C(=O)C(=O)CC([O-])=O',
'[O-]C(=O)C(O)(CC([O-])=O)CC([O-])=O',
'[O-]C(=O)C(CC([O-])=O)C(O)C([O-])=O',
'[O-]C(=O)CCC(=O)C([O-])=O',
'[O-]C(=O)CCC([O-])=O',
'[O-]C(=O)C=CC([O-])=O',
'[O-]C(=O)C(O)CC([O-])=O')


#Hand aligned truth set
truth_set = ((('OCOCOCCOO', 'OCOCOCCOO'),
  ('OCOCO----CCOO', 'OCOCOCCOOCCOO'),
  ('OCO-CO---C-COO', 'OCOCC-COOCOCOO'),
  ('OCOCO-C-COO', 'OCOC-CCOCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOCOCCOO')),
 [('OCOCOCCOOCCOO', 'OCOCO----CCOO'),
  ('OCOCOCCOOCCOO', 'OCOCOCCOOCCOO'),
  ('OCOCOCCOOC-COO', 'OCOC-CCOOCOCOO'),
  ('OCOCOCCOOC-COO', 'OCOC-C---COCOO'),
  ('OCOCOCCOOCCOO', 'OCOC----CCOO'),
  ('OCOCOCCOOCCOO', 'OCOC-----CCOO'),
  ('OCOCOCCOOCCOO', 'OCOCO----CCOO')],
 [('OCOCC-COOCOCOO', 'OCO-CO---C-COO'),
  ('OCOC-CCOOCOCOO', 'OCOCOCCOOC-COO'),
  ('OCOCCCOOCOCOO', 'OCOCCCOOCOCOO'),
  ('OCOCCCOOCOCOO', 'OCOCC---COCOO'),
  ('OCOCCCOOCOCOO', 'OCO-C---C-COO'),
  ('OCOCCCOOCOCOO', 'OCO-C---C-COO'),
  ('OCOC-CCOOCOCOO', 'OCOCOCCOO-----')],
 (('OCOC-CCOCOO', 'OCOCO-C-COO'),
  ('OCOC-C---COCOO', 'OCOCOCCOOC-COO'),
  ('OCOCC---COCOO', 'OCOCCCOOCOCOO'),
  ('OCOCCCOCOO', 'OCOCCCOCOO'),
  ('OCOCCCOCOO', 'OCOCC--COO'),
  ('OCOCCCOCOO', 'OCOCC--COO'),
  ('OCOC-CCOCOO', 'OCOCOC--COO')),
 (('OCOC-CCOO', 'OCOCOCCOO'),
  ('OCO---C--CCOO', 'OCOCOCCOOCCOO'),
  ('OCO-C---C-COO', 'OCOCCCOOCOCOO'),
  ('OCOCC--COO', 'OCOCCCOCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOC-CCOO', 'OCOCOCCOO')),
 (('OCOC-CCOO', 'OCOCOCCOO'),
  ('OCO---C--CCOO', 'OCOCOCCOOCCOO'),
  ('OCO-C---C-COO', 'OCOCCCOOCOCOO'),
  ('OCOCC--COO', 'OCOCCCOCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOC-CCOO', 'OCOCOCCOO')),
 (('OCOCOCCOO', 'OCOCOCCOO'),
  ('OCOCO----CCOO', 'OCOCOCCOOCCOO'),
  ('OCOCOCCOO-----', 'OCOC-CCOOCOCOO'),
  ('OCOCOC--COO', 'OCOC-CCOCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOCOCCOO')))

import pandas as pd
import time
alignment_log = open("Paired_Alignment_output.txt", "w")

start = time.time()
results = []
run = 1

# Loop over all combinations of parameters
for gap_open_penalty in range(-5, 6):
    for gap_extend_penalty in range(-5, 6):
        if gap_open_penalty <= gap_extend_penalty:
            similarities = []
            exact = []
            for j in range(0, len(canonized)):
                for i in range(0, len(canonized)):
                    # aligns
                    a = align(canonized[j], canonized[i])
                    alignment_log.write(f"\nParameter Test - Gap Open: {gap_open_penalty}, Gap Extend: {gap_extend_penalty}\n")
                    alignment_log.write(f"Original 1: {canonized[j]}\nOriginal 2: {canonized[i]}\n")
                    alignment_log.write(f"Aligned 1:  {a[3]}\nAligned 2:  {a[4]}\n")
                    alignment_log.write(f"Score: {a[2]}\n")
                    alignment_log.write("-" * 60 + "\n")
                    # compares truth with levenshtein
                    similarity1 = levenshtein_similarity(a[0], truth_set[j][i][0])
                    similarity2 = levenshtein_similarity(a[1], truth_set[j][i][1])
                    similarities.extend((similarity1, similarity2))

                    # compares with exact
                    similarity_a = string_similarity(a[0], truth_set[j][i][0])
                    similarity_b = string_similarity(a[1], truth_set[j][i][1])
                    exact.extend((similarity_a, similarity_b))

            levenshtein_average = sum(similarities) / len(similarities) 
            exact_average = sum(exact) / len(exact) 
            print(f"Gap Open: {gap_open_penalty}, Gap Extend: {gap_extend_penalty}, Levenshtein Avg: {levenshtein_average:.4f}, Exact Avg: {exact_average:.4f}")
            results.append([gap_open_penalty, gap_extend_penalty, levenshtein_average, exact_average])
            
df = pd.DataFrame(results, columns=[ 'Gap Open Penalty', 'Gap Extend Penalty', 'Levenshtein Average', 'Exact Average'])
df.to_csv("Paired_Alignment_parameter_results.csv", index=False)

best_params = df.sort_values(by='Levenshtein Average', ascending=False).iloc[0]
gap_open_penalty = best_params['Gap Open Penalty']
gap_extend_penalty = best_params['Gap Extend Penalty']
print(f"\nBest parameters found:\nGap Open = {gap_open_penalty}, Gap Extend = {gap_extend_penalty}")

optalignment_log = open("Paired_optimal_alignment_output.txt", "w")
optalignment_log.write(f"Optimal Parameters Used:\nGap Open: {gap_open_penalty}, Gap Extend: {gap_extend_penalty}\n\n")

for j in range(len(canonized)):
    for i in range(len(canonized)):
        a = align(canonized[j], canonized[i])

        # Optionally compute similarity
        levenshtein1 = levenshtein_similarity(a[0], truth_set[j][i][0])
        levenshtein2 = levenshtein_similarity(a[1], truth_set[j][i][1])

        # Write aligned output
        smiles1 = canonized[j]
        smiles2 = canonized[i]
        aligned1 = a[3]
        aligned2 = a[4]
        truth1 = truth_set[j][i][0]
        truth2 = truth_set[j][i][1]
        score_val = a[2]

        # Levenshtein and Exact
        lev1 = levenshtein_similarity(a[0], truth_set[j][i][0])
        lev2 = levenshtein_similarity(a[1], truth_set[j][i][1])
        exact1 = string_similarity(a[0], truth_set[j][i][0])
        exact2 = string_similarity(a[1], truth_set[j][i][1])

        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)

        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        tanimoto_fp = DataStructs.TanimotoSimilarity(fp1, fp2)

        tanimoto1 = tanimoto_aligned_characters(a[0], truth_set[j][i][0])
        tanimoto2 = tanimoto_aligned_characters(a[1], truth_set[j][i][1])
                
        optalignment_log.write(f"Original 1: {smiles1}\nOriginal 2: {smiles2}\n")
        optalignment_log.write(f"Aligned 1:  {aligned1}\nAligned 2:  {aligned2}\n")
        optalignment_log.write(f"Alignment Score: {score_val:.4f}\n")

        optalignment_log.write(f"Levenshtein Similarities: {lev1:.4f}, {lev2:.4f}\n")
        optalignment_log.write(f"Exact Match Similarities: {exact1:.4f}, {exact2:.4f}\n")
        optalignment_log.write(
            f"RDKit Fingerprint Tanimoto: "
            f"{round(tanimoto1, 4) if tanimoto1 is not None else 'Invalid'}, "
            f"{round(tanimoto2, 4) if tanimoto2 is not None else 'Invalid'}\n"
        )
        optalignment_log.write("-" * 75 + "\n")

optalignment_log.close()

stop = time.time()

df_sorted = df.sort_values(by='Levenshtein Average', ascending=False)

df_sorted.head(10)

#proof of concept using Pentose phosphate pathway
ppp = ['OC1OC(COP(O)(O)=O)C(O)C(O)C1O',
 'OC1C(COP(O)(O)=O)OC(=O)C(O)C1O',
 'OC(COP(O)(O)=O)C(O)C(O)C(O)C(O)=O',
 'OCC(=O)C(O)C(O)COP(O)(O)=O',
 'OC(C(O)COP(O)(O)=O)C(O)C=O',
 'OCC(=O)C(O)C(O)C(O)C(O)COP(O)(O)=O',
 'OC(C=O)C(O)COP(O)(O)=O',
 'OCC1(O)OC(COP(O)(O)=O)C(O)C1O',
 'OC1OC(COP(O)(O)=O)C(O)C(O)C1O']

data = []

gap_open_penalty = -2
gap_extend_penalty = -2
    
result = []

n = len(ppp)

for i in range(n):
    # Creating a new list starting from the current item and cycling through the list
    cycled_list = ppp[i:] + ppp[:i]

    # Calculate the maximum index for alignment (half the cycle length, rounded up)
    max_index = (len(cycled_list) + 1) // 2

    score_lst = []
    for j in range(max_index):
        # Calculate wrap-around index if needed
        wrap_around_index = len(cycled_list) - j if j != 0 else 0

        # Use direct or wrap-around index based on which is smaller
        align_index = j if j < wrap_around_index else wrap_around_index

        # Call align function with the chosen index
        a = align(cycled_list[0], cycled_list[align_index])
        score_lst.append(a[2])

    result.append(score_lst)

min_length = min(len(lst) for lst in result)

# Initialize a list to hold the sums
sums = [0] * min_length

# Sum up values at each index
for lst in result:
    for i in range(min_length):
        sums[i] += lst[i]

# Calculate the average for each index
averages = [sum_val / len(result) for sum_val in sums]

data.append(averages)

data = np.array(result)

# Calculate mean and standard deviation for each position
means = np.mean(data, axis=0)
std_devs = np.std(data, axis=0)

# Plotting Pentose Phosphate Pathway
plt.errorbar(range(len(means)), means, yerr=std_devs, fmt='o', capsize=5)
plt.plot(range(len(means)), means, label='Mean trend', linestyle='-', marker='o')
plt.title("Best Parameter alignment of Pentose Phosphate Pathway of Paired scoring")
plt.xlabel("Distance from Starting Position")
plt.ylabel("Average Score")
plt.grid(True)

# Set x-ticks to whole numbers
plt.xticks(range(len(means)))

plt.savefig('valid_paired_ppp_bars.png', bbox_inches='tight')
plt.show()

alignment_log.close()
