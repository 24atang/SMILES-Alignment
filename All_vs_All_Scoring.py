from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np    
import math
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

smiles = []

with open(r'smiles.txt', 'r') as fp:
    for line in fp:
        x = line[:-1]
        smiles.append(x)

#Returns all gast_charges in database into a list
gast_charges = []
for i in smiles:
    mol = Chem.MolFromSmiles(i)
    AllChem.ComputeGasteigerCharges(mol)
    charges = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms()]
    gast_charges.append(charges)
    
gast_charges = [item for sublist in gast_charges for item in sublist]
gast_charges = [x for x in gast_charges if not math.isnan(x)]
gast_charges = [x for x in gast_charges if not math.isinf(x)]

plt.hist(gast_charges, bins=50, edgecolor='black')
plt.title("Distribution of Gasteiger Charges")
plt.xlabel("Gasteiger Charge")
plt.ylabel("Frequency")
plt.grid(True)
plt.show()

sns.kdeplot(gast_charges, fill=True, color="skyblue", linewidth=1.5)
plt.xlabel("Gasteiger Charge")
plt.ylabel("Density")
plt.title("Density Plot of Gasteiger Charges")
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

sns.histplot(gast_charges, bins=50, kde=True, stat="density", color="skyblue", edgecolor="black")
plt.xlabel("Gasteiger Charge")
plt.ylabel("Density")
plt.title("Distribution of Gasteiger Charges")
sns.despine()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

#function that yields database of differences
def gast_diff(gast):
    gast_charge = []
    for i in range(0,len(gast)):
        for j in range(i+1,len(gast)):
            gast_charge.append(abs(gast[i]-gast[j]))
            
                
    return(gast_charge)

#creates the scoring matrix of all vs all
score = {}
total_count = len(gast_charges)  
for i in range(0,30):
    test = round(i*10**-1,1) #where 10^d, d is number of zeros after max of range
    count = len([x for x in gast_charges if x >= test])    
    prob = count/total_count
    score[test] = np.log2(prob)

with open('score_all_v_all.pkl', 'wb') as f:
    pickle.dump(score, f)
