import os
os.environ['RDKIT_LOG_LEVEL'] = 'ERROR'
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import sys
import contextlib
from Bio.KEGG import REST
import re
import pubchempy as pcp
from indigo import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import time

def safe_get_compounds(name, max_retries=5, delay=3):
    for attempt in range(max_retries):
        try:
            return pcp.get_compounds(name, 'name')
        except pcp.PubChemHTTPError as e:
            print(f"PubChem busy (attempt {attempt + 1}), retrying in {delay} seconds...")
            time.sleep(delay)
    print(f"Failed to retrieve: {name}")
    return []

#function that retreives the name of the compounds given a compound id
def get_compound_name(cpd_id):
    entry = REST.kegg_get(cpd_id).read()  # Fetch the compound entry from KEGG
    for line in entry.rstrip().split("\n"):  # Parse the entry line by line
        if line.startswith("NAME"):  # If the line contains the compound name
            name = line.split(" ", 1)[1]  # Extract the name
            return name.strip().rstrip(';')  # Remove leading/trailing whitespace and trailing semicolon
        
#Retreives all the canonical SMILES in a given pathway id, and returns them as list
def pathway(kegg_map):
    result = REST.kegg_link('compound',kegg_map).read()
    
    lines = result.strip().split('\n')  
    cpd_codes = [line.split('\t')[1] for line in lines]
    
    names = []
    for codes in cpd_codes:
        names.append(get_compound_name(codes))
        
    smiles = []

    # Fetch the SMILES strings
    for compound in names:
        results = safe_get_compounds(compound)
        if results:
            smiles.append(results[0].canonical_smiles)
            
    return(smiles)

#pathway for citrate cycle
pathway('M00009')

#SMILES list with co-enzymes removed
smiles_citrate = ('C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O',
 'C(C(C(C(=O)O)O)C(=O)O)C(=O)O',
 'C(CC(=O)O)C(=O)C(=O)O',
 'C(CC(=O)[O-])C(=O)[O-]',
 'C(=CC(=O)[O-])C(=O)[O-]',
 'C(C(C(=O)O)O)C(=O)O','C(C(=O)C(=O)O)C(=O)O')

#Re-canonization using indigo
canonized_citrate = []
indigo = Indigo()
for smile in smiles_citrate:
    mol1 = indigo.loadMolecule(smile)
    mol1.aromatize()
    canonized_citrate.append(mol1.canonicalSmiles())

canonized_citrate = ('[O-]C(=O)C(=O)CC([O-])=O',
    '[O-]C(=O)C(O)(CC([O-])=O)CC([O-])=O',
    '[O-]C(=O)C(CC([O-])=O)C(O)C([O-])=O',
    '[O-]C(=O)CCC(=O)C([O-])=O',
    '[O-]C(=O)CCC([O-])=O',
    '[O-]C(=O)C=CC([O-])=O',
    '[O-]C(=O)C(O)CC([O-])=O')

#Pathway for Pentose Phosphate Pathway
pathway('M00004')

paths = ['C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O',
 'C(C1C(C(C(C(=O)O1)O)O)O)OP(=O)(O)O',
 'C(C(C(C(C(C(=O)O)O)O)O)O)OP(=O)(O)O',
 'C(C(C(C(=O)CO)O)O)OP(=O)(O)O',
 'C(C(C(C(C=O)O)O)O)OP(=O)(O)O',         
 'C(C(C(C(C(C(=O)CO)O)O)O)O)OP(=O)(O)O',
 'C(C(C(C=O)O)O)OP(=O)(O)O',        
 'C(C1C(C(C(O1)(CO)O)O)O)OP(=O)(O)O',
 'C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O']

#Re-canonization using indigo
canonized_ppp = []
indigo = Indigo()
for smile in paths:
    mol1 = indigo.loadMolecule(smile)
    mol1.aromatize()
    canonized_ppp.append(mol1.canonicalSmiles())

print(canonized_ppp)


#pathway for Glycolysis (linear)
pathway('M00001')

glyc = ('C(C1C(C(C(C(O1)O)O)O)O)O',
 'C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O',
 'C(C(C(C(C(=O)CO)O)O)O)OP(=O)(O)O',
 'C(C(=O)COP(=O)(O)O)O',
 'C(C(C=O)O)OP(=O)(O)O',
 'C(C(C(=O)OP(=O)(O)O)O)OP(=O)(O)O',
 'C(C(C(=O)O)O)OP(=O)(O)O',
 'C(C(C(=O)O)OP(=O)(O)O)O',
 'C=C(C(=O)O)OP(=O)(O)O',
 'CC(=O)C(=O)[O-]')

#Re-canonization using indigo
canonized_glycolysis = []
indigo = Indigo()
for smile in glyc:
    mol1 = indigo.loadMolecule(smile)
    mol1.aromatize()
    canonized_glycolysis.append(mol1.canonicalSmiles())

print(canonized_glycolysis)

from Levenshtein import distance as levenshtein_distance

def generate_variants_with_substitution(parent_smiles, original_atom='O', new_atom='S', max_variants=10, visualize=False):
    variants = []

    rdmol = Chem.MolFromSmiles(parent_smiles)
    if rdmol is None:
        return variants  

    count = 0
    for atom in rdmol.GetAtoms():
        if atom.GetSymbol() == original_atom:
            editable = Chem.RWMol(rdmol)  
            editable.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(new_atom))
            variant_smiles = Chem.MolToSmiles(editable, canonical=False)

            if visualize:
                img = Draw.MolToImage(editable, size=(300, 300))
                img.show()

            try:
                mol = indigo.loadMolecule(variant_smiles)
                mol.aromatize()
                variants.append(mol.canonicalSmiles())
                count += 1
            except:
                continue

            if count >= max_variants:
                break

    return variants

from statistics import mean
from statistics import stdev

print("Levenshtein analysis on Krebs Cycle:")

def tanimoto_distance(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return None  

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)

    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return 1 - similarity  

for parent in smiles_citrate:
    # Canonize the parent using Indigo
    mol = indigo.loadMolecule(parent)
    mol.aromatize()
    canon_parent = mol.canonicalSmiles()
    variants = generate_variants_with_substitution(parent, new_atom='S')  

    # Compare variants to the canonicalized parent
    distances = [levenshtein_distance(canon_parent, v) for v in variants]
    print(f"\nParent: {canon_parent}")
    print("Variants:", variants)
    print("Levenshtein distances:", distances)
    if distances:
        print("Average Levenshtein distance:", mean(distances))
        print("Standard deviation:", stdev(distances) if len(distances) > 1 else 0)
    else:
        print("No distances to analyze.")

    tanimoto_distances = []
    for v in variants:
        d = tanimoto_distance(canon_parent, v)
        if d is not None:
            tanimoto_distances.append(d)

    if tanimoto_distances:
        print("Tanimoto distances:", tanimoto_distances)
        print("Average Tanimoto distance:", mean(tanimoto_distances))
        print("Standard deviation (Tanimoto):", stdev(tanimoto_distances) if len(tanimoto_distances) > 1 else 0)
    else:
        print("No valid Tanimoto distances.")
