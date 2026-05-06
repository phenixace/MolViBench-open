from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol1, mol2, consider_stereochemistry=True):
    try:
        mol1 = Chem.MolFromSmiles(mol1)
        mol2 = Chem.MolFromSmiles(mol2)
        if mol1 is None or mol2 is None:
            return None

        formula1 = rdMolDescriptors.CalcMolFormula(mol1)
        formula2 = rdMolDescriptors.CalcMolFormula(mol2)
        if formula1 != formula2:

        smiles1 = Chem.MolToSmiles(mol1, isomericSmiles=consider_stereochemistry, canonical=True)
        smiles2 = Chem.MolToSmiles(mol2, isomericSmiles=consider_stereochemistry, canonical=True)

        if smiles1 != smiles2:
        else:

    except Exception as e:
        print(e)
        return None
