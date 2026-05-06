from rdkit import Chem
from rdkit.Chem import Descriptors

def level_function(mol1, mol2):
    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return False

        formula1 = Descriptors.MolecularFormula(m1)
        formula2 = Descriptors.MolecularFormula(m2)
        if formula1 != formula2:
            return False

        smi1 = Chem.MolToSmiles(m1, isomericSmiles=False)
        smi2 = Chem.MolToSmiles(m2, isomericSmiles=False)

        return smi1 != smi2
    except Exception as e:
        print(e)
        return False
