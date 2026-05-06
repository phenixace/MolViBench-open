from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        return (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
    except Exception as e:
        print(e)
        return None
