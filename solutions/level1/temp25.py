from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        return AllChem.CalcNumAtomStereoCenters(mol)
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CC[C@H](F)C(=O)O'
    print(f'Output: {level_function(smiles)}')
