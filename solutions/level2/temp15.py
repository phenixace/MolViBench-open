from rdkit import Chem
from rdkit.Chem import AllChem
import random

def level_function(mol):



    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mol_h = Chem.AddHs(mol_obj)

        h_indices = [atom.GetIdx() for atom in mol_h.GetAtoms() if atom.GetAtomicNum() == 1]
        if not h_indices:
            return None

        h_idx = random.choice(h_indices)

        rw_mol = Chem.RWMol(mol_h)
        rw_mol.GetAtomWithIdx(h_idx).SetAtomicNum(9)
        Chem.SanitizeMol(rw_mol)

        rw_mol = Chem.RemoveHs(rw_mol)
        return Chem.MolToSmiles(rw_mol)
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'c1ccccc1'
    result = level_function(smiles)
    print(f'Output: {result}')
