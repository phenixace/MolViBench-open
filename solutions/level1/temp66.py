from rdkit import Chem
from rdkit.Chem import RDConfig
import os, sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        sa = sascorer.calculateScore(mol_obj)
        return round(sa, 4)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    print(f'Output: {level_function(smiles)}')
