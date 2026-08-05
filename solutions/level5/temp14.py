from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mols):

    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mw = Descriptors.MolWt(mol)
            if mw < 500:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "mw": round(mw, 2)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1', 'CC(=O)O']
    print(f'Output: {level_function(smiles_list)}')
