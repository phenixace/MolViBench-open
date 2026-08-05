from rdkit import Chem
from rdkit.Chem import Descriptors


def level_function(mols):

    try:
        results = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            logp = Descriptors.MolLogP(mol)
            if 0 <= logp <= 5:
                results.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "logp": round(logp, 2)
                })
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1', 'CCCCCCCCCCCCCCCC', 'O']
    print(f'Output: {level_function(smiles_list)}')
