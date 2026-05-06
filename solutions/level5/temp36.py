from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def level_function(mols):
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            mol_data.append({'smiles': Chem.MolToSmiles(mol), 'fp': fp})

        if len(mol_data) <= 10:
            return [d['smiles'] for d in mol_data]

        remaining = list(range(1, len(mol_data)))

        while len(selected) < 10 and remaining:
            best_idx = None
            best_min_dist = -1
            for idx in remaining:
                min_dist = min(
                    1 - DataStructs.TanimotoSimilarity(mol_data[idx]['fp'], mol_data[s]['fp'])
                    for s in selected
                )
                if min_dist > best_min_dist:
                    best_min_dist = min_dist
                    best_idx = idx
            if best_idx is not None:
                selected.append(best_idx)
                remaining.remove(best_idx)

        return [mol_data[i]['smiles'] for i in selected]
    except Exception as e:
        print(e)
        return None
