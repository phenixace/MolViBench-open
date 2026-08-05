from rdkit import Chem
from rdkit.Chem import Descriptors


def _synthetic_complexity_score(mol):

    num_rings = mol.GetRingInfo().NumRings()
    num_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    num_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in [1, 6])
    mw = Descriptors.MolWt(mol)
    return round(num_rings * 1.5 + num_chiral * 2.0 + num_hetero * 0.5 + mw / 200.0, 2)


def level_function(mols):

    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            score = _synthetic_complexity_score(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'complexity_score': score
            })

        mol_data.sort(key=lambda x: x['complexity_score'])
        return mol_data[:5]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles_list = ['CCO', 'c1ccccc1', 'CC(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'c1ccncc1', 'CC[C@H](F)C(=O)O']
    result = level_function(smiles_list)
    if result:
        for r in result:
            print(f"Output: {r['smiles']}{r['complexity_score']}")
