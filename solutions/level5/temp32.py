from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def level_function(mols):
    try:
        mol_data = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            qed = Descriptors.qed(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            mol_data.append({
                'smiles': Chem.MolToSmiles(mol),
                'qed': round(qed, 4),
                'tpsa': round(tpsa, 2)
            })

        for d in mol_data:
            d['score'] = round(d['qed'] - d['tpsa'] / 200, 4)

        mol_data.sort(key=lambda x: x['score'], reverse=True)
        return mol_data[:5]
    except Exception as e:
        print(e)
        return None
