from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_smi = Chem.MolToSmiles(mol_obj)


        ring_info = mol_obj.GetRingInfo()
        atom_rings = list(ring_info.AtomRings())

        derivatives = []


        stable_replacements = [(7, 'N'), (8, 'O'), (16, 'S')]
        for ring in atom_rings:
            for idx in ring:
                atom = mol_obj.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    for new_num, sym in stable_replacements:
                        try:
                            rw = Chem.RWMol(mol_obj)
                            rw.GetAtomWithIdx(idx).SetAtomicNum(new_num)
                            Chem.SanitizeMol(rw)
                            smi = Chem.MolToSmiles(rw)
                            if smi != orig_smi and smi not in [d['smiles'] for d in derivatives]:
                                qed = Descriptors.qed(rw)
                                derivatives.append({
                                    'smiles': smi,
                                    'qed': round(qed, 4),
                                    'modification': f'Replace ring C with {sym}'
                                })
                        except Exception:
                            continue

        derivatives.sort(key=lambda x: x['qed'], reverse=True)
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    result = level_function(smiles)
    print(f'Output: {result}')
