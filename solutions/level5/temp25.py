from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):
    """给定一个候选分子，生成更稳定的环结构衍生物。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_smi = Chem.MolToSmiles(mol_obj)

        # 策略: 1) 用杂原子替换环中原子  2) 芳环取代
        ring_info = mol_obj.GetRingInfo()
        atom_rings = list(ring_info.AtomRings())

        derivatives = []

        # 策略 1: 替换环上原子为更稳定的杂原子
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
                                    'modification': f'环上C→{sym}'
                                })
                        except Exception:
                            continue

        derivatives.sort(key=lambda x: x['qed'], reverse=True)
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccccc1"
    result = level_function(smiles)
    print(f"环结构衍生物: {result}")
