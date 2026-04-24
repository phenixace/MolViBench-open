from rdkit import Chem
from rdkit.Chem import Descriptors
import random


def level_function(mol):
    """给定一个候选分子，生成更小的衍生物以降低分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_mw = Descriptors.MolWt(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        derivatives = []

        # 策略: 删除末端原子/侧链
        terminal_atoms = [
            atom.GetIdx() for atom in mol_obj.GetAtoms()
            if atom.GetDegree() == 1 and atom.GetAtomicNum() != 1
        ]

        for idx in terminal_atoms:
            try:
                rw = Chem.RWMol(mol_obj)
                rw.RemoveAtom(idx)
                Chem.SanitizeMol(rw)
                smi = Chem.MolToSmiles(rw)
                if smi and smi != orig_smi and smi not in [d['smiles'] for d in derivatives]:
                    new_mol = Chem.MolFromSmiles(smi)
                    if new_mol:
                        new_mw = Descriptors.MolWt(new_mol)
                        if new_mw < orig_mw:
                            derivatives.append({
                                'smiles': smi,
                                'mw': round(new_mw, 2),
                                'mw_reduction': round(orig_mw - new_mw, 2)
                            })
            except Exception:
                continue

        derivatives.sort(key=lambda x: x['mw'])
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # 布洛芬
    result = level_function(smiles)
    print(f"更小的衍生物: {result}")
