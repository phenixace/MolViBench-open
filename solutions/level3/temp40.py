from rdkit import Chem
import random

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        rw = Chem.RWMol(mol_obj)
        num_atoms = rw.GetNumAtoms()
        if num_atoms == 0:
            return None

        idx = random.randint(0, num_atoms - 1)
        orig_num = rw.GetAtomWithIdx(idx).GetAtomicNum()
        candidates = [z for z in replacements if z != orig_num]
        if not candidates:
            return None

        new_num = random.choice(candidates)
        rw.GetAtomWithIdx(idx).SetAtomicNum(new_num)

        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            return None
    except Exception as e:
        print(e)
        return None
