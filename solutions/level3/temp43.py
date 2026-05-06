from rdkit import Chem
import random

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        ring_info = mol_obj.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        if not atom_rings:
            return None

        ring = list(random.choice(list(atom_rings)))

        rw = Chem.RWMol(mol_obj)

        target_idx = random.choice(ring)
        orig_num = rw.GetAtomWithIdx(target_idx).GetAtomicNum()
        candidates = [z for z in replacements if z != orig_num]
        if not candidates:
            return None

        new_num = random.choice(candidates)
        rw.GetAtomWithIdx(target_idx).SetAtomicNum(new_num)

        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            return None
    except Exception as e:
        print(e)
        return None
