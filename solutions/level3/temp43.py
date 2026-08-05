from rdkit import Chem
import random

def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None: return None

        ring_info = mol_obj.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        if not atom_rings: return None


        replacements = [6, 7, 8, 16]


        for _ in range(10):
            rw = Chem.RWMol(mol_obj)

            ring = list(random.choice(atom_rings))

            target_idx = random.choice(ring)
            atom = rw.GetAtomWithIdx(target_idx)

            orig_num = atom.GetAtomicNum()
            candidates = [z for z in replacements if z != orig_num]
            if not candidates: continue

            new_num = random.choice(candidates)
            atom.SetAtomicNum(new_num)


            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)


            atom.SetIsAromatic(False)

            try:

                rw.UpdatePropertyCache(strict=False)

                Chem.SanitizeMol(rw, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                return Chem.MolToSmiles(rw)
            except:
                continue

        return None
    except Exception as e:
        return None

if __name__ == '__main__':
    smiles = 'c1ccccc1'
    for i in range(10):
        res = level_function(smiles)
        print(f'Output: {i + 1}{res}')
