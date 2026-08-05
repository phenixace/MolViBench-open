from rdkit import Chem
import random

def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None: return None

        rw = Chem.RWMol(mol_obj)
        num_atoms = rw.GetNumAtoms()
        if num_atoms == 0: return None


        replacements = [6, 7, 8, 9, 16, 17]


        indices = list(range(num_atoms))
        random.shuffle(indices)

        for idx in indices:
            atom = rw.GetAtomWithIdx(idx)
            orig_num = atom.GetAtomicNum()


            candidates = [z for z in replacements if z != orig_num]
            if not candidates: continue

            new_num = random.choice(candidates)


            atom.SetAtomicNum(new_num)


            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(False)

            try:

                atom.UpdatePropertyCache()

                Chem.SanitizeMol(rw, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                return Chem.MolToSmiles(rw)
            except:

                atom.SetAtomicNum(orig_num)
                continue

        return None
    except Exception as e:
        return None

if __name__ == '__main__':
    smiles = 'c1ccccc1'
    for i in range(10):
        result = level_function(smiles)
        print(f'Output: {i + 1}{result}')
