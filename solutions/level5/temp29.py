from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_smi = Chem.MolToSmiles(mol_obj)
        orig_fp = AllChem.GetMorganFingerprintAsBitVect(mol_obj, 2, nBits=2048)

        ring_info = mol_obj.GetRingInfo()
        atom_rings = list(ring_info.AtomRings())

        if not atom_rings:
            return None


        heterocyclic_atoms = set()
        for ring in atom_rings:
            for idx in ring:
                if mol_obj.GetAtomWithIdx(idx).GetAtomicNum() != 6:
                    heterocyclic_atoms.update(ring)
                    break

        if not heterocyclic_atoms:
            return None



        replacements = [6, 7, 8, 16]
        results = []

        for idx in heterocyclic_atoms:
            atom = mol_obj.GetAtomWithIdx(idx)
            orig_num = atom.GetAtomicNum()
            if orig_num == 6:
                continue
            for new_num in replacements:
                if new_num == orig_num:
                    continue
                try:
                    rw = Chem.RWMol(mol_obj)

                    Chem.Kekulize(rw, clearAromaticFlags=True)
                    target_atom = rw.GetAtomWithIdx(idx)
                    target_atom.SetAtomicNum(new_num)

                    target_atom.SetNoImplicit(False)
                    target_atom.SetNumExplicitHs(0)

                    Chem.SanitizeMol(rw)
                    smi = Chem.MolToSmiles(rw)
                    if smi != orig_smi and smi not in [r['smiles'] for r in results]:
                        new_fp = AllChem.GetMorganFingerprintAsBitVect(rw, 2, nBits=2048)
                        sim = DataStructs.TanimotoSimilarity(orig_fp, new_fp)
                        if sim > 0.6:
                            results.append({
                                'smiles': smi,
                                'similarity': round(sim, 4)
                            })
                except Exception:
                    continue

        results.sort(key=lambda x: x['similarity'], reverse=True)
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccncc1CC(=O)O'
    result = level_function(smiles)
    print(f'Output: {result}')
