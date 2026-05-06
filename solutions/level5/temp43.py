from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        ring_info = mol_obj.GetRingInfo()
        atom_rings = list(ring_info.AtomRings())

        if not atom_rings:
            return None

        isomers = set()
        isomers.add(Chem.MolToSmiles(mol_obj))

        for ring in atom_rings:
            for idx in ring:
                orig_num = mol_obj.GetAtomWithIdx(idx).GetAtomicNum()
                for new_num in replacements:
                    if new_num == orig_num:
                        continue
                    try:
                        rw = Chem.RWMol(mol_obj)
                        rw.GetAtomWithIdx(idx).SetAtomicNum(new_num)
                        Chem.SanitizeMol(rw)
                        smi = Chem.MolToSmiles(rw)
                        isomers.add(smi)
                    except Exception:
                        continue

        best = None
        for smi in isomers:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            tpsa = rdMolDescriptors.CalcTPSA(m)
            if best is None or tpsa < best['tpsa']:
                best = {'smiles': smi, 'tpsa': round(tpsa, 2)}

        if best is None:
            return None

        best_mol = Chem.MolFromSmiles(best['smiles'])
        best['qed'] = round(Descriptors.qed(best_mol), 4)

        return {
            'num_isomers': len(isomers),
            'best': best
        }
    except Exception as e:
        print(e)
        return None
