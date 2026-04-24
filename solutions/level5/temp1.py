from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolDescriptors
import random


def level_function(mols):
    """给定一组已知活性分子，生成与其相似度 >0.7 的新分子候选。"""
    try:
        # Parse active molecules
        active_mols = []
        active_fps = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                active_mols.append(mol)
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                active_fps.append(fp)

        if not active_mols:
            return []

        # Define common substituents for modification
        substituents = ['C', 'O', 'N', 'F', 'Cl', 'CC', 'OC', 'NC']
        candidates = set()

        for mol in active_mols:
            smi = Chem.MolToSmiles(mol)
            rw_mol = Chem.RWMol(mol)

            # Strategy 1: Add atoms to each heavy atom
            for atom_idx in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetImplicitValence() > 0:
                    for sub_smi in ['C', 'N', 'O', 'F', 'Cl']:
                        try:
                            ed = Chem.RWMol(mol)
                            new_idx = ed.AddAtom(Chem.Atom(Chem.MolFromSmiles(sub_smi).GetAtomWithIdx(0).GetAtomicNum()))
                            ed.AddBond(atom_idx, new_idx, Chem.BondType.SINGLE)
                            try:
                                Chem.SanitizeMol(ed)
                                new_smi = Chem.MolToSmiles(ed)
                                if new_smi and new_smi != smi:
                                    candidates.add(new_smi)
                            except Exception:
                                pass
                        except Exception:
                            pass

            # Strategy 2: Replace atoms with different elements
            for atom_idx in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(atom_idx)
                orig_num = atom.GetAtomicNum()
                for new_num in [6, 7, 8]:
                    if new_num != orig_num:
                        try:
                            ed = Chem.RWMol(mol)
                            ed.GetAtomWithIdx(atom_idx).SetAtomicNum(new_num)
                            try:
                                Chem.SanitizeMol(ed)
                                new_smi = Chem.MolToSmiles(ed)
                                if new_smi and new_smi != smi:
                                    candidates.add(new_smi)
                            except Exception:
                                pass
                        except Exception:
                            pass

        # Filter candidates by similarity > 0.7
        results = []
        for cand_smi in candidates:
            cand_mol = Chem.MolFromSmiles(cand_smi)
            if cand_mol is None:
                continue
            cand_fp = AllChem.GetMorganFingerprintAsBitVect(cand_mol, 2, nBits=2048)
            max_sim = max(DataStructs.TanimotoSimilarity(cand_fp, afp) for afp in active_fps)
            if max_sim > 0.7:
                results.append((cand_smi, round(max_sim, 4)))

        results.sort(key=lambda x: x[1], reverse=True)
        return results
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = [
        "c1ccccc1",
        "CC(=O)Oc1ccccc1OC(C)=O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    ]
    result = level_function(smiles_list)
    print(f"result: {result[:10] if result else result}")
