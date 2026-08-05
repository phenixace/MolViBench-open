from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdMolDescriptors


def _synthetic_complexity_score(mol):

    num_rings = mol.GetRingInfo().NumRings()
    num_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    num_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in [1, 6])
    mw = Descriptors.MolWt(mol)
    score = num_rings * 1.5 + num_chiral * 2.0 + num_hetero * 0.5 + mw / 200.0
    return round(score, 2)


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        orig_score = _synthetic_complexity_score(mol_obj)
        orig_smi = Chem.MolToSmiles(mol_obj)

        derivatives = []


        for idx in range(mol_obj.GetNumAtoms()):
            atom = mol_obj.GetAtomWithIdx(idx)
            if atom.GetDegree() == 1 and atom.GetAtomicNum() != 1:
                try:
                    rw = Chem.RWMol(mol_obj)
                    rw.RemoveAtom(idx)
                    Chem.SanitizeMol(rw)
                    smi = Chem.MolToSmiles(rw)
                    if smi and smi != orig_smi:
                        new_mol = Chem.MolFromSmiles(smi)
                        if new_mol:
                            new_score = _synthetic_complexity_score(new_mol)
                            if new_score < orig_score and smi not in [d['smiles'] for d in derivatives]:
                                derivatives.append({
                                    'smiles': smi,
                                    'complexity_score': new_score,
                                    'reduction': round(orig_score - new_score, 2)
                                })
                except Exception:
                    continue


        for idx in range(mol_obj.GetNumAtoms()):
            atom = mol_obj.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() not in [1, 6] and not atom.GetIsAromatic():
                try:
                    rw = Chem.RWMol(mol_obj)
                    rw.GetAtomWithIdx(idx).SetAtomicNum(6)
                    Chem.SanitizeMol(rw)
                    smi = Chem.MolToSmiles(rw)
                    if smi and smi != orig_smi:
                        new_mol = Chem.MolFromSmiles(smi)
                        if new_mol:
                            new_score = _synthetic_complexity_score(new_mol)
                            if new_score < orig_score and smi not in [d['smiles'] for d in derivatives]:
                                derivatives.append({
                                    'smiles': smi,
                                    'complexity_score': new_score,
                                    'reduction': round(orig_score - new_score, 2)
                                })
                except Exception:
                    continue

        derivatives.sort(key=lambda x: x['complexity_score'])
        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'
    result = level_function(smiles)
    print(f'Output: {result}')
