from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import random


def level_function(mol, seed=42):

    try:
        random.seed(seed)
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        def lipinski_score(m):

            mw = Descriptors.MolWt(m)
            logp = Crippen.MolLogP(m)
            hbd = rdMolDescriptors.CalcNumHBD(m)
            hba = rdMolDescriptors.CalcNumHBA(m)
            violations = 0
            if mw >= 500: violations += 1
            if logp >= 5: violations += 1
            if hbd > 5: violations += 1
            if hba > 10: violations += 1
            return violations, {"MW": round(mw, 2), "LogP": round(logp, 2),
                               "HBD": hbd, "HBA": hba}

        current_mol = mol_obj
        current_violations, current_props = lipinski_score(current_mol)
        history = [{"round": 0, "smiles": Chem.MolToSmiles(current_mol),
                    "violations": current_violations, **current_props}]

        replacement_atoms = [6, 7, 8, 9]

        for round_num in range(1, 11):
            rw = Chem.RWMol(current_mol)
            num_atoms = rw.GetNumAtoms()
            if num_atoms == 0:
                break


            atom_idx = random.randint(0, num_atoms - 1)
            original_num = rw.GetAtomWithIdx(atom_idx).GetAtomicNum()


            candidates = [a for a in replacement_atoms if a != original_num]
            random.shuffle(candidates)

            mutated = False
            for new_num in candidates:
                try:
                    test_mol = Chem.RWMol(current_mol)
                    test_mol.GetAtomWithIdx(atom_idx).SetAtomicNum(new_num)
                    Chem.SanitizeMol(test_mol)
                    v, props = lipinski_score(test_mol)
                    if v <= current_violations:
                        current_mol = test_mol.GetMol()
                        current_violations = v
                        current_props = props
                        mutated = True
                        break
                except Exception:
                    continue

            history.append({
                "round": round_num,
                "smiles": Chem.MolToSmiles(current_mol),
                "violations": current_violations,
                "mutated": mutated,
                **current_props
            })

        return {
            "final_smiles": Chem.MolToSmiles(current_mol),
            "final_violations": current_violations,
            "history": history
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(CCCCCCC)cc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['final_smiles']}{result['final_violations']}")
