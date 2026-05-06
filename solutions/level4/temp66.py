from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

def level_function(fragment_smiles_list, num_rounds=3):
    try:
        frags = []
        for smi in fragment_smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                frags.append(mol)

        if not frags:
            return None

        current_mol = frags[0]
        mw_curve = [{"round": 0, "MW": round(Descriptors.MolWt(current_mol), 2),
                     "smiles": Chem.MolToSmiles(current_mol)}]

        for round_num in range(1, num_rounds + 1):
            frag_idx = round_num % len(frags)
            frag = frags[frag_idx]

            combined = Chem.CombineMols(current_mol, frag)
            rw = Chem.RWMol(combined)

            n1 = current_mol.GetNumAtoms()
            attached = False
            for i in range(n1):
                atom1 = rw.GetAtomWithIdx(i)
                if atom1.GetImplicitValence() > 0:
                    for j in range(n1, rw.GetNumAtoms()):
                        atom2 = rw.GetAtomWithIdx(j)
                        if atom2.GetImplicitValence() > 0:
                            try:
                                rw.AddBond(i, j, Chem.BondType.SINGLE)
                                Chem.SanitizeMol(rw)
                                current_mol = rw.GetMol()
                                attached = True
                                break
                            except Exception:
                                rw = Chem.RWMol(combined)
                                continue
                    if attached:
                        break

            if not attached:
                current_mol = combined

            mw = Descriptors.MolWt(current_mol)
            mw_curve.append({
                "round": round_num,
                "MW": round(mw, 2),
                "smiles": Chem.MolToSmiles(current_mol),
                "fragment_added": Chem.MolToSmiles(frag)
            })

        return {
            "num_rounds": num_rounds,
            "mw_growth_curve": mw_curve,
            "final_MW": mw_curve[-1]["MW"],
            "MW_increase": round(mw_curve[-1]["MW"] - mw_curve[0]["MW"], 2)
        }
    except Exception as e:
        print(e)
        return None
