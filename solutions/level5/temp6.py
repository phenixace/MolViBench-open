from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def level_function(mol):
    """给定一个蛋白结合口袋中的片段（fragment），生成 fragment-growing 的衍生物。"""
    try:
        frag = Chem.MolFromSmiles(mol)
        if frag is None:
            return None

        original_smi = Chem.MolToSmiles(frag)
        derivatives = set()

        # Growth groups to add
        growth_groups = [
            ('methyl', 6, None),
            ('amino', 7, None),
            ('hydroxyl', 8, None),
            ('fluoro', 9, None),
            ('chloro', 17, None),
        ]

        # Reaction-based growth for more complex groups
        growth_reactions = [
            ('ethyl', '[*:1]([H])>>[*:1]CC'),
            ('methoxy', '[cH1:1]>>[c:1]OC'),
            ('acetyl', '[cH1:1]>>[c:1]C(C)=O'),
            ('carboxyl', '[cH1:1]>>[c:1]C(=O)O'),
            ('amide', '[cH1:1]>>[c:1]C(=O)N'),
            ('sulfonamide', '[cH1:1]>>[c:1]S(=O)(=O)N'),
            ('cyano', '[cH1:1]>>[c:1]C#N'),
            ('phenyl', '[cH1:1]>>[c:1]-c1ccccc1'),
            ('morpholine', '[cH1:1]>>[c:1]N1CCOCC1'),
            ('piperidine', '[cH1:1]>>[c:1]N1CCCCC1'),
        ]

        # Simple atom additions at available positions
        for atom_idx in range(frag.GetNumAtoms()):
            atom = frag.GetAtomWithIdx(atom_idx)
            if atom.GetNumImplicitHs() > 0:
                for name, atomic_num, _ in growth_groups:
                    try:
                        rw_mol = Chem.RWMol(frag)
                        new_idx = rw_mol.AddAtom(Chem.Atom(atomic_num))
                        rw_mol.AddBond(atom_idx, new_idx, Chem.BondType.SINGLE)
                        Chem.SanitizeMol(rw_mol)
                        new_smi = Chem.MolToSmiles(rw_mol)
                        if new_smi and new_smi != original_smi:
                            derivatives.add(new_smi)
                    except Exception:
                        pass

        # Reaction-based growth
        for name, rxn_smarts in growth_reactions:
            try:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((frag,))
                for prod_set in products:
                    for prod in prod_set:
                        try:
                            Chem.SanitizeMol(prod)
                            new_smi = Chem.MolToSmiles(prod)
                            if new_smi and new_smi != original_smi:
                                derivatives.add(new_smi)
                        except Exception:
                            pass
            except Exception:
                pass

        # Return derivatives with MW
        result = []
        for smi in derivatives:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                mw = round(Descriptors.MolWt(m), 2)
                result.append((smi, mw))

        result.sort(key=lambda x: x[1])
        return result
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccncc1")
    print(f"result: {result}")
