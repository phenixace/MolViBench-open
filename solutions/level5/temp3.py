from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """给定一个 scaffold（核心骨架），生成 10 个不同侧链修饰的分子。"""
    try:
        scaffold = Chem.MolFromSmiles(mol)
        if scaffold is None:
            return None

        # 10 different side chains
        side_chains = {
            'methyl': 6,        # C -> CH3
            'ethyl': None,      # handled via SMILES
            'propyl': None,
            'hydroxyl': 8,      # O -> OH
            'amino': 7,         # N -> NH2
            'fluoro': 9,        # F
            'chloro': 17,       # Cl
            'methoxy': None,
            'carboxyl': None,
            'cyano': None,
        }

        # Simple atom additions
        simple_atoms = {
            'methyl': ('C', 6),
            'hydroxyl': ('O', 8),
            'amino': ('N', 7),
            'fluoro': ('F', 9),
            'chloro': ('Cl', 17),
        }

        # Complex group additions via reaction SMARTS
        complex_groups = {
            'ethyl': '[cH1:1]>>[c:1]CC',
            'propyl': '[cH1:1]>>[c:1]CCC',
            'methoxy': '[cH1:1]>>[c:1]OC',
            'carboxyl': '[cH1:1]>>[c:1]C(=O)O',
            'cyano': '[cH1:1]>>[c:1]C#N',
        }

        original_smi = Chem.MolToSmiles(scaffold)
        derivatives = []

        # Find first atom with implicit H for simple additions
        target_idx = None
        for atom_idx in range(scaffold.GetNumAtoms()):
            atom = scaffold.GetAtomWithIdx(atom_idx)
            if atom.GetNumImplicitHs() > 0:
                target_idx = atom_idx
                break

        # Simple atom additions
        if target_idx is not None:
            for name, (symbol, atomic_num) in simple_atoms.items():
                try:
                    rw_mol = Chem.RWMol(scaffold)
                    new_idx = rw_mol.AddAtom(Chem.Atom(atomic_num))
                    rw_mol.AddBond(target_idx, new_idx, Chem.BondType.SINGLE)
                    Chem.SanitizeMol(rw_mol)
                    new_smi = Chem.MolToSmiles(rw_mol)
                    if new_smi and new_smi != original_smi:
                        derivatives.append((name, new_smi))
                except Exception:
                    pass

        # Complex group additions via reactions
        for name, rxn_smarts in complex_groups.items():
            try:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((scaffold,))
                if products:
                    prod = products[0][0]
                    Chem.SanitizeMol(prod)
                    new_smi = Chem.MolToSmiles(prod)
                    if new_smi and new_smi != original_smi:
                        derivatives.append((name, new_smi))
            except Exception:
                pass

        # If we have fewer than 10, try adding to other positions
        if len(derivatives) < 10:
            for atom_idx in range(scaffold.GetNumAtoms()):
                if len(derivatives) >= 10:
                    break
                atom = scaffold.GetAtomWithIdx(atom_idx)
                if atom.GetNumImplicitHs() > 0 and atom_idx != target_idx:
                    try:
                        rw_mol = Chem.RWMol(scaffold)
                        new_idx = rw_mol.AddAtom(Chem.Atom(6))
                        rw_mol.AddBond(atom_idx, new_idx, Chem.BondType.SINGLE)
                        Chem.SanitizeMol(rw_mol)
                        new_smi = Chem.MolToSmiles(rw_mol)
                        if new_smi and new_smi != original_smi:
                            derivatives.append((f'methyl_pos{atom_idx}', new_smi))
                    except Exception:
                        pass

        return derivatives[:10]
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    result = level_function("c1ccccc1")
    print(f"result: {result}")
