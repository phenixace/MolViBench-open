from rdkit import Chem
import random


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        rw = Chem.RWMol(mol_obj)
        bonds = list(rw.GetBonds())
        if not bonds:
            return None

        bond = random.choice(bonds)
        bond_idx = bond.GetIdx()
        current_type = bond.GetBondType()


        bond_types = [Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]
        candidates = [bt for bt in bond_types if bt != current_type]
        if not candidates:
            return None

        new_type = random.choice(candidates)
        rw.GetBondWithIdx(bond_idx).SetBondType(new_type)

        try:
            Chem.SanitizeMol(rw)
            return Chem.MolToSmiles(rw)
        except Exception:
            return None
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'CCCC'
    for _ in range(5):
        print(f'Output: {level_function(smiles)}')
