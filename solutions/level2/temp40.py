from rdkit import Chem
from rdkit.Chem import rdchem
import random

def level_function(mol):
    try:
        mol = Chem.MolFromSmiles(mol)
        rwmol = Chem.RWMol(mol)
        bonds = list(rwmol.GetBonds())
        if not bonds:
            return None
        bond = random.choice(bonds)
        bond_types = [
            rdchem.BondType.SINGLE,
            rdchem.BondType.DOUBLE,
            rdchem.BondType.TRIPLE,
        ]
        current_type = bond.GetBondType()
        new_types = [bt for bt in bond_types if bt != current_type]
        new_type = random.choice(new_types)
        bond.SetBondType(new_type)
        try:
            Chem.SanitizeMol(rwmol)
            return Chem.MolToSmiles(rwmol)
        except Exception:
            bond.SetBondType(current_type)
            return Chem.MolToSmiles(mol)
    except Exception as e:
        print(e)
        return None
