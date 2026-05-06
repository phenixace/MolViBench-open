from rdkit import Chem

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        pattern = Chem.MolFromSmarts("n")
        has_aromatic_n = mol_obj.HasSubstructMatch(pattern)

        if not has_aromatic_n:
            return None

        rw = Chem.RWMol(mol_obj)
        for atom in rw.GetAtoms():
            if atom.GetIsAromatic() and atom.GetAtomicNum() == 7:
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)

        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        product_smiles = Chem.MolToSmiles(rw)

        total_charge = sum(atom.GetFormalCharge() for atom in rw.GetAtoms())

        return {
            "has_aromatic_n": has_aromatic_n,
            "product": product_smiles,
            "total_charge": total_charge
        }
    except Exception as e:
        print(e)
        return None
