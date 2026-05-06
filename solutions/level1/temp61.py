from rdkit import Chem

def level_function(mol):
    try:
        if not isinstance(mol, str) or len(mol.strip()) == 0:
            return {"valid": False, "error": "Empty or non-string input"}
        mol_obj = Chem.MolFromSmiles(mol, sanitize=False)
        if mol_obj is None:
            return {"valid": False, "error": f"Cannot parse SMILES: '{mol}'"}
        try:
            Chem.SanitizeMol(mol_obj)
        except Exception as san_err:
            return {"valid": False, "error": f"Sanitization failed: {str(san_err)}"}
        return {"valid": True, "canonical_smiles": Chem.MolToSmiles(mol_obj)}
    except Exception as e:
        return {"valid": False, "error": str(e)}
