from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """模拟从分子上脱除 Boc 保护基。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Boc deprotection: R-NHC(=O)OC(C)(C)C -> R-NH2
        rxn = AllChem.ReactionFromSmarts(
            '[N:1]C(=O)OC(C)(C)C>>[NH2:1]'
        )

        products = rxn.RunReactants((mol_obj,))
        if products:
            prod = products[0][0]
            try:
                Chem.SanitizeMol(prod)
                return Chem.MolToSmiles(prod)
            except Exception:
                pass

        return None  # No Boc group found
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(C)(C)OC(=O)Nc1ccccc1"
    print(f"Boc脱除: {level_function(smiles)}")
