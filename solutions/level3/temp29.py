from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):
    """模拟消除反应。"""
    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None
        elimination_patterns = [
            '[C:1]([H])[C:2][F:3]>>[C:1]=[C:2]',
            '[C:1]([H])[C:2][Cl:3]>>[C:1]=[C:2]',
            '[C:1]([H])[C:2][Br:3]>>[C:1]=[C:2]',
            '[C:1]([H])[C:2][I:3]>>[C:1]=[C:2]',
        ]
        result_smiles = []
        for pattern in elimination_patterns:
            rxn = AllChem.ReactionFromSmarts(pattern)
            products = rxn.RunReactants((molecule,))
            if products:
                for product_set in products:
                    for product in product_set:
                        try:
                            Chem.SanitizeMol(product)
                            smi = Chem.MolToSmiles(product)
                            if smi not in result_smiles:
                                result_smiles.append(smi)
                        except Exception:
                            continue
        return result_smiles if result_smiles else None
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CCCBr"
    print(f"消除反应产物: {level_function(smiles)}")
