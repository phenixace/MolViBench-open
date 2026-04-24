from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


def level_function(mol):
    """给定分子 → 判断是否含酯基 → 若有 → 断裂为两部分 → 计算片段分子量。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        # Step 1: 判断是否含酯基
        pattern = Chem.MolFromSmarts('[C:1](=O)[O:2][C:3]')
        has_ester = mol_obj.HasSubstructMatch(pattern)

        if not has_ester:
            return None

        # Step 2: 断裂酯基 (R-COO-R' → R-COOH + R'-OH)
        rxn = AllChem.ReactionFromSmarts(
            '[C:1](=O)[O:2][C:3]>>[C:1](=O)O.[C:3]O'
        )
        products = rxn.RunReactants((mol_obj,))
        if not products:
            return None

        # Step 3: 计算每个片段的分子量
        fragments = []
        for product in products[0]:
            Chem.SanitizeMol(product)
            smi = Chem.MolToSmiles(product)
            mw = rdMolDescriptors.CalcExactMolWt(product)
            fragments.append({
                "smiles": smi,
                "molecular_weight": round(mw, 4)
            })

        return {
            "has_ester": has_ester,
            "fragments": fragments
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "CC(=O)OCC"
    print(f"result: {level_function(smiles)}")
