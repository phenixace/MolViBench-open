from rdkit import Chem
from rdkit.Chem import SaltRemover


def level_function(mol):
    """从含盐 SMILES（如 [Na+].[Cl-].CCO）中提取最大有机片段。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        remover = SaltRemover.SaltRemover()
        stripped = remover.StripMol(mol_obj)
        if stripped.GetNumAtoms() == 0:
            # If all fragments removed, return largest fragment
            frags = Chem.GetMolFrags(mol_obj, asMols=True, sanitizeFrags=True)
            if not frags:
                return None
            largest = max(frags, key=lambda m: m.GetNumAtoms())
            return Chem.MolToSmiles(largest)
        return Chem.MolToSmiles(stripped)
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "[Na+].[Cl-].CCO"
    print(f"最大有机片段: {level_function(smiles)}")
