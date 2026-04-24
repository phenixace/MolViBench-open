from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def level_function(reaction_smiles):
    """计算一个化学反应的反应差异指纹（reaction difference fingerprint）。"""
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if rxn is None:
            return None

        # Generate reaction fingerprint (difference fingerprint)
        fp = AllChem.CreateDifferenceFingerprintForReaction(rxn)

        # Convert to explicit bit representation
        info = {}
        non_zero = fp.GetNonzeroElements()

        return {
            "num_nonzero_elements": len(non_zero),
            "nonzero_elements": dict(list(non_zero.items())[:20])  # first 20
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    rxn_smi = "CCO.CC(=O)O>>CC(=O)OCC.O"
    result = level_function(rxn_smi)
    if result:
        print(f"非零元素数: {result['num_nonzero_elements']}")
