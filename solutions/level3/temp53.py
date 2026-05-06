from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def level_function(reaction_smiles):
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smiles, useSmiles=True)
        if rxn is None:
            return None

        fp = AllChem.CreateDifferenceFingerprintForReaction(rxn)

        info = {}
        non_zero = fp.GetNonzeroElements()

        return {
            "num_nonzero_elements": len(non_zero),
        }
    except Exception as e:
        print(e)
        return None
