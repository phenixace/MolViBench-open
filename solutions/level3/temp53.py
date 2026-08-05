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
            "nonzero_elements": dict(list(non_zero.items())[:20])
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    rxn_smi = 'CCO.CC(=O)O>>CC(=O)OCC.O'
    result = level_function(rxn_smi)
    if result:
        print(f"Output: {result['num_nonzero_elements']}")
