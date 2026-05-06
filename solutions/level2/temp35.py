from rdkit import Chem
import random

def level_function():
    try:
        fragments = [
            "C", "CC", "CCC", "c1ccccc1", "C=O", "C(=O)O", "CCO", "CN",
            "C=C", "C#C", "C(F)(F)F", "c1ccncc1", "C1CCCCC1", "C(=O)N",
            "OC", "NC", "SC", "ClC", "BrC", "FC", "c1ccc(O)cc1",
            "c1ccc(N)cc1", "CC(=O)O", "CCCC", "CCCCC", "c1ccoc1",
            "c1ccsc1", "CC=CC", "C(O)C", "C1CC1"
        ]
        smiles = random.choice(fragments)
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol)
        return None
    except Exception as e:
        print(e)
        return None
