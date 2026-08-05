from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(filename):
    try:
        mol = Chem.MolFromMolFile(filename, removeHs=True)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=10)
    Chem.MolToMolFile(mol, "test_input.mol")
    result = level_function("test_input.mol")
    print(f"Output: {result}")
