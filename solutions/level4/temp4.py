from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):
     pass
    try:
        mol = Chem.MolFromSmiles(mol_smi)
        if mol is None: return None

        halogen_pattern = Chem.MolFromSmarts('[F,Cl,Br,I]')
        if not mol.HasSubstructMatch(halogen_pattern):
            return None

        product = Chem.DeleteSubstructs(mol, halogen_pattern)

        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)

        logp = rdMolDescriptors.CalcCrippenDescriptors(product)[0]

        return {
            "has_halogen": True,
            "original_smiles": mol_smi,
            "product": product_smiles,
            "logp": round(logp, 4)
        }
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == "__main__":
    smiles = "c1ccc(Cl)cc1"
    result = level_function(smiles)
    print(f"Output: {result}")
