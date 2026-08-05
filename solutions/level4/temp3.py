from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def level_function(mol_smi):
    try:
        mol = Chem.MolFromSmiles(mol_smi)
        if mol is None: return None



        pattern = Chem.MolFromSmarts('[NH2]')
        if not mol.HasSubstructMatch(pattern):
            return None





        rxn = AllChem.ReactionFromSmarts('[*:1][NH2:2]>>[*:1][NH:2]C(=O)C')
        products = rxn.RunReactants((mol,))

        if not products:
            return None


        product = products[0][0]
        Chem.SanitizeMol(product)
        product_smiles = Chem.MolToSmiles(product)


        tpsa = rdMolDescriptors.CalcTPSA(product)

        return {
            "has_amino": True,
            "product": product_smiles,
            "tpsa": round(tpsa, 2)
        }
    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == '__main__':
    test_smi = 'Nc1ccccc1'
    result = level_function(test_smi)
    print(f'Output: {result}')
