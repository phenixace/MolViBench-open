from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None
        ortho_para_directors = [
            Chem.MolFromSmarts('[OH]'),
            Chem.MolFromSmarts('[NH2]'),
            Chem.MolFromSmarts('[OR]'),
            Chem.MolFromSmarts('[NR2]'),
            Chem.MolFromSmarts('[CH3]'),
            Chem.MolFromSmarts('[F]'),
            Chem.MolFromSmarts('[Cl]'),
            Chem.MolFromSmarts('[Br]'),
        ]
        meta_directors = [
            Chem.MolFromSmarts('[N+](=O)[O-]'),
            Chem.MolFromSmarts('C(=O)O'),
            Chem.MolFromSmarts('C#N'),
            Chem.MolFromSmarts('C(=O)'),
            Chem.MolFromSmarts('S(=O)(=O)'),
        ]
        ring_info = molecule.GetRingInfo()
        aromatic_atoms = [a.GetIdx() for a in molecule.GetAtoms() if a.GetIsAromatic()]
        if not aromatic_atoms:
            return None
        directing = "ortho/para"
        for pat in meta_directors:
            if pat and molecule.HasSubstructMatch(pat):
                directing = "meta"
                break
        for pat in ortho_para_directors:
            if pat and molecule.HasSubstructMatch(pat):
                directing = "ortho/para"
                break
        bromination_smarts = '[cH:1]>>[c:1]Br'
        rxn = AllChem.ReactionFromSmarts(bromination_smarts)
        products = rxn.RunReactants((molecule,))
        product_smiles = []
        if products:
            for product_set in products:
                for product in product_set:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product)
                        if smi not in product_smiles:
                            product_smiles.append(smi)
                    except Exception:
                        continue
        return {
            "directing_effect": directing,
            "possible_products": product_smiles,
            "num_positions": len(product_smiles),
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccc(O)cc1'
    result = level_function(smiles)
    print(f'Output: {result}')
