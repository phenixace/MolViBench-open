from rdkit import Chem

def level_function(filename):
    try:
        supplier = Chem.SDMolSupplier(filename)
        smiles_list = []
        for mol in supplier:
            if mol is not None:
                smiles_list.append(Chem.MolToSmiles(mol))
        return smiles_list
    except Exception as e:
        print(e)
        return None
