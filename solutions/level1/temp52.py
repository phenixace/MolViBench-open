from rdkit import Chem
from rdkit.Chem import SaltRemover


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        remover = SaltRemover.SaltRemover()
        stripped = remover.StripMol(mol_obj)
        if stripped.GetNumAtoms() == 0:

            frags = Chem.GetMolFrags(mol_obj, asMols=True, sanitizeFrags=True)
            if not frags:
                return None
            largest = max(frags, key=lambda m: m.GetNumAtoms())
            return Chem.MolToSmiles(largest)
        return Chem.MolToSmiles(stripped)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = '[Na+].[Cl-].CCO'
    print(f'Output: {level_function(smiles)}')
