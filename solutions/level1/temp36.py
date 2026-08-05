from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def level_function(mol1, mol2, consider_stereochemistry=True):









    try:
        mol1 = Chem.MolFromSmiles(mol1)
        mol2 = Chem.MolFromSmiles(mol2)
        if mol1 is None or mol2 is None:
            return None


        formula1 = rdMolDescriptors.CalcMolFormula(mol1)
        formula2 = rdMolDescriptors.CalcMolFormula(mol2)
        if formula1 != formula2:
            return False


        smiles1 = Chem.MolToSmiles(mol1, isomericSmiles=consider_stereochemistry, canonical=True)
        smiles2 = Chem.MolToSmiles(mol2, isomericSmiles=consider_stereochemistry, canonical=True)

        if smiles1 != smiles2:
            return True
        else:
            return False

    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles1 = 'CC[C@H](F)C(=O)O'
    smiles2 = 'CC[C@@H](F)C(=O)O'
    smiles3 = 'CCC(F)C(=O)O'
    print('Output')
    print(f'Output: {smiles1}{smiles2}{level_function(smiles1, smiles2, consider_stereochemistry=True)}')
    print(f'Output: {smiles1}{smiles3}{level_function(smiles1, smiles3, consider_stereochemistry=True)}')
    print('Output')
    print(f'Output: {smiles1}{smiles2}{level_function(smiles1, smiles2, consider_stereochemistry=False)}')
    print(f'Output: {smiles1}{smiles3}{level_function(smiles1, smiles3, consider_stereochemistry=False)}')
