from rdkit import Chem
from collections import Counter

def level_function(reactants, products):

    def get_mol_stats(smiles_list):
        total_atoms = Counter()
        total_charge = 0

        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None, None

            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()

                total_atoms[symbol] += 1
                total_atoms['H'] += atom.GetTotalNumHs()

                total_charge += atom.GetFormalCharge()
        return total_atoms, total_charge

    r_atoms, r_charge = get_mol_stats(reactants)
    p_atoms, p_charge = get_mol_stats(products)

    if r_atoms is None or p_atoms is None:
        return None


    return r_atoms == p_atoms and r_charge == p_charge

if __name__ == '__main__':
    reactants = ['CC(=O)O', 'CCO']
    products = ['CC(=O)OCC', 'O']
    print(f'Output: {level_function(reactants, products)}')
