from rdkit import Chem
from rdkit.Chem import rdchem

def get_metal_atomic_numbers():





    metals = set()

    metals.update([3, 11, 19, 37, 55, 87])

    metals.update([4, 12, 20, 38, 56, 88])

    metals.update(range(21, 31))
    metals.update(range(39, 49))
    metals.update(range(72, 81))

    metals.update(range(57, 72))

    metals.update(range(89, 104))

    metals.update([13, 31, 49, 50, 81, 82, 83, 84])
    return metals

def level_function(mol):



    try:
        mol = Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol
        if mol is None:
            return None

        metals = get_metal_atomic_numbers()
        metal_atoms = []

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in metals:
                metal_atoms.append(atom.GetSymbol())

        return True if metal_atoms else False

    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles1 = 'PC(N)C[C@H](F)C(=O)O'
    smiles2 = '[Na+].[O-]C(=O)C'
    smiles3 = '[Cu+2].[O-]C(=O)C'
    print(f'Output: {smiles1}{level_function(smiles1)}')
    print(f'Output: {smiles2}{level_function(smiles2)}')
    print(f'Output: {smiles3}{level_function(smiles3)}')
