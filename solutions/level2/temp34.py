from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None

        if mol.GetNumHeavyAtoms() > 50:
            return None
        mol = Chem.AddHs(mol)
        res = AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=5)
        if res == -1:
            return None
        AllChem.MMFFOptimizeMolecule(mol, maxIters=100)
        conf = mol.GetConformer()
        coords = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords.append((atom.GetSymbol(), round(pos.x, 4), round(pos.y, 4), round(pos.z, 4)))
        return coords
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('CCO')
    print(f'Output: {result}')
