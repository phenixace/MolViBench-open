from rdkit import Chem
from rdkit.Chem import AllChem

def level_function(mol):



    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mol_h = Chem.AddHs(mol_obj)
        if mol_h.GetNumHeavyAtoms() > 50:
            return None

        result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG(), maxAttempts=5)
        if result == -1:
            return None

        opt_result = AllChem.UFFOptimizeMolecule(mol_h, maxIters=100)

        conf = mol_h.GetConformer()
        coords = []
        for i in range(mol_h.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append((mol_h.GetAtomWithIdx(i).GetSymbol(), pos.x, pos.y, pos.z))
        return coords
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CCO'
    result = level_function(smiles)
    if result:
        print('Output')
        for atom_symbol, x, y, z in result:
            print(f'Output: {atom_symbol}{x:.4f}{y:.4f}{z:.4f}')
