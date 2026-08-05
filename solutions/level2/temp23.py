from rdkit import Chem
from rdkit.Chem import AllChem


def level_function(mol):

    try:
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if mol.GetNumHeavyAtoms() > 50:
            return None
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=AllChem.ETKDG(), maxAttempts=5)
        if len(conf_ids) == 0:
            return None
        best_energy = float('inf')
        best_conf_id = -1
        for conf_id in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conf_id)
            if ff is None:
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            if ff is None:
                continue
            energy = ff.CalcEnergy()
            if energy < best_energy:
                best_energy = energy
                best_conf_id = conf_id
        if best_conf_id == -1:
            return None
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        return (smiles, best_energy)
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1CCO'
    result = level_function(smiles)
    print(f'Output: {result}')
