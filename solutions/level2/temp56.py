from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition

def level_function(mols, core_smiles):
    try:
        core = Chem.MolFromSmiles(core_smiles)
        if core is None:
            return None

        mol_objs = []
        for smi in mols:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                mol_objs.append(m)

        if not mol_objs:
            return None

        rg_params = rdRGroupDecomposition.RGroupDecompositionParameters()
        rg_params.removeHydrogensPostMatch = True

        rg = rdRGroupDecomposition.RGroupDecomposition(core, rg_params)
        for m in mol_objs:
            rg.Add(m)
        rg.Process()

        columns = rg.GetRGroupsAsColumns()
        result = {}
        for key, mol_list in columns.items():
            result[key] = [Chem.MolToSmiles(m) for m in mol_list]

        return result
    except Exception as e:
        print(e)
        return None
