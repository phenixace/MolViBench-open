from rdkit import Chem
from rdkit.Chem import rdFMCS

def level_function(mols):
    try:
        mol_objs = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mol_objs.append(mol)

        if len(mol_objs) < 2:
            return None

        mcs_result = rdFMCS.FindMCS(
            mol_objs,
            threshold=0.8,
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            timeout=60
        )

        if mcs_result.canceled:
            mcs_result = rdFMCS.FindMCS(
                mol_objs,
                threshold=0.7,
                timeout=60
            )

        mcs_smarts = mcs_result.smartsString
        if not mcs_smarts:
            return None

        mcs_mol = Chem.MolFromSmarts(mcs_smarts)
        mcs_smiles = None
        if mcs_mol is not None:
            try:
                mcs_smiles = Chem.MolToSmiles(mcs_mol)
            except Exception:
                mcs_smiles = None

        return {
            'mcs_smarts': mcs_smarts,
            'mcs_smiles': mcs_smiles,
            'num_atoms': mcs_result.numAtoms,
            'num_bonds': mcs_result.numBonds,
        }
    except Exception as e:
        print(e)
        return None
