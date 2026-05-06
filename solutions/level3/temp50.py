from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mw = Descriptors.MolWt(mol_obj)
        logp = Descriptors.MolLogP(mol_obj)
        hbd = rdMolDescriptors.CalcNumHBD(mol_obj)
        hba = rdMolDescriptors.CalcNumHBA(mol_obj)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_obj)
        tpsa = rdMolDescriptors.CalcTPSA(mol_obj)

        passes = (mw <= 300 and logp <= 3 and hbd <= 3 and
                  hba <= 3 and rot_bonds <= 3 and tpsa <= 60)

        return {
            "passes_ro3": passes,
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "HBD": hbd,
            "HBA": hba,
            "RotBonds": rot_bonds,
            "TPSA": round(tpsa, 2)
        }
    except Exception as e:
        print(e)
        return None
