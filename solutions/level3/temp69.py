from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

def level_function(mol):
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None

        mw = Descriptors.MolWt(mol_obj)
        logp = Crippen.MolLogP(mol_obj)
        tpsa = Descriptors.TPSA(mol_obj)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_obj)
        num_atoms = mol_obj.GetNumAtoms()
        mr = Crippen.MolMR(mol_obj)

        veber = rot_bonds <= 10 and tpsa <= 140

        ghose = (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and
                 40 <= num_atoms <= 480 and 20 <= mr <= 130)

        egan = tpsa <= 131.6 and logp <= 5.88

        return {
            "properties": {
                "MW": round(mw, 2),
                "LogP": round(logp, 2),
                "TPSA": round(tpsa, 2),
                "RotBonds": rot_bonds,
                "NumAtoms": num_atoms,
                "MR": round(mr, 2)
            },
            "Veber": veber,
            "Ghose": ghose,
            "Egan": egan,
            "passes_all": veber and ghose and egan
        }
    except Exception as e:
        print(e)
        return None
