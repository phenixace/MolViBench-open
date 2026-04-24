from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


def level_function(mol):
    """判断分子是否符合 Ghose 过滤规则（160≤MW≤480，-0.4≤LogP≤5.6，40≤原子数≤480，20≤MR≤130）。"""
    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None
        mw = Descriptors.MolWt(mol_obj)
        logp = Crippen.MolLogP(mol_obj)
        num_atoms = mol_obj.GetNumAtoms()  # heavy atoms
        mr = Crippen.MolMR(mol_obj)
        passes = (160 <= mw <= 480 and
                  -0.4 <= logp <= 5.6 and
                  40 <= num_atoms <= 480 and
                  20 <= mr <= 130)
        return {
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "NumAtoms": num_atoms,
            "MR": round(mr, 2),
            "passes_Ghose": passes
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles = "c1ccc(NC(=O)c2ccccc2)cc1"
    print(f"Ghose过滤: {level_function(smiles)}")
