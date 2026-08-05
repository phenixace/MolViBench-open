from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None


        num_atoms = mol_obj.GetNumAtoms()
        has_many_atoms = num_atoms > 20

        if not has_many_atoms:
            return None


        scaffold = MurckoScaffold.GetScaffoldForMol(mol_obj)
        Chem.SanitizeMol(scaffold)
        scaffold_smiles = Chem.MolToSmiles(scaffold)


        mol_wt = rdMolDescriptors.CalcExactMolWt(scaffold)

        return {
            "num_atoms": num_atoms,
            "product": scaffold_smiles,
            "molecular_weight": round(mol_wt, 4)
        }
    except Exception as e:
        print(e)
        return None

if __name__ == '__main__':
    smiles = 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'
    print(f'Output: {level_function(smiles)}')
