from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def level_function(mol_smi, new_scaffold_smiles):




    try:
        mol = Chem.MolFromSmiles(mol_smi)

        new_scaffold = Chem.MolFromSmiles(new_scaffold_smiles)

        if mol is None or new_scaffold is None:
            return None


        scaffold = MurckoScaffold.GetScaffoldForMol(mol)



        side_chains = Chem.ReplaceCore(mol, scaffold)

        if side_chains is None:

            return {"hopped_molecules": [Chem.MolToSmiles(new_scaffold)]}




        query = Chem.MolFromSmiles('[*]')




        frags = Chem.GetMolFrags(side_chains, asMols=True)

        current_mol = new_scaffold
        for frag in frags:

            replaced = AllChem.ReplaceSubstructs(current_mol, query, frag, replacementConnectionPoint=0)
            if replaced:
                current_mol = replaced[0]


        final_mol = Chem.DeleteSubstructs(current_mol, query)


        final_smi = Chem.MolToSmiles(final_mol)
        if '.' in final_smi:

            final_smi = max(final_smi.split('.'), key=len)

        return {
            "original_scaffold": Chem.MolToSmiles(scaffold),
            "new_scaffold": Chem.MolToSmiles(new_scaffold),
            "hopped_molecules": [final_smi]
        }

    except Exception as e:
        print(f"Error logic: {e}")
        return None

if __name__ == '__main__':
    mol_smi = 'CC(=O)Nc1ccccc1'
    new_scaf = '[*:1]c1ccncc1'
    result = level_function(mol_smi, new_scaf)
    if result:
        print(f"Output: {result['original_scaffold']}")
        print(f"Output: {result['new_scaffold']}")
        print(f"Output: {result['hopped_molecules']}")
