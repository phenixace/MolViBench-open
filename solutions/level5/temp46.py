from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


def level_function(mols):

    try:

        mol_objs = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol_objs.append(mol)

        if len(mol_objs) < 2:
            return None


        mcs_result = rdFMCS.FindMCS(mol_objs, timeout=10)
        if mcs_result.numAtoms == 0:
            return None

        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        mcs_smarts = mcs_result.smartsString



        template = Chem.MolFromSmiles(Chem.MolToSmiles(
            Chem.MolFromSmarts(mcs_smarts)
        )) if mcs_mol else None

        new_mols = set()
        if template:
            rxns = [
                '[cH:1]>>[c:1]C',
                '[cH:1]>>[c:1]F',
                '[cH:1]>>[c:1]O',
            ]
            for rxn_smarts in rxns:
                rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                products = rxn.RunReactants((template,))
                for ps in products:
                    for p in ps:
                        try:
                            Chem.SanitizeMol(p)
                            new_mols.add(Chem.MolToSmiles(p))
                        except Exception:
                            continue


        stable_mols = []
        for smi in new_mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol_h = Chem.AddHs(mol)
            if mol_h.GetNumHeavyAtoms() > 50:
                continue
            result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDG(), maxAttempts=5)
            if result != -1:
                try:
                    energy = AllChem.UFFOptimizeMolecule(mol_h, maxIters=100)
                    if energy == 0:
                        stable_mols.append(smi)
                except Exception:
                    continue

        return {
            'mcs': mcs_smarts,
            'num_generated': len(new_mols),
            'stable_molecules': stable_mols
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles_list = ['c1ccc(O)cc1', 'c1ccc(N)cc1', 'c1ccc(C)cc1']
    result = level_function(smiles_list)
    print(f'Output: {result}')
