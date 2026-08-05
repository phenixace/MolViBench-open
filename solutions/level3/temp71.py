from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import random


def level_function(smiles_list, test_ratio=0.2, random_seed=42):

    try:
        random.seed(random_seed)


        scaffold_to_mols = {}
        valid_smiles = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            canonical = Chem.MolToSmiles(mol)
            valid_smiles.append(canonical)

            scaffold = MurckoScaffold.MurckoScaffoldSmiles(
                mol=mol, includeChirality=False
            )
            if scaffold not in scaffold_to_mols:
                scaffold_to_mols[scaffold] = []
            scaffold_to_mols[scaffold].append(canonical)

        if not valid_smiles:
            return None


        scaffolds = sorted(scaffold_to_mols.keys(),
                          key=lambda s: len(scaffold_to_mols[s]),
                          reverse=True)


        n_total = len(valid_smiles)
        n_test_target = int(n_total * test_ratio)


        scaffold_list = list(scaffolds)
        random.shuffle(scaffold_list)

        test_set = []
        train_set = []
        test_scaffolds = set()

        for scaffold in scaffold_list:
            mols = scaffold_to_mols[scaffold]
            if len(test_set) < n_test_target:
                test_set.extend(mols)
                test_scaffolds.add(scaffold)
            else:
                train_set.extend(mols)

        return {
            "total_molecules": n_total,
            "num_scaffolds": len(scaffolds),
            "train_size": len(train_set),
            "test_size": len(test_set),
            "train_smiles": train_set,
            "test_smiles": test_set,
            "test_scaffolds": sorted(list(test_scaffolds))
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    mols = ['c1ccccc1O', 'c1ccccc1N', 'c1ccccc1F', 'c1ccncc1O', 'c1ccncc1N', 'CCCCO', 'CCCCN', 'CCCC(=O)O', 'c1ccc2ccccc2c1', 'c1ccc2ccccc2c1O']
    result = level_function(mols, test_ratio=0.3, random_seed=42)
    if result:
        print(f"Output: {result['train_size']}{result['test_size']}")
        print(f"Output: {result['num_scaffolds']}")
