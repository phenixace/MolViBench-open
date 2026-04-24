from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def level_function(mols):
    """给定一组分子，计算 scaffold 多样性指数。"""
    try:
        scaffolds = set()
        valid_count = 0

        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            valid_count += 1
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)
                scaffolds.add(scaffold_smi)
            except Exception:
                continue

        if valid_count == 0:
            return None

        diversity_index = len(scaffolds) / valid_count

        return {
            'num_molecules': valid_count,
            'num_unique_scaffolds': len(scaffolds),
            'scaffold_diversity_index': round(diversity_index, 4),
            'scaffolds': list(scaffolds)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["c1ccccc1CC", "c1ccc(CC)cc1", "c1ccncc1",
                   "CCO", "c1ccc2c(c1)cccc2"]
    result = level_function(smiles_list)
    print(f"Scaffold 多样性: {result}")
