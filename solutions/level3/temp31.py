from rdkit import Chem
from rdkit.Chem import AllChem, FindMolChiralCenters


def level_function(mol):
    """预测反应的 stereoselectivity。"""
    try:
        molecule = Chem.MolFromSmiles(mol)
        if molecule is None:
            return None

        existing_centers = FindMolChiralCenters(molecule, includeUnassigned=True)

        potential_centers = []
        for atom in molecule.GetAtoms():
            if atom.GetDegree() >= 3:
                neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                if len(set(neighbors)) >= 2 or atom.GetDegree() == 4:
                    if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                        potential_centers.append({
                            "atom_idx": atom.GetIdx(),
                            "atom_symbol": atom.GetSymbol(),
                            "neighbors": neighbors,
                        })

        result = {
            "existing_stereocenters": [
                {"atom_idx": idx, "chirality": chiral}
                for idx, chiral in existing_centers
            ],
            "potential_new_stereocenters": potential_centers,
            "num_possible_stereoisomers": 2 ** len(potential_centers) if potential_centers else 1,
        }

        return result
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    mol = "CC(=O)CC"
    print(f"Stereoselectivity 预测: {level_function(mol)}")
