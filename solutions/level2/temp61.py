from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def level_function(mol):

    try:
        mol_obj = Chem.MolFromSmiles(mol)
        if mol_obj is None:
            return None



        node_features = []
        for atom in mol_obj.GetAtoms():
            features = [
                atom.GetAtomicNum(),
                atom.GetDegree(),
                atom.GetFormalCharge(),
                atom.GetTotalNumHs(),
                int(atom.GetIsAromatic()),
                atom.GetHybridization().real
            ]
            node_features.append(features)


        edge_index = [[], []]
        edge_attr = []
        bond_type_map = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.AROMATIC: 4
        }
        for bond in mol_obj.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            bt = bond_type_map.get(bond.GetBondType(), 0)

            edge_index[0].extend([i, j])
            edge_index[1].extend([j, i])
            edge_attr.extend([bt, bt])

        return {
            "num_nodes": mol_obj.GetNumAtoms(),
            "node_features": node_features,
            "edge_index": edge_index,
            "edge_attr": edge_attr,
            "num_edges": len(edge_attr)
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    smiles = 'c1ccccc1'
    result = level_function(smiles)
    if result:
        print(f"Output: {result['num_nodes']}{result['num_edges']}")
        print(f"Output: {result['num_nodes']}")
