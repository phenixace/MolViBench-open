from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdFMCS


def level_function(mol1, mol2):

    try:
        m1 = Chem.MolFromSmiles(mol1)
        m2 = Chem.MolFromSmiles(mol2)
        if m1 is None or m2 is None:
            return None


        formula1 = Chem.rdMolDescriptors.CalcMolFormula(m1)
        formula2 = Chem.rdMolDescriptors.CalcMolFormula(m2)
        are_isomers = formula1 == formula2

        if are_isomers:

            energies = {}
            for name, mol in [("mol1", m1), ("mol2", m2)]:
                mol_h = Chem.AddHs(mol)
                if mol_h.GetNumHeavyAtoms() > 50:
                    energies[name] = None
                    continue
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=10, params=params, maxAttempts=5)
                min_energy = float('inf')
                for cid in cids:
                    try:
                        AllChem.MMFFOptimizeMolecule(mol_h, confId=cid, maxIters=100)
                        ff = AllChem.MMFFGetMoleculeForceField(
                            mol_h, AllChem.MMFFGetMoleculeProperties(mol_h), confId=cid)
                        if ff:
                            e = ff.CalcEnergy()
                            min_energy = min(min_energy, e)
                    except Exception:
                        continue
                energies[name] = round(min_energy, 2) if min_energy < float('inf') else None

            return {
                "are_isomers": True,
                "formula": formula1,
                "mol1_min_energy": energies["mol1"],
                "mol2_min_energy": energies["mol2"],
                "lower_energy": "mol1" if (energies["mol1"] or float('inf')) < (energies["mol2"] or float('inf')) else "mol2"
            }
        else:

            mcs = rdFMCS.FindMCS([m1, m2], timeout=10)
            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
            mcs_smiles = None
            if mcs_mol:
                try:

                    match = m1.GetSubstructMatch(mcs_mol)
                    if match:
                        mcs_smiles = Chem.MolFragmentToSmiles(m1, match)
                except Exception:
                    pass

            return {
                "are_isomers": False,
                "formula1": formula1,
                "formula2": formula2,
                "MCS_smarts": mcs.smartsString,
                "MCS_smiles": mcs_smiles,
                "MCS_numAtoms": mcs.numAtoms
            }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    print('Output:', level_function('CCO', 'COC'))
    print('Output:', level_function('c1ccccc1', 'CCO'))
