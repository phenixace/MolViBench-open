import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors


def level_function(seed_smiles, population_size=10, n_generations=5, seed=42):

    try:
        np.random.seed(seed)


        _sascorer = None
        try:
            import os as _os
            from rdkit.Chem import RDConfig
            import importlib.util as _ilu
            _sa_path = _os.path.join(RDConfig.RDContribDir, 'SA_Score', 'sascorer.py')
            _spec = _ilu.spec_from_file_location("sascorer", _sa_path)
            _sascorer = _ilu.module_from_spec(_spec)
            _spec.loader.exec_module(_sascorer)
        except Exception:
            _sascorer = None

        def get_sa_score(mol):

            try:
                if _sascorer is not None:
                    return _sascorer.calculateScore(mol)
                return Descriptors.BertzCT(mol) / 100.0
            except Exception:
                return 5.0

        def fitness(mol):

            try:
                qed = Descriptors.qed(mol)
                logp = Crippen.MolLogP(mol)
                sa = get_sa_score(mol)


                if 1.0 <= logp <= 3.0:
                    logp_score = 1.0
                else:
                    logp_score = max(0, 1.0 - abs(logp - 2.0) / 5.0)


                sa_score = max(0, 1.0 - (sa - 1.0) / 9.0)

                return 0.4 * qed + 0.3 * logp_score + 0.3 * sa_score
            except Exception:
                return 0.0

        def mutate(mol):

            try:
                rwmol = Chem.RWMol(mol)
                n_atoms = rwmol.GetNumAtoms()
                if n_atoms < 2:
                    return mol

                mutation_type = np.random.choice(['add_atom', 'remove_atom', 'change_bond', 'change_atom'])

                if mutation_type == 'add_atom':
                    atom_types = [6, 7, 8, 9, 17]
                    new_atom = Chem.Atom(atom_types[np.random.randint(len(atom_types))])
                    idx = rwmol.AddAtom(new_atom)
                    attach = np.random.randint(n_atoms)
                    rwmol.AddBond(attach, idx, Chem.BondType.SINGLE)

                elif mutation_type == 'remove_atom' and n_atoms > 5:
                    idx = np.random.randint(n_atoms)
                    atom = rwmol.GetAtomWithIdx(idx)
                    if atom.GetDegree() == 1:
                        rwmol.RemoveAtom(idx)

                elif mutation_type == 'change_bond':
                    bonds = list(rwmol.GetBonds())
                    if bonds:
                        bond = bonds[np.random.randint(len(bonds))]
                        if not bond.IsInRing():
                            bt = np.random.choice([Chem.BondType.SINGLE, Chem.BondType.DOUBLE])
                            bond.SetBondType(bt)

                elif mutation_type == 'change_atom':
                    idx = np.random.randint(n_atoms)
                    atom_types = [6, 7, 8]
                    rwmol.GetAtomWithIdx(idx).SetAtomicNum(atom_types[np.random.randint(len(atom_types))])

                Chem.SanitizeMol(rwmol)
                return rwmol.GetMol()
            except Exception:
                return mol

        def crossover(mol1, mol2):

            try:
                combined = Chem.CombineMols(mol1, mol2)
                rwmol = Chem.RWMol(combined)
                n1 = mol1.GetNumAtoms()
                n2 = mol2.GetNumAtoms()
                if n1 > 0 and n2 > 0:
                    idx1 = np.random.randint(n1)
                    idx2 = n1 + np.random.randint(n2)
                    rwmol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
                Chem.SanitizeMol(rwmol)
                return rwmol.GetMol()
            except Exception:
                return mol1


        seed_mol = Chem.MolFromSmiles(seed_smiles)
        if seed_mol is None:
            return None

        population = [seed_mol]
        for _ in range(population_size - 1):
            m = mutate(Chem.RWMol(seed_mol))
            if m is not None:
                population.append(m)
            else:
                population.append(seed_mol)

        history = []
        for gen in range(n_generations):

            scored = [(mol, fitness(mol)) for mol in population]
            scored.sort(key=lambda x: x[1], reverse=True)

            best_mol, best_fit = scored[0]
            best_smi = Chem.MolToSmiles(best_mol)

            history.append({
                "generation": gen,
                "best_fitness": round(best_fit, 4),
                "best_smiles": best_smi,
                "mean_fitness": round(float(np.mean([s[1] for s in scored])), 4),
            })


            survivors = []
            seen = set()
            for mol, _ in scored:
                smi = Chem.MolToSmiles(mol)
                if smi not in seen:
                    seen.add(smi)
                    survivors.append(mol)
                if len(survivors) >= population_size // 2:
                    break

            if not survivors:
                survivors = [scored[0][0]]


            new_pop = list(survivors)
            while len(new_pop) < population_size:
                if np.random.random() < 0.7:

                    parent = survivors[np.random.randint(len(survivors))]
                    child = mutate(parent)
                else:

                    p1 = survivors[np.random.randint(len(survivors))]
                    p2 = survivors[np.random.randint(len(survivors))]
                    child = crossover(p1, p2)
                if child is not None:
                    new_pop.append(child)
                else:
                    new_pop.append(survivors[0])

            population = new_pop[:population_size]


        final_scored = [(mol, fitness(mol)) for mol in population]
        final_scored.sort(key=lambda x: x[1], reverse=True)

        top_results = []
        seen_smiles = set()
        for mol, fit in final_scored:
            smi = Chem.MolToSmiles(mol)
            if smi in seen_smiles:
                continue
            seen_smiles.add(smi)
            top_results.append({
                "smiles": smi,
                "fitness": round(fit, 4),
                "QED": round(Descriptors.qed(mol), 4),
                "LogP": round(Crippen.MolLogP(mol), 4),
                "SA_score": round(get_sa_score(mol), 4),
            })
            if len(top_results) >= 10:
                break

        return {
            "seed": seed_smiles,
            "n_generations": n_generations,
            "population_size": population_size,
            "initial_fitness": round(history[0]["best_fitness"], 4),
            "final_fitness": round(history[-1]["best_fitness"], 4),
            "improvement": round(history[-1]["best_fitness"] - history[0]["best_fitness"], 4),
            "top_molecules": top_results,
            "evolution_history": history[::5],
        }
    except Exception as e:
        print(e)
        return None


if __name__ == '__main__':
    result = level_function('c1ccc(NC(=O)C)cc1', population_size=30, n_generations=20)
    if result:
        print(f"Output: {result['initial_fitness']}{result['final_fitness']}")
        for m in result['top_molecules'][:3]:
            print(f"Output: {m['smiles']}{m['QED']}{m['LogP']}{m['SA_score']}")
