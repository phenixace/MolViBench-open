from rdkit import Chem


def level_function(reaction_smiles):
    """判断给定反应的类型（亲核取代/亲电加成/消除/重排）。"""
    try:
        # Parse reaction SMILES (reactants>>products)
        parts = reaction_smiles.split(">>")
        if len(parts) != 2:
            return None

        reactant_smis = parts[0].split(".")
        product_smis = parts[1].split(".")

        reactants = [Chem.MolFromSmiles(s) for s in reactant_smis if Chem.MolFromSmiles(s)]
        products = [Chem.MolFromSmiles(s) for s in product_smis if Chem.MolFromSmiles(s)]

        if not reactants or not products:
            return None

        # Analyze reaction characteristics
        # Count atoms and bonds in reactants vs products
        r_heavy = sum(m.GetNumHeavyAtoms() for m in reactants)
        p_heavy = sum(m.GetNumHeavyAtoms() for m in products)

        r_bonds = sum(m.GetNumBonds() for m in reactants)
        p_bonds = sum(m.GetNumBonds() for m in products)

        r_rings = sum(Chem.rdMolDescriptors.CalcNumRings(m) for m in reactants)
        p_rings = sum(Chem.rdMolDescriptors.CalcNumRings(m) for m in products)

        # Check for leaving groups / halides in reactants
        halide_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
        r_has_halide = any(m.HasSubstructMatch(halide_pattern) for m in reactants)
        p_has_halide = any(m.HasSubstructMatch(halide_pattern) for m in products)

        # Check for double bonds
        db_pattern = Chem.MolFromSmarts("[#6]=[#6]")
        r_has_db = any(m.HasSubstructMatch(db_pattern) for m in reactants)
        p_has_db = any(m.HasSubstructMatch(db_pattern) for m in products)

        # Classification heuristics
        if r_heavy == p_heavy and len(reactants) == 1 and len(products) == 1:
            return "重排反应 (Rearrangement)"
        elif r_has_halide and not p_has_halide and len(products) >= len(reactants):
            return "亲核取代反应 (Nucleophilic Substitution)"
        elif not r_has_db and p_has_db and p_bonds < r_bonds:
            return "消除反应 (Elimination)"
        elif r_has_db and not p_has_db:
            return "亲电加成反应 (Electrophilic Addition)"
        elif r_has_halide and not p_has_halide:
            return "亲核取代反应 (Nucleophilic Substitution)"
        elif len(reactants) > len(products):
            return "加成反应 (Addition)"
        elif len(reactants) < len(products):
            return "消除反应 (Elimination)"
        else:
            return "取代反应 (Substitution)"
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    # SN2: bromoethane + hydroxide -> ethanol + bromide
    rxn = "CCBr.O>>CCO.Br"
    print(f"反应类型: {level_function(rxn)}")
