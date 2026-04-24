from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def level_function(mols):
    """使用 UMAP 对一组分子指纹进行降维可视化。"""
    try:
        from sklearn.manifold import TSNE  # Fallback if umap not available

        fps = []
        valid_smiles = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(np.array(fp))
                valid_smiles.append(Chem.MolToSmiles(mol))

        if len(fps) < 3:
            return None

        X = np.array(fps)

        # Try UMAP first, fall back to t-SNE
        try:
            import umap
            reducer = umap.UMAP(n_components=2, random_state=42)
            embedding = reducer.fit_transform(X)
        except ImportError:
            perplexity = min(30, len(fps) - 1)
            reducer = TSNE(n_components=2, random_state=42, perplexity=perplexity)
            embedding = reducer.fit_transform(X)

        result = []
        for i, smi in enumerate(valid_smiles):
            result.append({
                "smiles": smi,
                "x": round(float(embedding[i, 0]), 4),
                "y": round(float(embedding[i, 1]), 4)
            })
        return result
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["c1ccccc1", "CCO", "CC(=O)O", "c1ccc(O)cc1", "CCCC",
                   "CC(C)C", "c1ccncc1", "CC=O", "CCC(=O)O", "CCCCO"]
    result = level_function(smiles_list)
    if result:
        for r in result[:3]:
            print(r)
