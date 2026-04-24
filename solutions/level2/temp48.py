from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def level_function(mols, perplexity=5, n_iter=500, lr=100.0):
    """对分子做 t-SNE 降维可视化（基于指纹）。"""
    try:
        fps = []
        valid_smiles = []
        for smi in mols:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(list(fp))
            valid_smiles.append(smi)
        if len(fps) < 2:
            return None

        X = np.array(fps, dtype=float)
        n = X.shape[0]
        perplexity = min(perplexity, n - 1)

        # ---------- 简易 t-SNE 实现 (Barnes-Hut 近似省略) ----------
        # 1. 计算高维相似度 (对称化的条件概率)
        def _pairwise_sq_dist(M):
            sum_sq = np.sum(M ** 2, axis=1)
            return sum_sq[:, None] + sum_sq[None, :] - 2 * M @ M.T

        D = _pairwise_sq_dist(X)
        P = np.zeros((n, n))
        for i in range(n):
            lo, hi = 1e-10, 1e4
            for _ in range(50):
                sigma = (lo + hi) / 2
                p_row = np.exp(-D[i] / (2 * sigma ** 2))
                p_row[i] = 0
                s = p_row.sum()
                if s == 0:
                    break
                p_row /= s
                entropy = -np.sum(p_row[p_row > 0] * np.log2(p_row[p_row > 0]))
                if entropy > np.log2(perplexity):
                    hi = sigma
                else:
                    lo = sigma
            P[i] = p_row
        P = (P + P.T) / (2 * n)
        P = np.maximum(P, 1e-12)

        # 2. 梯度下降
        Y = np.random.randn(n, 2) * 0.01
        for it in range(n_iter):
            d_low = _pairwise_sq_dist(Y)
            Q = 1.0 / (1.0 + d_low)
            np.fill_diagonal(Q, 0)
            Q_sum = Q.sum()
            if Q_sum == 0:
                break
            Q = Q / Q_sum
            Q = np.maximum(Q, 1e-12)
            PQ = P - Q
            grad = np.zeros_like(Y)
            for i in range(n):
                diff = Y[i] - Y
                grad[i] = 4.0 * np.sum((PQ[i] * Q[i])[:, None] * diff, axis=0)
            Y -= lr * grad

        return Y.tolist()
    except Exception as e:
        print(e)
        return None


if __name__ == "__main__":
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CCCC", "c1ccc(O)cc1",
                   "CC(C)O", "c1ccncc1", "CCN"]
    result = level_function(smiles_list)
    print(f"t-SNE 降维结果: {result}")
