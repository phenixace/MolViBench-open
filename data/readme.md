# RDKitBench: Molecular Vibe Coding Benchmark

> 系统性评估 LLM 在分子科学编程任务上的能力

**总题量：358 题（L1: 75 | L2: 72 | L3: 72 | L4: 75 | L5: 64）**

---

## Level 1: 基础级问题（分子表示和属性计算）  [75题]

**核心能力：** 分子表示转换、属性计算、子结构识别、错误处理

**评估方式：** 输入任意 SMILES → 输出确定性数值/布尔/字符串，精确匹配或数值容差比较

### 原有类型
- 分子描述符计算（MW, LogP, HBD, HBA, TPSA, QED, MR 等）
- 分子图转换（邻接矩阵、原子特征矩阵、molblock）
- 子结构匹配（芳香环、官能团、杂环等）
- 格式转换（SMILES, InChI, 2D SVG）
- 原子/键/环统计
- Lipinski Rule of Five 判断
- 手性分子判断

### 🆕 补充类型
- **SELFIES 互转：** SMILES ↔ SELFIES 表示转换
- **盐/混合物处理：** 从含盐 SMILES 中提取最大有机片段
- **立体化学标注：** R/S 绝对构型、E/Z 双键构型输出
- **同位素处理：** 含氘标记分子的精确分子量计算
- **Scaffold 提取：** Murcko scaffold 提取
- **电荷与原子环境：** Gasteiger 偏电荷、杂化类型 (sp/sp2/sp3)
- **复杂度指标：** BertzCT 复杂度、Balaban J 指数
- **大环识别：** 判断是否含 ≥12 元环
- **Error Handling：** 非法 SMILES 输入验证与错误信息
- **SMARTS 输出：** 分子的 SMARTS 模式表示
- **表面积计算：** Labute ASA
- **原子贡献分解：** Crippen MR 原子贡献
- **合成可及性：** SA Score 计算 → `float`
- **Fsp3：** sp3 碳原子占比 → `float`
- **重原子计数：** 非氢原子数 → `int`
- **芳香原子比例：** 芳香原子/总重原子 → `float`
- **环系统统计：** 脂肪环数、芳香环数 → `int`
- **额外过滤规则：** Veber 规则、Ghose 过滤 → `bool`
- **原子团统计：** NHOHCount、NOCount → `int`

---

## Level 2: 中级问题（分子变换和相似性）  [72题]

**核心能力：** 分子指纹、相似度计算、结构修饰、构象生成、文件I/O

**评估方式：** 输入固定分子/分子库 → 输出确定性列表/数值，集合比较或排序列表比较

### 原有类型
- 相似度计算（Tanimoto, Dice, Cosine, Hamming）
- Morgan 指纹生成
- 构建分子衍生物（取代、替换）
- 3D 构象生成与优化（UFF, MMFF94）
- 文件 I/O（SDF, PDB, MOL）
- 聚类与降维（PCA, t-SNE, K-means）

### 🆕 补充类型
- **多类型指纹：** MACCS Keys、RDKit 拓扑指纹、Atom Pair、Topological Torsion
- **最大公共子结构 (MCS)：** 两分子 MCS 计算
- **R-Group 分解：** 共同骨架分子的 R-group decomposition
- **分子对齐：** 基于 MCS 的 2D 对齐与可视化
- **分子片段化：** BRICS 切割、RECAP 切割
- **分子标准化：** 去盐、标准化互变异构体、标准化电荷
- **GNN 特征生成：** 转为 PyTorch Geometric 兼容格式
- **反应原子映射：** 给定反应 SMARTS 生成 atom mapping
- **批量处理：** SDF 批量读取 → DataFrame 描述符计算
- **UMAP 降维：** 指纹 UMAP 可视化
- **组合化学枚举：** 骨架 + R-group 列表 → 枚举所有组合产物 → `list[SMILES]`
- **相似性搜索：** query + 分子库 → Tanimoto Top-K → `list[(SMILES, float)]`
- **子结构数据库搜索：** SMARTS + 分子库 → 返回匹配分子及原子索引
- **肽序列转换：** 氨基酸单字母序列 → 线性肽 SMILES
- **反应枚举：** 两组反应物 + SMARTS → 枚举所有产物 → `list[SMILES]`
- **Butina 聚类：** 给定距离阈值 → 返回聚类编号 → `list[int]`
- **2D 药效团指纹：** Pharm2D fingerprint → 确定性 bit vector
- **子结构原子索引：** 分子 + SMARTS → 返回所有匹配原子索引元组

---

## Level 3: 高级问题（推理与复杂操作）  [72题]

**核心能力：** 反应模拟、立体化学、预测、高级化学变换

**评估方式：** 输入固定分子 → 输出确定性产物 SMILES/警报名称列表/指标值

### 原有类型
- 反应模拟（SMARTS 定义反应、经典有机反应）
- 异构体生成与判断（对映体、非对映、互变异构）
- 分子扰动（原子/键/侧链/环级别）
- 代谢产物/毒性/药效团预测
- PAINS/REACH/Rule of Three 过滤

### 🆕 补充类型
- **逆合成分析：** 给定产物，找可能的前体反应物
- **反应分类：** 判断反应类型（取代/加成/消除/重排）
- **反应指纹：** 计算 reaction difference fingerprint
- **保护基操作：** Boc 保护基的加入与脱除
- **Scaffold Hopping：** 保留侧链，替换核心骨架
- **Matched Molecular Pair (MMP)：** 找出两分子间的结构变换
- **ADMET 子结构规则：** BBB 通透性、CYP450 抑制预测
- **现代偶联反应：** Suzuki 偶联、Click Chemistry、Buchwald-Hartwig
- **大环处理：** 大环分子构象生成与环张力评估
- **共价抑制剂：** Michael 加成反应模拟
- **Brenk 结构警报：** 检测并列出所有匹配警报名称 → `list[str]`
- **hERG 心脏毒性：** 基于子结构规则检测 → `bool`
- **Ames 致突变性：** 结构警报检测 → `bool`
- **NIH MLSMR 警报：** 扩展版 PAINS 检测 → `bool`
- **多规则过滤：** Veber + Ghose + Egan 同时应用 → `dict{rule: bool}`
- **药效团相似度：** 2D Pharm2D 指纹相似度 → `float`
- **Scaffold 划分：** 基于 Murcko scaffold 的训练/测试集划分 → `list[list]`
- **生成质量评估：** 计算 validity/uniqueness/novelty → `tuple[float]`

---

## Level 4: 多步推理问题  [75题]

**核心能力：** 条件分支、迭代优化、回溯搜索、错误恢复、流水线处理

**评估方式：** 输入固定分子/库 → 输出确定性结果（分类统计/Top-K 列表/存活数）

### 原有类型 — if-then-modify-compute 模式 (50题)

给定分子 → 判断是否含 X → 若有 → 做 Y → 计算 Z

示例：
- 给定分子 → 判断是否含芳环 → 若有 → 替换一个氢为羟基 → 计算 LogP
- 给定分子 → 判断是否含羧基 → 若有 → 转化为酯 → 计算分子量

### 🆕 补充类型 — 多样化推理模式 (20题)

**条件分支型：**
- 给定分子 → 计算 MW → 若 MW>300 走 BRICS 分解路线，否则走片段生长路线 → 计算 LogP
- 给定分子 → 计算 LogP → 分三种情况（>3 / <0 / 正常）做不同修饰 → 调整到目标区间
- 给定分子列表 → PAINS 过滤 → 对剩余分子计算 QED → 返回 Top-3

**迭代优化型：**
- 反复添加极性基团 → 每轮检查 TPSA → 直到落入目标区间
- 迭代替换侧链 → QED 不再提升时自动停止
- 随机突变 + Lipinski 指导的贪心优化 × 10轮

**回溯搜索型：**
- 从苯出发逐步加取代基 → MW 超目标范围时回退尝试其他基团
- 片段库两两 BRICS 拼接 → 过滤无效 → 选最相似目标

**错误恢复型：**
- 解析含错误的 SMILES 列表 → 对无效的尝试修复 → 统计成功/失败
- 顺序尝试多个反应 SMARTS → 失败跳过记录 → 收集所有成功产物

**多分子协同型：**
- 两个分子互换侧链 → 生成杂交分子 → 计算与亲本的相似度
- 分子集合构建相似度图 → 找最相似对 → 计算其 MCS

**流水线处理型：**
- CSV 读取 → 过滤 → 去重 → 描述符计算 → 排序 → 导出
- SDF 读取 → 标准化 → 指纹 → 聚类 → 选代表 → 导出

**循环生成型：**
- 种子分子 → 三轮迭代（甲基取代 → 羟基取代 → 卤素取代）→ 逐轮选最优
- 片段库迭代拼接 × 3轮 → 记录 MW 增长曲线

**复杂条件型：**
- 多条件 Lipinski 决策树 → 逐条判断并输出具体数值
- 判断异构体关系 → 分支处理（是→比较构象能量 / 否→计算 MCS）
- 生成衍生物 → Lipinski 过滤 → SA Score 排序 → 输出 Top-3

**数据库搜索 + 过滤型：**
- 分子列表 → Veber + Ghose + Brenk 多规则分类（safe/risky/reject）→ 统计数量
- query + 库 → Tanimoto Top-10 → PAINS + Brenk 检查 → 返回 clean hits
- 骨架 + 3组 R-group → 组合枚举 → Lipinski → QED 排序 → Top-5
- 分子集 → Butina 聚类 → 每类提取 MCS → scaffold 频率表
- 分子库 → 级联过滤 (Lipinski→Veber→PAINS→Brenk) → 每级存活数
- query + 库 → 子结构搜索 → 描述符计算 → SA Score 排序 → Top-5

---

## Level 5: 分子发现与优化  [64题]

**核心能力：** 端到端药物发现流程、多目标优化、库设计、ML 驱动优化

**评估方式：** 给定固定输入集 → 输出确定性报告/指标值/候选列表

### 阶段一：候选生成（Hit Identification）[10题]

- 给定一组已知活性分子，生成与其相似度 >0.7 的新分子候选
- 给定一个目标 SMILES，生成所有一取代衍生物
- 给定一个 scaffold，生成 10 个不同侧链修饰的分子
- 给定一个已知药效团，生成符合要求的分子
- 给定一组分子，提取共同子结构 scaffold
- 给定蛋白结合口袋片段，生成 fragment-growing 衍生物
- 给定药物候选分子，生成环开/环合异构体
- 给定小分子，生成立体异构体并保留可行构象
- 给定 SMILES，枚举所有卤素取代体
- 给定一组分子，基于结构相似性构建候选库

### 阶段二：初步筛选（Filtering）[10题]

- Lipinski Rule of Five 筛选
- TPSA / LogP / MW / 旋转键 / QED 阈值筛选
- PAINS / 毒性子结构过滤
- CNS 药物规则 / Rule of Three 筛选

### 阶段三：优化（Lead Optimization）[10题]

- 极性/疏水性衍生物生成
- LogP / TPSA 目标值优化
- 分子量调控、QED 维持
- 合成复杂性降低
- 杂环核心替换 / 芳环取代方式探索

### 阶段四：候选对比与多目标优化  [10题]

- 相似度 vs QED 二维图
- Pareto 优化（多目标）
- 合成复杂度评分对比
- 结构多样性选择
- 聚类 + 代表分子选择
- Scaffold 多样性指数
- 相似度-药物性热图

### 阶段五：从 Hit 到 Lead 的完整流程  [10题]

- Scaffold 提取 → 衍生物 → Lipinski 筛选
- Fragment growing → PAINS 过滤 → LogP 优化
- 环异构体 → TPSA 选择 → QED 计算
- 卤代衍生物 → MW 筛选 → 合成复杂度
- 多目标优化（QED, LogP, TPSA）→ 选最佳 lead

### 🆕 阶段六：高级分析与 ML 驱动优化  [10题]

- **MaxMin 多样性选择：** 选择结构多样性最大的子集
- **BRICS 重组设计：** 切割 → 重组 → 过滤 → 评估药物性
- **Scaffold Morphing：** 系统性核心环替换 → 性质变化评估
- **QSAR 建模：** 分子指纹/描述符 → 回归/分类模型 → 活性预测
- **MMP 分析：** 找出对活性影响最大的结构变换
- **聚焦库设计：** 基于靶标活性物生成相似且多样的候选
- **遗传算法优化：** 多目标进化同时优化 QED、LogP、SA Score
- **虚拟筛选 Pipeline：** 读取库 → 描述符 → 相似性筛选 → ADMET → 排序 → 报告

### 🆕 阶段七：可评估的端到端流程  [5题]

- **组合库设计 Pipeline：** 骨架 + R-group → 枚举 → 级联过滤 → MaxMin 多样性选择 → 候选报告
- **相似性搜索 VS：** query → Top-100 搜索 → ADMET 过滤 → Butina 聚类 → 每类选最优
- **生成质量评估：** 生成 SMILES + 参考集 → validity/uniqueness/novelty/diversity → 性质分布
- **QSAR Pipeline (scaffold-split)：** 划分 → Morgan 指纹 → RF 回归 → R²/RMSE
- **Applicability Domain：** 训练集指纹 → 边界计算 → 新分子域内/域外判断

---

## 难度体系设计理念

| Level | Bloom's Taxonomy | 核心测试能力 | 推理类型 |
|-------|-----------------|-------------|---------|
| L1 | Remember & Understand | API 调用、格式转换 | 单步操作 |
| L2 | Apply | 工具组合、中等变换 | 两步组合 |
| L3 | Analyze | 化学推理、复杂操作 | 多步+领域知识 |
| L4 | Analyze & Evaluate | 条件推理、迭代优化 | 多模式推理 |
| L5 | Create | 端到端设计、优化 | 综合创造 |

## 评估维度

### 评估框架设计原则

**核心原则：给定任意合法输入，输出必须确定性可验证。**

| 输出类型 | 评估方法 | 示例 | Level 分布 |
|---------|---------|------|----------|
| `int` | 精确匹配 | 重原子数、环数、NHOHCount | L1-L2 |
| `float` | `abs(pred - ref) < ε` | SA Score、LogP、Fsp3、相似度 | L1-L3 |
| `bool` | 精确匹配 | Lipinski、Veber、Ghose、hERG | L1-L3 |
| `str` (SMILES) | canonical SMILES 等价性 | InChI、scaffold、反应产物 | L1-L3 |
| `list[str]` | 排序后集合比较 | 组合枚举产物、搜索结果 | L2-L4 |
| `list[int]` | 精确匹配 | 聚类编号、匹配原子索引 | L2-L4 |
| `dict` | key-value 逐项比较 | 多规则过滤结果、级联存活数 | L3-L5 |
| `tuple[float]` | 元素级容差 | validity/uniqueness/novelty | L3-L5 |

### 三维评估指标

| 维度 | 说明 | 权重 |
|------|------|------|
| **Executability** | 代码是否能无错运行 | 必要条件 |
| **Chemical Correctness** | 输出与 ground truth 一致（按上表规则比较） | 核心指标 |
| **Code Quality** | 可读性、鲁棒性、效率、Error handling | 辅助指标 |

### 库依赖决策

**仅使用 RDKit 及其生态库**（不引入 OpenBabel），原因：
1. 单库约束确保 output 确定性（不同库的 canonical SMILES 可能不同）
2. 降低复现门槛（无需 OpenBabel 的复杂安装）
3. 允许的辅助库：numpy, pandas, scikit-learn, matplotlib

## 与现有 Benchmark 的差异

| Benchmark | 测试什么 | RDKitBench 差异 |
|-----------|---------|----------------|
| HumanEval / MBPP | 通用编程能力 | ✅ 领域特定（化学信息学） |
| ChemBench / MoleculeNet | 化学知识/预测准确性 | ✅ 测试编程实现能力 |
| SWE-bench | 软件工程能力 | ✅ 科学计算场景 |
| SciBench | 科学推理 | ✅ 编程+化学双重评估 |
| **RDKitBench** | **Molecular Vibe Coding** | 🆕 首个分子科学编程 benchmark |
