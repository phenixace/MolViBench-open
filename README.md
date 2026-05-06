# MolViBench: Evaluating LLMs on Molecular Vibe Coding

[![arXiv](https://img.shields.io/badge/arXiv-2605.02351-b31b1b.svg)](https://arxiv.org/abs/2605.02351v2)

**MolViBench** is the first benchmark tailored for *Molecular Vibe Coding* — a paradigm where chemists interact with LLMs to generate executable programs for molecular tasks.

MolViBench comprises **358 curated tasks** across five cognitive levels (following Bloom's taxonomy), ranging from single-API recall to end-to-end virtual screening pipeline design, spanning **12 real-world drug discovery workflows**.

> **Paper:** [MolViBench: Evaluating LLMs on Molecular Vibe Coding](https://arxiv.org/abs/2605.02351v2)
>
> **Authors:** Jiatong Li, Yuxuan Ren, Weida Wang, Changmeng Zheng, Xiao-yong Wei†, Qing Li, Yatao Bian
>
> The Hong Kong Polytechnic University · National University of Singapore · Fudan University

## Overview

- **358 tasks** (EN/CN) across **5 difficulty levels**
- **Reference solutions** for all tasks
- **Automated evaluation framework** with type-aware output comparison and AST-based API-semantic fallback analysis

| Level | Tasks | Cognitive Level | Core Capability |
|-------|-------|-----------------|-----------------|
| L1 | 75 | Remember & Understand | Molecular representation & property calculation |
| L2 | 72 | Apply | Molecular transformation & similarity |
| L3 | 72 | Analyze | Reasoning & complex operations |
| L4 | 75 | Analyze & Evaluate | Multi-step reasoning with control flow |
| L5 | 64 | Create | End-to-end molecular discovery pipelines |

## Installation

```bash
pip install -r requirements.txt
```

All tasks are restricted to **RDKit** as the primary cheminformatics library, with NumPy, pandas, scikit-learn, matplotlib, and selfies as permitted auxiliary dependencies.

## Evaluation

Validate that all reference solutions can run:

```bash
python test.py --validate-solutions --lang en
```

Generate a ground-truth cache:

```bash
python test.py --generate-ground-truth --lang en --n-test 10
```

Evaluate all levels (expects predictions in `predictions/`):

```bash
python test.py --all --lang en
```

Evaluate a single level or question:

```bash
python test.py --level 1 --lang en
python test.py --level 1 --question 5 --lang en
```

## Data

- Task CSVs: `data/en/level{1..5}.csv`, `data/cn/level{1..5}.csv`
- Each task specifies a natural-language instruction; the model must generate a self-contained `level_function()` that processes molecular inputs and returns a chemically meaningful result.

## Inference Paradigms

Three paradigms for querying LLMs (all support OpenAI-compatible APIs):

| Script | Paradigm | Description |
|--------|----------|-------------|
| `inference_dg.py` | Direct Generation (DG) | Single-pass: task → code |
| `inference_ir.py` | Incremental Repair (IR) | Iterative: task → code → execute → error feedback → fix (up to N rounds) |
| `inference_ac.py` | Agent Collaboration (AC) | Two-agent: Coder generates code, Tester designs tests & judges correctness |

```bash
# DG — Direct Generation
python inference_dg.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o

# IR — Incremental Repair (3 rounds)
python inference_ir.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o --max-rounds 3

# AC — Agent Collaboration
python inference_ac.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o
```

### Hyperparameters

| Parameter | Direct | Incremental Repair | Agent Collaboration |
|-----------|--------|--------------------|---------------------|
| Temperature | 0.0 | 0.0 | 0.0 |
| Max tokens | 4096 | 4096 | 4096 |
| API timeout (s) | 120 | 120 | 120 |
| Execution timeout (s) | — | 30 | 30 |
| Max rounds | — | 3 | 3 |

## Citation

If you find MolViBench useful in your research, please cite:

```bibtex
@article{li2026molvibench,
  title={MolViBench: Evaluating LLMs on Molecular Vibe Coding},
  author={Li, Jiatong and Ren, Yuxuan and Wang, Weida and Zheng, Changmeng and Wei, Xiao-yong and Li, Qing and Bian, Yatao},
  journal={arXiv preprint arXiv:2605.02351},
  year={2026}
}
```

## License

MIT License
