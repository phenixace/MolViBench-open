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

### Scoring cascade

The evaluator applies the following decision order to each generated program:

1. A program that does not execute fails and cannot be rescued by API overlap.
2. Every executable output first receives type-aware value comparison.
3. If the expert solution contains a stochastic construct, a value mismatch is
   checked structurally (container/schema, numeric fields, and molecular-string
   validity). This branch does not use API overlap.
4. Only an executable result that remains intrinsically non-comparable (for
   example, image/SVG, conformer, RDKit `Mol`, array, or file/in-memory
   representation differences) can enter the strict AST/API-semantic fallback
   at `tau=0.5`. Stage-2 passes plus this narrow branch define **Main Pass@1**.
5. API coverage at `tau=0.7` is also reported for all executable Stage-2
   failures as a broader workflow-coverage diagnostic; it does not change Main
   Pass@1.

The JSON report preserves the legacy Stage-2 `pass@1` field and additionally
reports `main_pass@1` (alias `pass@1_combined`) and `fallback_pass@1` (alias
`pass@1_fallback`). Per-question records state whether exact, structural, strict
API, or diagnostic evaluation was used.

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

# DG temperature sweep (reads the key without exposing it in the process list)
python inference_dg.py --base-url https://api.openai.com/v1 --api-key-file key.txt \
  --model gpt-4o --temperatures 0.3 0.7 1.0 --seed 42

# IR — Incremental Repair (3 rounds)
python inference_ir.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o --max-rounds 3

# AC — Agent Collaboration
python inference_ac.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o
```

### API-semantic evaluation

The API-semantic judge makes one API request per solution and reports two
decisions separately: whether the candidate solves the task, and whether it is
functionally equivalent to the reference implementation.

Its defaults are locked to `op-4.7` at `https://apicursor.com/v1`, read the key
from `key.txt`, and select only
`predictions/prev_predictions/*_en_basic` directories while excluding every
model name containing `glm`. Inspect the exact scope without making requests:

```bash
python api_semantic_eval.py --dry-run
```

Run with 30 concurrent single-solution requests:

```bash
python api_semantic_eval.py --yes --concurrency 30
```

Select the paper's English Incremental Repair or Agent Collaboration outputs
without changing the one-solution-per-request behavior:

```bash
python api_semantic_eval.py --variant ir --yes --concurrency 30
python api_semantic_eval.py --variant ac --yes --concurrency 30
```

Every completed judgment is appended to a JSONL checkpoint. Re-running the
same command resumes unfinished/API-error items without requesting successful
items again. The consolidated report is written to
`results/api_semantic_prev_predictions_en_basic_op-4.7.json`, with a compact
per-model comparison in the adjacent `_summary.csv` file.

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
