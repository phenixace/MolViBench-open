# MolViBench

MolViBench is a **molecular vibe coding benchmark** for evaluating LLMs on generating executable RDKit-based code from natural-language chemistry tasks.

## Overview

- **358 bilingual tasks** (CN/EN) across **5 difficulty levels**
- **Reference solutions** for all tasks
- **Automated evaluation framework** with deterministic verification

| Level | Tasks | Core Capability | Bloom's Taxonomy |
|-------|-------|----------------|-----------------|
| L1 | 75 | Molecular representation & property calculation | Remember & Understand |
| L2 | 72 | Molecular transformation & similarity | Apply |
| L3 | 72 | Reasoning & complex operations | Analyze |
| L4 | 75 | Multi-step reasoning | Analyze & Evaluate |
| L5 | 64 | Molecular discovery & optimization | Create |

## Quick Start

```bash
pip install -r requirements.txt
```

## Evaluation

Validate that all reference solutions can run:

```bash
python test.py --validate-solutions --lang cn
```

Generate a ground-truth cache:

```bash
python test.py --generate-ground-truth --lang cn --n-test 10
```

Evaluate all levels (expects predictions in `predictions/`):

```bash
python test.py --all --lang cn
```

Evaluate a single level or question:

```bash
python test.py --level 1 --lang cn
python test.py --level 1 --question 5 --lang cn
```

## Data

- Task CSVs: `data/cn/level{1..5}.csv`, `data/en/level{1..5}.csv`
- Dataset description: `data/readme.md`
- Croissant metadata: `croissant.json`

## Inference

Three paradigms for querying LLMs (all support OpenAI-compatible APIs):

| Script | Paradigm | Description |
|--------|----------|-------------|
| `inference_dg.py` | Direct Generation (DG) | Single-pass: task → code |
| `inference_sr.py` | Self-Repair (SR) | Iterative: task → code → execute → error feedback → fix (up to N rounds) |
| `inference_ac.py` | Agent Collaboration (AC) | Two-agent: Coder generates code, Tester designs tests & judges correctness |

```bash
# DG — direct generation
python inference_dg.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o

# SR — self-repair (3 rounds)
python inference_sr.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o --max-rounds 3

# AC — agent collaboration
python inference_ac.py --base-url https://api.openai.com/v1 --api-key sk-xxx --model gpt-4o
```

## Citation

If you use MolViBench in your research, please cite:

```bibtex
@inproceedings{molvibench2025,
  title={MolViBench: A Molecular Vibe Coding Benchmark for Evaluating LLMs},
  author={Anonymous},
  booktitle={NeurIPS 2025 Datasets and Benchmarks Track},
  year={2025}
}
```

## License

MIT License
