"""
MolViBench Evaluation Framework
================================

Prompt templates for querying LLMs.
Supports multiple prompt strategies for ablation studies.
"""

# ============================================================
# System Prompts
# ============================================================

SYSTEM_PROMPT_EN = """You are an expert cheminformatics programmer. 
Write Python code using RDKit to solve molecular computing tasks.
Your code must define a function called `level_function` that takes the specified input and returns the result.
Only use these libraries: rdkit, numpy, pandas, scikit-learn, matplotlib, selfies.
Include error handling with try/except blocks.
Do NOT include any import statements for libraries other than the ones listed above.
Return only the Python code, no explanations."""

SYSTEM_PROMPT_CN = """你是一名化学信息学编程专家。
请使用 RDKit 编写 Python 代码解决分子计算任务。
你的代码必须定义一个名为 `level_function` 的函数，接受指定输入并返回结果。
仅使用以下库：rdkit, numpy, pandas, scikit-learn, matplotlib, selfies。
包含 try/except 错误处理。
不要导入上述列表之外的库。
只返回 Python 代码，不要解释。"""


# ============================================================
# User Prompt Templates
# ============================================================

def make_prompt_basic(question: str, lang: str = "en") -> str:
    """
    Basic prompt: question only, no examples.
    This is the hardest setting — tests pure understanding.
    """
    if lang == "cn":
        return f"""请编写一个 Python 函数 `level_function` 来完成以下任务：

任务：{question}

要求：
1. 函数名必须为 `level_function`
2. 使用 RDKit 库
3. 包含 try/except 错误处理
4. 输入无效时返回 None

请直接给出完整的 Python 代码。"""
    else:
        return f"""Write a Python function `level_function` to accomplish the following task:

Task: {question}

Requirements:
1. Function must be named `level_function`
2. Use the RDKit library
3. Include try/except error handling
4. Return None for invalid input

Provide the complete Python code."""


def make_prompt_with_signature(question: str, signature: str, return_type: str,
                               lang: str = "en") -> str:
    """
    Signature prompt: includes function signature and expected return type.
    Medium difficulty — tests implementation ability.
    """
    if lang == "cn":
        return f"""请编写一个 Python 函数来完成以下任务：

任务：{question}

函数签名：{signature}
返回类型：{return_type}

要求：
1. 使用 RDKit 库
2. 包含 try/except 错误处理，异常时返回 None
3. 输出必须是确定性的（相同输入 → 相同输出）

请直接给出完整的 Python 代码。"""
    else:
        return f"""Write a Python function to accomplish the following task:

Task: {question}

Function signature: {signature}
Return type: {return_type}

Requirements:
1. Use the RDKit library
2. Include try/except error handling, return None on exception
3. Output must be deterministic (same input → same output)

Provide the complete Python code."""


def make_prompt_with_example(question: str, signature: str, return_type: str,
                             example_input: str, example_output: str,
                             lang: str = "en") -> str:
    """
    Example prompt: includes function signature, return type, and one I/O example.
    Easiest setting — tests code generation from spec.
    """
    if lang == "cn":
        return f"""请编写一个 Python 函数来完成以下任务：

任务：{question}

函数签名：{signature}
返回类型：{return_type}

示例：
  输入：{example_input}
  输出：{example_output}

要求：
1. 使用 RDKit 库
2. 包含 try/except 错误处理，异常时返回 None
3. 输出必须是确定性的

请直接给出完整的 Python 代码。"""
    else:
        return f"""Write a Python function to accomplish the following task:

Task: {question}

Function signature: {signature}
Return type: {return_type}

Example:
  Input: {example_input}
  Output: {example_output}

Requirements:
1. Use the RDKit library
2. Include try/except error handling, return None on exception
3. Output must be deterministic

Provide the complete Python code."""


# ============================================================
# Prompt strategy registry
# ============================================================

PROMPT_STRATEGIES = {
    "basic": make_prompt_basic,
    "signature": make_prompt_with_signature,
    "example": make_prompt_with_example,
}


# ============================================================
# Code extraction from LLM response
# ============================================================

def extract_code_from_response(response: str) -> str:
    """
    Extract Python code from LLM response.
    Handles markdown code blocks and plain code.
    """
    # Try to extract from ```python ... ``` blocks
    import re

    patterns = [
        r'```python\s*\n(.*?)```',
        r'```py\s*\n(.*?)```',
        r'```\s*\n(.*?)```',
    ]

    for pattern in patterns:
        matches = re.findall(pattern, response, re.DOTALL)
        if matches:
            # Return the longest match (most likely the full code)
            return max(matches, key=len).strip()

    # If no code blocks, check if the entire response looks like code
    lines = response.strip().split('\n')
    code_lines = [l for l in lines if not l.startswith('#') or l.strip().startswith('import') or l.strip().startswith('from') or l.strip().startswith('def ')]

    if any('def ' in l for l in lines) and any('import' in l or 'from ' in l for l in lines):
        return response.strip()

    return response.strip()
