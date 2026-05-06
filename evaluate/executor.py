"""
MolViBench Evaluation Framework
================================

Safe code execution with timeout and sandboxing.
Executes LLM-generated code in isolated subprocess.
"""

import subprocess
import sys
import os
import json
import tempfile
import time
import traceback
import importlib.util
from typing import Any, Optional, Tuple


# Maximum execution time per test case (seconds)
DEFAULT_TIMEOUT = 60


def load_function_from_file(filepath: str, func_name: str = "level_function") -> Optional[callable]:
    """
    Dynamically load a function from a Python file.
    Returns the function object, or None if not found.
    """
    try:
        spec = importlib.util.spec_from_file_location("dynamic_module", filepath)
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        func = getattr(module, func_name, None)
        return func
    except Exception as e:
        return None


def execute_function_direct(filepath: str, func_name: str, args: list, kwargs: dict,
                            timeout: int = DEFAULT_TIMEOUT) -> dict:
    """
    Execute a function from a file directly in-process.
    Fast but no isolation — use for trusted solution code.

    Returns:
        {
            "success": bool,
            "output": Any,
            "error": str or None,
            "time_s": float
        }
    """
    func = load_function_from_file(filepath, func_name)
    if func is None:
        return {
            "success": False,
            "output": None,
            "error": f"Function '{func_name}' not found in {filepath}",
            "time_s": 0.0,
        }

    start = time.time()
    try:
        output = func(*args, **kwargs)
        elapsed = time.time() - start
        return {
            "success": True,
            "output": output,
            "error": None,
            "time_s": round(elapsed, 4),
        }
    except Exception as e:
        elapsed = time.time() - start
        return {
            "success": False,
            "output": None,
            "error": f"{type(e).__name__}: {e}",
            "time_s": round(elapsed, 4),
        }


def execute_function_subprocess(filepath: str, func_name: str, args: list, kwargs: dict,
                                timeout: int = DEFAULT_TIMEOUT) -> dict:
    """
    Execute a function from a file in an isolated subprocess.
    Safe for untrusted LLM-generated code. Uses JSON for I/O.

    Returns:
        {
            "success": bool,
            "output": Any (JSON-serializable),
            "error": str or None,
            "time_s": float,
            "stdout": str,
            "stderr": str,
            "syntax_error": bool,
            "runtime_error": bool,
            "timeout": bool,
        }
    """
    # Create a wrapper script that imports and calls the function
    wrapper_code = f'''
import sys, json, traceback, os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(r"{filepath}"))

try:
    import importlib.util
    spec = importlib.util.spec_from_file_location("target_module", r"{filepath}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
except SyntaxError as e:
    result = {{"success": False, "output": None, "error": f"SyntaxError: {{e}}", "syntax_error": True, "runtime_error": False}}
    print("__RESULT__" + json.dumps(result, default=str))
    sys.exit(0)
except Exception as e:
    result = {{"success": False, "output": None, "error": f"ImportError: {{e}}", "syntax_error": False, "runtime_error": True}}
    print("__RESULT__" + json.dumps(result, default=str))
    sys.exit(0)

func = getattr(module, "{func_name}", None)
if func is None:
    # Try common alternative names
    for alt in ["level_function", "solution", "solve", "main"]:
        func = getattr(module, alt, None)
        if func is not None:
            break

if func is None:
    result = {{"success": False, "output": None, "error": "Function '{func_name}' not found", "syntax_error": False, "runtime_error": True}}
    print("__RESULT__" + json.dumps(result, default=str))
    sys.exit(0)

# Deserialize arguments
args = json.loads(r"""{json.dumps(args, default=str)}""")
kwargs = json.loads(r"""{json.dumps(kwargs, default=str)}""")

try:
    output = func(*args, **kwargs)
    result = {{"success": True, "output": output, "error": None, "syntax_error": False, "runtime_error": False}}
except Exception as e:
    tb = traceback.format_exc()
    result = {{"success": False, "output": None, "error": f"{{type(e).__name__}}: {{e}}\\n{{tb}}", "syntax_error": False, "runtime_error": True}}

print("__RESULT__" + json.dumps(result, default=str))
'''

    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, encoding='utf-8') as f:
        f.write(wrapper_code)
        wrapper_path = f.name

    start = time.time()
    try:
        result = subprocess.run(
            [sys.executable, wrapper_path],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=os.path.dirname(filepath),
            env={**os.environ, "PYTHONIOENCODING": "utf-8"},
        )
        elapsed = time.time() - start

        stdout = result.stdout
        stderr = result.stderr

        # Extract result from stdout
        result_line = None
        other_stdout = []
        for line in stdout.split('\n'):
            if line.startswith("__RESULT__"):
                result_line = line[len("__RESULT__"):]
            else:
                other_stdout.append(line)

        if result_line:
            parsed = json.loads(result_line)
            return {
                "success": parsed.get("success", False),
                "output": parsed.get("output"),
                "error": parsed.get("error"),
                "time_s": round(elapsed, 4),
                "stdout": '\n'.join(other_stdout).strip(),
                "stderr": stderr.strip(),
                "syntax_error": parsed.get("syntax_error", False),
                "runtime_error": parsed.get("runtime_error", False),
                "timeout": False,
            }
        else:
            return {
                "success": False,
                "output": None,
                "error": f"No result returned. stderr: {stderr[:500]}",
                "time_s": round(elapsed, 4),
                "stdout": stdout.strip(),
                "stderr": stderr.strip(),
                "syntax_error": False,
                "runtime_error": True,
                "timeout": False,
            }

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        return {
            "success": False,
            "output": None,
            "error": f"Timeout after {timeout}s",
            "time_s": round(elapsed, 4),
            "stdout": "",
            "stderr": "",
            "syntax_error": False,
            "runtime_error": False,
            "timeout": True,
        }
    except Exception as e:
        elapsed = time.time() - start
        return {
            "success": False,
            "output": None,
            "error": f"Execution error: {e}",
            "time_s": round(elapsed, 4),
            "stdout": "",
            "stderr": str(e),
            "syntax_error": False,
            "runtime_error": True,
            "timeout": False,
        }
    finally:
        try:
            os.unlink(wrapper_path)
        except OSError:
            pass


def check_syntax(filepath: str) -> Tuple[bool, str]:
    """Check if a Python file has valid syntax."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            source = f.read()
        compile(source, filepath, 'exec')
        return True, "OK"
    except SyntaxError as e:
        return False, f"Line {e.lineno}: {e.msg}"
    except Exception as e:
        return False, str(e)


def check_function_exists(filepath: str, func_name: str = "level_function") -> Tuple[bool, str]:
    """Check if the expected function is defined in the file."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            source = f.read()

        # Simple AST check
        import ast
        tree = ast.parse(source)
        func_names = [
            node.name for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ]

        if func_name in func_names:
            return True, f"Found '{func_name}'"

        # Check alternatives
        alternatives = ["level_function", "solution", "solve"]
        found = [n for n in func_names if n in alternatives]
        if found:
            return True, f"Found alternative: '{found[0]}' (expected '{func_name}')"

        return False, f"Function '{func_name}' not found. Defined: {func_names}"
    except Exception as e:
        return False, str(e)
