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
import inspect
from typing import Any, Optional, Tuple


# Maximum execution time per test case (seconds)
DEFAULT_TIMEOUT = 60


def _seed_random_generators(seed: Optional[int]) -> None:
    """Seed common process-global RNGs before loading/executing a solution."""
    if seed is None:
        return
    import random
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except Exception:
        pass


def _with_seed_kwarg(func: callable, args: list, kwargs: dict,
                      seed: Optional[int]) -> dict:
    """Pass the case seed when a function explicitly exposes a seed parameter."""
    if seed is None or "seed" in kwargs:
        return kwargs
    try:
        signature = inspect.signature(func)
        parameter = signature.parameters.get("seed")
        if parameter is None or parameter.kind not in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        ):
            return kwargs
        if "seed" in signature.bind_partial(*args, **kwargs).arguments:
            return kwargs
    except (TypeError, ValueError):
        return kwargs
    seeded_kwargs = dict(kwargs)
    seeded_kwargs["seed"] = seed
    return seeded_kwargs


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
                            timeout: int = DEFAULT_TIMEOUT,
                            seed: Optional[int] = None) -> dict:
    """
    Execute a function from a file directly in-process.
    Uses threading with daemon threads — if timeout, we abandon the thread and move on.
    The daemon thread will be killed when the main process exits.
    """
    _seed_random_generators(seed)
    func = load_function_from_file(filepath, func_name)
    if func is None:
        return {
            "success": False,
            "output": None,
            "error": f"Function '{func_name}' not found in {filepath}",
            "time_s": 0.0,
        }
    call_kwargs = _with_seed_kwarg(func, args, kwargs, seed)

    import threading

    result_container = [None]
    error_container = [None]

    def _run():
        try:
            result_container[0] = func(*args, **call_kwargs)
        except Exception as e:
            error_container[0] = f"{type(e).__name__}: {e}"

    start = time.time()
    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    thread.join(timeout=timeout)
    elapsed = time.time() - start

    if thread.is_alive():
        # Thread is still running — timeout. Daemon thread will die with process.
        return {
            "success": False,
            "output": None,
            "error": f"Timeout after {timeout}s",
            "time_s": round(elapsed, 4),
            "timeout": True,
        }

    if error_container[0] is not None:
        return {
            "success": False,
            "output": None,
            "error": error_container[0],
            "time_s": round(elapsed, 4),
        }

    return {
        "success": True,
        "output": result_container[0],
        "error": None,
        "time_s": round(elapsed, 4),
    }


def _kill_proc_tree(pid: int):
    """Kill a process and all its children. Works on Windows and Unix."""
    try:
        import signal
        if sys.platform == "win32":
            # taskkill /T kills the entire process tree
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(pid)],
                capture_output=True, timeout=10,
            )
        else:
            os.killpg(os.getpgid(pid), signal.SIGKILL)
    except Exception:
        # Last resort: just kill the process itself
        try:
            os.kill(pid, 9)
        except Exception:
            pass


def execute_function_subprocess(filepath: str, func_name: str, args: list, kwargs: dict,
                                timeout: int = DEFAULT_TIMEOUT,
                                seed: Optional[int] = None) -> dict:
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
import sys, json, traceback, os, inspect
import random
random.seed({seed!r})
try:
    import numpy as np
    np.random.seed({seed!r})
except Exception:
    pass

# Add parent directory to path
sys.path.insert(0, os.path.dirname(r"{filepath}"))

def smart_serialize(obj):
    """Convert non-JSON-serializable objects to informative representations."""
    if obj is None:
        return None
    if isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, (list, tuple)):
        return [smart_serialize(x) for x in obj]
    if isinstance(obj, dict):
        return {{k: smart_serialize(v) for k, v in obj.items()}}
    if isinstance(obj, set):
        return sorted([smart_serialize(x) for x in obj], key=lambda x: str(x))
    if isinstance(obj, bytes):
        return {{"__type__": "bytes", "length": len(obj)}}

    type_name = type(obj).__name__
    module_name = type(obj).__module__ or ""
    full_type = f"{{module_name}}.{{type_name}}" if module_name else type_name

    # PIL Image
    if "PIL" in module_name or type_name == "Image":
        try:
            return {{"__type__": "PIL.Image", "mode": obj.mode, "size": list(obj.size)}}
        except Exception:
            return {{"__type__": "PIL.Image"}}

    # numpy array
    if type_name == "ndarray":
        try:
            return {{"__type__": "numpy.ndarray", "shape": list(obj.shape),
                     "dtype": str(obj.dtype), "values": obj.tolist()}}
        except Exception:
            return {{"__type__": "numpy.ndarray", "shape": list(obj.shape)}}

    # pandas DataFrame
    if type_name == "DataFrame":
        try:
            return {{"__type__": "pandas.DataFrame",
                     "shape": list(obj.shape),
                     "columns": list(obj.columns),
                     "data": obj.head(20).to_dict(orient="records")}}
        except Exception:
            return {{"__type__": "pandas.DataFrame"}}

    # pandas Series
    if type_name == "Series":
        try:
            return {{"__type__": "pandas.Series",
                     "length": len(obj),
                     "values": obj.head(20).tolist()}}
        except Exception:
            return {{"__type__": "pandas.Series"}}

    # RDKit Mol object
    if "rdkit" in module_name.lower() or "Chem" in module_name:
        try:
            from rdkit import Chem
            if hasattr(obj, "GetNumAtoms"):
                smi = Chem.MolToSmiles(obj)
                return {{"__type__": "rdkit.Mol", "smiles": smi,
                         "num_atoms": obj.GetNumAtoms()}}
        except Exception:
            pass
        return {{"__type__": full_type, "str": str(obj)[:200]}}

    # RDKit enum types (BondType, HybridizationType, etc.)
    if "rdkit" in module_name.lower():
        try:
            return {{"__type__": full_type, "name": obj.name if hasattr(obj, "name") else str(obj)}}
        except Exception:
            return {{"__type__": full_type, "str": str(obj)[:200]}}

    # Fallback: try str()
    try:
        return {{"__type__": full_type, "str": str(obj)[:500]}}
    except Exception:
        return {{"__type__": full_type}}

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
if {seed!r} is not None and "seed" not in kwargs:
    try:
        signature = inspect.signature(func)
        parameter = signature.parameters.get("seed")
        bound = signature.bind_partial(*args, **kwargs)
        if parameter is not None and parameter.kind in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        ) and "seed" not in bound.arguments:
            kwargs["seed"] = {seed!r}
    except (TypeError, ValueError):
        pass

try:
    output = func(*args, **kwargs)
    serialized = smart_serialize(output)
    result = {{"success": True, "output": serialized, "error": None, "syntax_error": False, "runtime_error": False}}
except Exception as e:
    tb = traceback.format_exc()
    result = {{"success": False, "output": None, "error": f"{{type(e).__name__}}: {{e}}\\n{{tb}}", "syntax_error": False, "runtime_error": True}}

print("__RESULT__" + json.dumps(result, default=str))
'''

    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False, encoding='utf-8') as f:
        f.write(wrapper_code)
        wrapper_path = f.name

    start = time.time()
    proc = None
    try:
        # Use Popen + CREATE_NEW_PROCESS_GROUP on Windows so we can
        # reliably kill the entire process tree on timeout.
        kwargs_popen = dict(
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=os.path.dirname(filepath),
            env={**os.environ, "PYTHONIOENCODING": "utf-8"},
        )
        if sys.platform == "win32":
            kwargs_popen["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP

        proc = subprocess.Popen(
            [sys.executable, wrapper_path],
            **kwargs_popen,
        )

        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            # Kill the entire process tree
            _kill_proc_tree(proc.pid)
            try:
                stdout, stderr = proc.communicate(timeout=5)
            except Exception:
                stdout, stderr = "", ""
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

        elapsed = time.time() - start

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
        # Make sure the process is dead
        if proc and proc.poll() is None:
            _kill_proc_tree(proc.pid)
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
