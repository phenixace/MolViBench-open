import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from evaluate.evaluator import EvalConfig, evaluate_single


class EvaluatorFallbackTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        self.solution = self.root / "solution.py"
        self.prediction = self.root / "prediction.py"
        self.prediction.write_text(
            "def level_function():\n    return None\n",
            encoding="utf-8",
        )
        self.config = EvalConfig(str(self.root))
        self.config.safe_mode = False

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_solution(self, body="return None", prefix=""):
        self.solution.write_text(
            f"{prefix}def level_function():\n    {body}\n",
            encoding="utf-8",
        )

    def _execute(self, solution_output, prediction_output, prediction_ok=True):
        def fake_execute(path, *_args, **_kwargs):
            if Path(path) == self.solution:
                return {
                    "success": True,
                    "output": solution_output,
                    "time_s": 0.01,
                }
            return {
                "success": prediction_ok,
                "output": prediction_output if prediction_ok else None,
                "error": "candidate failed" if not prediction_ok else "",
                "time_s": 0.01,
            }

        return fake_execute

    @patch("evaluate.evaluator.api_fallback_check")
    @patch("evaluate.evaluator.infer_test_inputs", return_value=[([], {})])
    @patch("evaluate.evaluator.execute_function_direct")
    def test_stochastic_value_mismatch_uses_structure_not_api(
        self,
        execute,
        _inputs,
        api_check,
    ):
        self._write_solution(
            body="return {'score': float(np.random.rand())}",
            prefix="import numpy as np\n",
        )
        execute.side_effect = self._execute(
            {"score": 0.1},
            {"score": 0.9},
        )

        result = evaluate_single(
            str(self.solution),
            str(self.prediction),
            self.config,
        )

        self.assertTrue(result["nondeterministic"])
        self.assertTrue(result["main_pass"])
        self.assertEqual(result["details"][0]["method"], "structural")
        self.assertFalse(result["api_fallback"]["triggered"])
        api_check.assert_not_called()

    @patch("evaluate.evaluator.api_fallback_check")
    @patch("evaluate.evaluator.infer_test_inputs", return_value=[([], {})])
    @patch("evaluate.evaluator.execute_function_direct")
    def test_non_executable_prediction_cannot_be_rescued(
        self,
        execute,
        _inputs,
        api_check,
    ):
        self._write_solution(body="return 1")
        execute.side_effect = self._execute(1, None, prediction_ok=False)

        result = evaluate_single(
            str(self.solution),
            str(self.prediction),
            self.config,
        )

        self.assertFalse(result["executable"])
        self.assertFalse(result["main_pass"])
        self.assertFalse(result["api_fallback"]["triggered"])
        self.assertFalse(result["api_fallback_broad"]["triggered"])
        api_check.assert_not_called()

    @patch("evaluate.evaluator.api_fallback_check", return_value=(True, "ok"))
    @patch("evaluate.evaluator.infer_test_inputs", return_value=[([], {})])
    @patch("evaluate.evaluator.execute_function_direct")
    def test_comparable_wrong_value_is_diagnostic_only(
        self,
        execute,
        _inputs,
        api_check,
    ):
        self._write_solution(body="return 1")
        execute.side_effect = self._execute(1, 2)

        result = evaluate_single(
            str(self.solution),
            str(self.prediction),
            self.config,
        )

        self.assertFalse(result["main_pass"])
        self.assertFalse(result["api_fallback"]["triggered"])
        self.assertTrue(result["diagnostic_fallback_pass"])
        self.assertEqual(api_check.call_count, 1)
        self.assertEqual(
            api_check.call_args.kwargs["min_overlap_ratio"],
            0.7,
        )

    @patch("evaluate.evaluator.api_fallback_check", return_value=(True, "ok"))
    @patch("evaluate.evaluator.infer_test_inputs", return_value=[([], {})])
    @patch("evaluate.evaluator.execute_function_direct")
    def test_non_comparable_executable_output_can_enter_strict_fallback(
        self,
        execute,
        _inputs,
        api_check,
    ):
        self._write_solution(body="return 'CCO'")
        execute.side_effect = self._execute(
            "CCO",
            {"__type__": "PIL.Image", "size": [100, 100]},
        )

        result = evaluate_single(
            str(self.solution),
            str(self.prediction),
            self.config,
        )

        self.assertTrue(result["api_fallback"]["triggered"])
        self.assertTrue(result["main_pass"])
        thresholds = [
            call.kwargs["min_overlap_ratio"] for call in api_check.call_args_list
        ]
        self.assertEqual(thresholds, [0.5, 0.7])


if __name__ == "__main__":
    unittest.main()
