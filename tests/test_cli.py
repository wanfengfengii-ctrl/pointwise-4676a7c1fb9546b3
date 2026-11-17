import json
import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


class TestCLI(unittest.TestCase):
    def test_persistent_public_roundtrip(self):
        root = Path(__file__).resolve().parents[1]
        with TemporaryDirectory(prefix="remap-cli-") as directory:
            db = str(Path(directory) / "cli.sqlite")

            def call(*args, status=0):
                completed = subprocess.run(
                    [sys.executable, "-B", "-m", "remap", "--db", db, *args],
                    cwd=root, text=True, capture_output=True, timeout=30)
                self.assertEqual(completed.returncode, status, completed.stderr)
                return json.loads(completed.stdout)

            revision = call("source", str(root / "samples/source.json"))["source_revision"]
            self.assertEqual(call("head", "unit-observation")["source_revision"], revision)
            data = json.loads((root / "samples/plan.json").read_text(encoding="utf-8"))
            data["source_revision"] = revision
            plan_file = Path(directory) / "plan.json"
            plan_file.write_text(json.dumps(data), encoding="utf-8")
            plan = call("plan", str(plan_file))["plan_id"]
            result = call("execute", str(plan), "cli-request")
            self.assertEqual(result["transfers"][0]["mass"], "11/384")
            self.assertEqual(result["source_ledger"][0]["uncovered_mass"], "5/384")
            self.assertEqual(call("result", str(result["id"])), result)
            self.assertEqual(call("execute", str(plan), "cli-request"), result)
            self.assertEqual(len(call("catalog")["results"]), 1)
            error = call("result", "999", status=2)
            self.assertIn("error", error)
            bad = Path(directory) / "bad.json"
            bad.write_text('{"name":"bad","cells":[]}', encoding="utf-8")
            self.assertIn("error", call("source", str(bad), status=2))
            self.assertEqual(len(call("catalog")["sources"]), 1)
