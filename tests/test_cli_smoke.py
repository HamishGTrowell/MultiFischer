"""End-to-end smoke test using the distributed example dataset."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class CLISmokeTests(unittest.TestCase):
    def test_single_pair_example_run(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "results"
            environment = os.environ.copy()
            environment["MPLBACKEND"] = "Agg"
            environment["PYTHONUTF8"] = "1"

            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "multifischer",
                    "--input",
                    str(REPOSITORY_ROOT / "example.csv"),
                    "--out",
                    str(output),
                    "--singlepair",
                    "365",
                    "385",
                ],
                cwd=REPOSITORY_ROOT,
                env=environment,
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )

            self.assertEqual(
                result.returncode,
                0,
                msg=f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}",
            )
            for filename in (
                "365_385nm_qy_ratio_sensitivity.csv",
                "365_385nm_qy_ratio_sensitivity.png",
                "run_settings.json",
                "multifischer.log",
            ):
                self.assertTrue((output / filename).is_file(), filename)


if __name__ == "__main__":
    unittest.main()
