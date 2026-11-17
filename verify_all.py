"""Single verify entry point: original suite plus the independent audit."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(argv):
    return subprocess.run([sys.executable, "-B", *argv], cwd=ROOT).returncode


def main():
    status = run(["-m", "unittest", "discover", "-s", "tests", "-v"])
    if status:
        return status
    return run(["verify_flux.py", "--check"])


if __name__ == "__main__":
    raise SystemExit(main())
