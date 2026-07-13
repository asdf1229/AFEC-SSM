#!/usr/bin/env python3
"""Check DecQ's lattice/output structure on the five-edge diamond fixture."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


EXPECTED = {
    "Local Patterns": 14,
    "Global Patterns": 14,
    "Minimal Patterns": 8,
    "Global-Pattern Rows": 18,
    "Unique Results": 2,
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args()
    if not args.executable.is_file():
        parser.error(f"missing executable: {args.executable}")

    root = Path(__file__).resolve().parent / "sasum_oracle" / "paper_example"
    completed = subprocess.run(
        [
            str(args.executable.resolve()),
            "-d",
            str(root / "graph_g.txt"),
            "-q",
            str(root / "query_graph" / "q.txt"),
            "-t",
            "2",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=args.timeout,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{args.executable} exited with {completed.returncode}:\n{completed.stdout}"
        )
    if "(reached)" in completed.stdout:
        raise RuntimeError("DecQ reached a configured resource limit")

    failures: list[str] = []
    for label, expected in EXPECTED.items():
        match = re.search(rf"^{re.escape(label)}:\s+(\d+)\s*$", completed.stdout, re.MULTILINE)
        observed = int(match.group(1)) if match else None
        if observed != expected:
            failures.append(f"{label}: observed {observed}, expected {expected}")
    if failures:
        raise RuntimeError("DecQ structural check failed:\n- " + "\n- ".join(failures))

    print("PASS: DecQ lattice=14, minimal=8, pattern rows=18, unique mappings=2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
