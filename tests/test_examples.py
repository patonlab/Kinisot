#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The worked examples must reproduce their committed expected outputs.

Replays every `run <case> ...` line of examples/run_examples.sh in-process
and compares the result lines with examples/<case>/expected_output.dat, so
the README and example write-ups cannot drift from what Kinisot prints.
"""

import os
import re
import shlex

import pytest
from conftest import datapath

from kinisot import cli

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")
RESULT = re.compile(r"^  (KIE|EQE) @ .*$", re.M)


def example_commands():
    with open(os.path.join(EXAMPLES, "run_examples.sh")) as handle:
        script = handle.read()
    # walk the script, expanding `for atoms in ...` loops and tracking the directory of each run
    commands = []
    directory = "gaussian"
    loop_values = None
    for line in script.splitlines():
        line = line.strip()
        if line.startswith("cd "):
            directory = next((d for d in ("orca", "xtb", "ase") if "/" + d in line), "gaussian")
        elif line.startswith("for atoms in "):
            loop_values = line[len("for atoms in ") :].split(";")[0].split()
        elif line == "done":
            loop_values = None
        elif line.startswith("run "):
            case, args = line.split(None, 1)[1].split(None, 1)
            for value in loop_values if '"$atoms"' in args else [None]:
                commands.append((case, directory, args.replace('"$atoms"', value) if value else args))
    return commands


COMMANDS = example_commands()


def test_all_example_commands_found():
    assert len(COMMANDS) == 32
    assert sum(1 for case, _, _ in COMMANDS if case == "claisen") == 12
    assert sum(1 for _, directory, _ in COMMANDS if directory == "xtb") == 4


@pytest.mark.parametrize("case", sorted({case for case, _, _ in COMMANDS}))
def test_examples_reproduce_expected_output(case, tmp_path, monkeypatch):
    output = str(tmp_path / "output.dat")
    for command_case, directory, args in COMMANDS:
        if command_case == case:
            monkeypatch.chdir(datapath(directory))
            assert cli.main(shlex.split(args) + ["--quiet", "--output", output]) == 0
    with open(os.path.join(EXAMPLES, case, "expected_output.dat")) as handle:
        expected_text = handle.read()
    with open(output) as handle:
        produced_text = handle.read()
    produced = [m.group(0) for m in RESULT.finditer(produced_text)]
    expected = [m.group(0) for m in RESULT.finditer(expected_text)]
    assert produced == expected
    # the species header lines are recorded too
    assert produced_text.count("Species:") == expected_text.count("Species:")


def test_api_example_script_runs():
    import subprocess
    import sys

    proc = subprocess.run([sys.executable, os.path.join(EXAMPLES, "api_example.py")], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert "4            1.0127     1.0366     0.9790     1.0297     1.0330" in proc.stdout


def test_benchmark_runner(tmp_path, monkeypatch):
    import importlib.util
    import json
    import sys

    spec = importlib.util.spec_from_file_location("bench_run", os.path.join(ROOT, "benchmarks", "run.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cases = module.load_cases(["claisen"])
    assert len(cases) == 1 and cases[0]["reference"]["doi"] == "10.1021/ja992372h"
    rows = module.run_case(cases[0])
    assert [r["position"] for r in rows][:2] == ["C1", "C2"]
    assert rows[3]["computed"] == pytest.approx(1.032993, abs=1e-6)  # C4, Bell-corrected, from the Claisen example
    assert all(r["experimental"] is None for r in rows)  # placeholders until entered from the paper
    text = module.format_case(cases[0], rows)
    assert "no experimental value" in text and "enter them in `case.json`" in text
    # a filled-in value produces a deviation and a mean absolute deviation
    cases[0]["kies"][3]["experimental"], cases[0]["kies"][3]["uncertainty"] = 1.030, 0.002
    rows = module.run_case(cases[0])
    assert rows[3]["deviation"] == pytest.approx(1.032993 - 1.030, abs=1e-5)
    assert "Mean absolute deviation over 1 measured positions" in module.format_case(cases[0], rows)
    del sys, json
