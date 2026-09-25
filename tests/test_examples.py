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
    # expand the Claisen loop by hand; every other command is a plain `run` line
    loop = re.search(r"for atoms in (.*?); do\n\s*run claisen (.*?)\n", script)
    commands = [("claisen", "gaussian", loop.group(2).replace('"$atoms"', atoms)) for atoms in loop.group(1).split()]
    directory = "gaussian"
    for line in script.splitlines():
        line = line.strip()
        if line.startswith("cd "):
            directory = "orca" if "orca" in line else "ase" if "ase" in line else "gaussian"
        elif line.startswith("run ") and '"$atoms"' not in line:
            case, args = line.split(None, 1)[1].split(None, 1)
            commands.append((case, directory, args))
    return commands


COMMANDS = example_commands()


def test_all_example_commands_found():
    assert len(COMMANDS) == 28
    assert sum(1 for case, _, _ in COMMANDS if case == "claisen") == 12


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
