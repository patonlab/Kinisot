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
RESULT = re.compile(r"^  (KIE|EQE)( \((ensemble|series|channels)\))? @ .*$", re.M)


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
            directory = next((d for d in ("orca", "xtb", "ase", "conformers") if "/" + d in line), "gaussian")
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
    assert len(COMMANDS) == 35
    assert sum(1 for case, _, _ in COMMANDS if case == "claisen") == 12
    assert sum(1 for _, directory, _ in COMMANDS if directory == "xtb") == 4
    assert sum(1 for _, directory, _ in COMMANDS if directory == "conformers") == 3


@pytest.mark.parametrize("case", sorted({case for case, _, _ in COMMANDS}))
def test_examples_reproduce_expected_output(case, tmp_path, monkeypatch):
    output = str(tmp_path / "output.dat")
    for command_case, directory, args in COMMANDS:
        if command_case == case:
            monkeypatch.chdir(os.path.join(EXAMPLES, directory) if directory == "conformers" else datapath(directory))
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
    # Meyer, DelMonte & Singleton 1999, Table 4: relative to C5, two experiments compared through their mean
    reference = module.compute_kie(
        rct=module.resolve(cases[0], "reactants"),
        ts=module.resolve(cases[0], "transition_structure"),
        iso="4",
        temperature=393.0,
        scale=0.961,
        reference="5",
    )
    assert rows[3]["position"] == "C4" and rows[3]["experimental"] == [1.035, 1.033]
    assert rows[3]["computed"] == pytest.approx(reference.kie_tunnel_relative, abs=1e-9)
    assert rows[3]["computed"] == pytest.approx(1.0310, abs=1e-4)
    assert rows[3]["deviation"] == pytest.approx(reference.kie_tunnel_relative - 1.034, abs=1e-9)
    measured = [r for r in rows if r["deviation"] is not None]
    assert len(measured) == 5 and sum(abs(r["deviation"]) for r in measured) / 5 == pytest.approx(0.0009, abs=1e-4)
    text = module.format_case(cases[0], rows)
    assert "1.035 ± 0.002, 1.033 ± 0.002" in text and "Mean absolute deviation over 5 measured positions" in text
    assert "no experimental value" in text  # H7,H8: no 2H KIE was measured at 120 C
    # a per-entry reference overrides the case's
    cases[0]["kies"][3]["reference"] = "6"
    assert module.run_case(cases[0])[3]["computed"] != pytest.approx(rows[3]["computed"], abs=1e-4)
    # the Diels-Alder case against Singleton & Thomas 1995 (relative to the methyl group; 2H to the mean of its Hs)
    (case,) = module.load_cases(["diels_alder"])
    rows = {r["position"]: r for r in module.run_case(case)}
    assert rows["C1"]["computed"] == pytest.approx(1.0216, abs=1e-4) and rows["C1"]["experimental"] == 1.022
    assert rows["C4"]["computed"] == pytest.approx(1.0172, abs=1e-4)
    assert rows["H1Z (inside)"]["deviation"] == pytest.approx(0.018, abs=0.001)
    deviations = [abs(r["deviation"]) for r in rows.values()]
    assert len(deviations) == 9 and sum(deviations) / 9 == pytest.approx(0.0031, abs=1e-4)
    # Shi epoxidation (Singleton & Wang 2005): Gaussian jobs at the SI geometries reproduce the paper's six QUIVER
    # predictions for TS 10 to the three decimals given (the note column), relative to the mean of the meta carbons
    (case,) = module.load_cases(["shi_epoxidation"])
    rows = module.run_case(case)
    for row in rows:
        assert round(row["computed"], 3) == float(row["note"].split()[1]), row["position"]
    assert rows[0]["computed"] == pytest.approx(1.0221, abs=1e-4) and rows[0]["experimental"] == [1.022, 1.02]
    assert sum(abs(r["deviation"]) for r in rows) / 6 == pytest.approx(0.0012, abs=1e-4)
    # cases from the PyQuiverHS SI: an EQE (product file) and a small reaction-coordinate frequency (imag_cutoff)
    (case,) = module.load_cases(["tetramethylcyclohexane_eie"])
    (row,) = module.run_case(case)
    assert row["computed"] == pytest.approx(1.0417, abs=1e-4) and row["deviation"] == pytest.approx(-0.0003, abs=1e-4)
    (case,) = module.load_cases(["biaryl_diketone"])
    assert module.run_case(case)[0]["computed"] == pytest.approx(1.0752, abs=1e-4)
    del case["imag_cutoff"]  # the 45i cm-1 mode is not a reaction coordinate at the default 50 cm-1 cutoff
    with pytest.raises(Exception, match="imaginary frequency beyond the 50.0 cm-1 cutoff"):
        module.run_case(case)
    # a case with measurements but no structures yet (Baeyer-Villiger) is listed, not computed
    (case,) = module.load_cases(["baeyer_villiger"])
    rows = module.run_case(case)
    assert all(r["computed"] is None and r["deviation"] is None for r in rows)
    text = module.format_case(case, rows)
    assert "Not computed yet" in text and "1.0096 ± 0.0006" in text and "1.001 ± 0.002, 1.000 ± 0.002" in text
    # a partly filled case (transition structure but no reactants, or no iso labels) is still only listed
    case["transition_structure"] = ["ts.out"]
    assert all(r["computed"] is None for r in module.run_case(case))
    assert "Not computed yet" in module.format_case(case, rows)
    (case,) = module.load_cases(["claisen"])
    case["kies"][0]["iso"] = None
    rows = module.run_case(case)
    assert rows[0]["computed"] is None and rows[1]["computed"] is not None
    del sys, json
