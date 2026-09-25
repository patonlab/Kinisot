#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Command-line behaviour: exit codes, the results file, flags."""

import os
import subprocess
import sys

import pytest
from conftest import datapath, write_minimum, write_ts

from kinisot import Kinisot

GS = datapath("gaussian/claisen_gs.out")
TS = datapath("gaussian/claisen_ts.out")
TMCH = datapath("gaussian/tetramethylcyclohexane.out")
CLAISEN = ["--rct", GS, "--ts", TS, "--iso", "5", "-t", "393", "-s", "0.961"]


@pytest.fixture
def run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def _run(args):
        return Kinisot.main(args), tmp_path

    return _run


def read(tmp_path, name="Kinisot_output.dat"):
    with open(os.path.join(str(tmp_path), name)) as handle:
        return handle.read()


def test_kie_run(run, capsys):
    rc, tmp_path = run(CLAISEN)
    assert rc == 0
    out = read(tmp_path)
    assert "Species: %s isotopologue: 5" % GS in out
    assert "KIE @ 393.0 K" in out
    # golden values from the characterization tests, printed to 6 decimals
    assert "1.000176   0.999077   1.000680   1.001962   1.001895   1.000044   1.001940" in out
    # per-species rows carry each species' own Bigeleisen-Mayer factors
    assert "claisen_gs: iso @ 5" in out and "claisen_ts: iso @ 5" in out
    assert "imaginary 463.9i; 35 kept" in out
    assert "Results appended to Kinisot_output.dat" in capsys.readouterr().out


def test_eqe_run_and_auto_scaling(run):
    rc, tmp_path = run(["--rct", TMCH, "--prd", TMCH, "--iso", "24,25,26", "--iso", "28,29,30", "-t", "300"])
    assert rc == 0
    out = read(tmp_path)
    assert "Found vibrational scaling factor 0.977 for RB3LYP/6-31G(d)" in out
    assert "EQE @ 300.0 K" in out


def test_multiple_reactant_files(run):
    rc, tmp_path = run(
        [
            "--rct",
            datapath("gaussian/dienophile.out"),
            "--rct",
            datapath("gaussian/diene.out"),
            "--ts",
            datapath("gaussian/DATS.out"),
            "--iso",
            "0",
            "--iso",
            "6",
            "--iso",
            "15",
            "-s",
            "0.963",
        ]
    )
    assert rc == 0
    out = read(tmp_path)
    assert "dienophile + diene: iso @ 0 / 6" in out
    assert "1.018631" in out and "1.022411" in out


def test_results_append_unless_overwrite(run):
    run(CLAISEN)
    _, tmp_path = run(CLAISEN)
    assert read(tmp_path).count("KINISOT.py v") == 2
    run(CLAISEN + ["--overwrite"])
    assert read(tmp_path).count("KINISOT.py v") == 1


def test_output_path_and_quiet(run, capsys):
    rc, tmp_path = run(CLAISEN + ["-o", "claisen.dat", "--quiet"])
    assert rc == 0
    assert "KIE @ 393.0 K" in read(tmp_path, "claisen.dat")
    assert capsys.readouterr().out == ""


def test_cutoff_alias_still_accepted(run):
    rc, _ = run(CLAISEN + ["--cutoff", "50"])
    assert rc == 0


def test_version(run, capsys):
    with pytest.raises(SystemExit) as exc:
        run(["--version"])
    assert exc.value.code == 0
    assert "Kinisot " + Kinisot.__version__ in capsys.readouterr().out


@pytest.mark.parametrize(
    "args",
    [
        ["--rct", GS, "--iso", "5"],  # neither --ts nor --prd
        ["--rct", GS, "--ts", TS, "--prd", TS, "--iso", "5"],  # both
        ["--rct", GS, "--ts", TS, "--iso", "5", "--iso", "5", "--iso", "5"],  # label count
        ["--rct", GS, "--ts", TS, "--iso", "5", "-t", "0"],
        ["--rct", GS, "--ts", TS, "--iso", "5", "-s", "-1"],
        ["--rct", GS, "--ts", TS, "--iso", "5", "--no-such-flag"],
    ],
)
def test_usage_errors_exit_2(run, args):
    with pytest.raises(SystemExit) as exc:
        run(args)
    assert exc.value.code == 2


def test_input_error_exit_1_and_recorded(run, capsys):
    rc, tmp_path = run(["--rct", GS, "--ts", TS, "--iso", "99"])
    assert rc == 1
    assert "atom 99 is out of range" in capsys.readouterr().err
    assert "ERROR: " in read(tmp_path)


def test_unwritable_output_exit_1(run, capsys):
    rc, _ = run(CLAISEN + ["-o", "no_such_dir/out.dat"])
    assert rc == 1
    assert "cannot open the results file" in capsys.readouterr().err


def test_parse_error_exit_1(run, capsys):
    rc, _ = run(["--rct", "missing.out", "--ts", TS, "--iso", "5"])
    assert rc == 1
    assert "cannot read missing.out" in capsys.readouterr().err


def test_warning_and_level_mismatch_are_printed(run, tmp_path):
    rct = write_minimum(tmp_path / "rct.out", level="RM062X")
    ts = write_ts(tmp_path / "ts.out")
    rc, _ = run(["--rct", rct, "--ts", ts, "--iso", "1"])
    assert rc == 0
    out = read(tmp_path)
    assert "WARNING: the files were not computed at the same level of theory" in out
    assert "RM062X/6-31G(d)" in out and "Vib. scale factor = 1.0" in out


def test_python_m_kinisot(tmp_path):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in [os.path.dirname(os.path.dirname(os.path.abspath(__file__))), env.get("PYTHONPATH", "")] if p
    )
    proc = subprocess.run(
        [sys.executable, "-m", "kinisot"] + CLAISEN, cwd=str(tmp_path), env=env, capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr
    assert "KIE @ 393.0 K" in proc.stdout
    assert os.path.exists(os.path.join(str(tmp_path), "Kinisot_output.dat"))
