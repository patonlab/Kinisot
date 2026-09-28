#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""JSON job files (kinisot --job; kinisot/jobs.py): the command line gives what compute_kie and channels give."""

import json
import os
import shutil
import warnings

import pytest

from kinisot import Conformers, KinisotInputError, KinisotWarning, Series, channels, cli, compute_kie
from kinisot.jobs import load_job

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SHI = os.path.join(ROOT, "benchmarks", "shi_epoxidation")
NITRO = os.path.join(ROOT, "benchmarks", "nitroarene_phosphetane")
REACTANT = os.path.join(SHI, "methylstyrene.log")
TS = {"A": os.path.join(SHI, "ts10.log"), "B": os.path.join(SHI, "ensemble", "ts_B.log"),
      "G": os.path.join(SHI, "ensemble", "ts_G.log")}  # fmt: skip
SETTINGS = {"temperature": 273.15, "scale": 0.9614, "project": True}
OPTIONS = dict(temperature=273.15, scale=0.9614, project=True)


def quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)
        return function(*args, **kwargs)


def run(tmp_path, monkeypatch, job, *extra):
    """Write the job next to the working directory, with file names relative to it, and run it."""
    monkeypatch.chdir(tmp_path)
    with open("job.json", "w") as handle:
        json.dump(job, handle)
    code = quiet(cli.main, ["--job", "job.json", "-q", "--json", "out.json", *extra])
    if code != 0:
        return code, None
    with open("out.json") as handle:
        return code, json.load(handle)


def rel(tmp_path, path):
    """``path`` relative to the job's directory; absolute when there is no relative path (Windows, another drive)."""
    try:
        return os.path.relpath(path, str(tmp_path))
    except ValueError:
        return path


def test_file_names_are_relative_to_the_job_file(tmp_path, monkeypatch):
    """Files next to the job (in a subdirectory), run from another directory: the same on every platform."""
    data = tmp_path / "job" / "data"
    data.mkdir(parents=True)
    for name in ("gs_1.hessian.json", "ts_chair.hessian.json"):
        shutil.copy(os.path.join(ROOT, "examples", "conformers", name), str(data / name))
    job = {"temperature": 393, "scale": 1.0, "reactants": ["data/gs_1.hessian.json"],
           "transition_structure": ["data/ts_chair.hessian.json"], "iso": "4"}  # fmt: skip
    with open(str(tmp_path / "job" / "job.json"), "w") as handle:
        json.dump(job, handle)
    monkeypatch.chdir(tmp_path)
    assert quiet(cli.main, ["--job", os.path.join("job", "job.json"), "-q", "--json", "out.json"]) == 0
    with open("out.json") as handle:
        result = json.load(handle)
    expected = quiet(compute_kie, rct=str(data / "gs_1.hessian.json"), ts=str(data / "ts_chair.hessian.json"), iso="4",
                     temperature=393, scale=1.0)  # fmt: skip
    assert result["kie_tunnel"] == pytest.approx(expected.kie_tunnel, abs=1e-13)


def test_series_job(tmp_path, monkeypatch):
    job = dict(
        SETTINGS,
        reactants=[rel(tmp_path, REACTANT)],
        series={"steps": [rel(tmp_path, TS["A"]), [rel(tmp_path, TS["B"]), rel(tmp_path, TS["G"])]],
                "free_energies": [0.0, 0.4]},
        isotopologues=[
            {"name": "C-beta", "iso": ["2", "8", "8"], "reference": ["11", "4", "4"]},
            {"name": "C-alpha", "iso": ["1", "7", "7"]},
        ],
    )  # fmt: skip
    code, data = run(tmp_path, monkeypatch, job, "--csv", "rows.csv")
    assert code == 0 and [d["isotopologue"] for d in data] == ["C-beta", "C-alpha"]
    series = Series([TS["A"], [TS["B"], TS["G"]]], free_energies=[0.0, 0.4])
    expected = quiet(compute_kie, rct=REACTANT, ts=series, iso=["2", "8", "8"], reference=["11", "4", "4"], **OPTIONS)
    assert data[0]["kie_tunnel"] == pytest.approx(expected.kie_tunnel, abs=1e-13)
    assert data[0]["kie_tunnel_relative"] == pytest.approx(expected.kie_tunnel_relative, abs=1e-13)
    assert data[1]["reference"] is None and data[0]["series"] is True
    with open("Kinisot_output.dat") as handle:
        text = handle.read()
    assert "Isotopologue: C-alpha" in text and "KIE (series) @ 273.15 K" in text and "Steps in series" in text
    assert "Reference atoms (step 1): methylstyrene C11 -> 13C; ts10 C4 -> 13C" in text
    assert "Labelled atoms (step 2): methylstyrene C1 -> 13C; ts_B (+1 conformer) C7 -> 13C" in text
    with open("rows.csv") as handle:
        assert handle.readline().startswith("isotopologue,kind")


def test_channels_job(tmp_path, monkeypatch):
    phno2, ts1b = (rel(tmp_path, os.path.join(NITRO, f)) for f in ("phno2.out", "ts1b.out"))
    job = {
        "temperature": 393, "scale": 0.9614, "project": True, "tunneling": "none",
        "channels": [
            {"name": "attacked", "reactants": [phno2], "transition_structure": [ts1b]},
            {"name": "spectator", "reactants": [phno2], "transition_structure": [ts1b], "amount": 1},
        ],
        "shares": [1, 1],
        "iso": [["12:15N,13:18O", "1:15N,2:18O"], ["12:15N,14:18O", "1:15N,3:18O"]],
        "reference": [["12:15N", "1:15N"], ["12:15N", "1:15N"]],
    }  # fmt: skip
    code, data = run(tmp_path, monkeypatch, job)
    assert code == 0 and data["channels"] is True and data["names"] == ["attacked", "spectator"]
    assert data["kie_tunnel_relative"] == pytest.approx(1.0304, abs=1e-4)
    with open("Kinisot_output.dat") as handle:
        text = handle.read()
    assert "KIE (channels) @ 393.0 K" in text and "o spectator: ts1b" in text
    assert "Labelled atoms (spectator): phno2 N12 -> 15N, O14 -> 18O; ts1b N1 -> 15N, O3 -> 18O" in text


def test_conformer_and_eqe_jobs(tmp_path, monkeypatch):
    job = dict(
        SETTINGS,
        reactants=[rel(tmp_path, REACTANT)],
        transition_structure=[{"files": [rel(tmp_path, TS["A"]), rel(tmp_path, TS["B"])], "free_energies": [0, 0.3],
                               "degeneracy": [1, 2]}],
        iso=["2", "8"], weights="user", temperature=[250, 300],
    )  # fmt: skip
    code, data = run(tmp_path, monkeypatch, job)
    assert code == 0 and [d["temperature"] for d in data] == [250.0, 300.0]
    expected = quiet(compute_kie, rct=REACTANT, ts=Conformers([TS["A"], TS["B"]], [0, 0.3], [1, 2]), iso=["2", "8"],
                     weights="user", **dict(OPTIONS, temperature=300.0))  # fmt: skip
    assert data[1]["kie_tunnel"] == pytest.approx(expected.kie_tunnel, abs=1e-13) and data[1]["ensemble"]

    tmch = rel(tmp_path, os.path.join(ROOT, "benchmarks", "tetramethylcyclohexane_eie", "tetramethylcyclohexane.log"))
    job = {"temperature": "290.15", "scale": 1.0, "reactants": [tmch], "product": [tmch],
           "iso": ["28,29,30", "24,25,26"]}  # fmt: skip
    code, data = run(tmp_path, monkeypatch, job)
    assert code == 0 and data["kind"] == "EQE" and data["kie"] == pytest.approx(1.042, abs=2e-3)


def test_job_errors(tmp_path, monkeypatch):
    base = dict(SETTINGS, reactants=[rel(tmp_path, REACTANT)], transition_structure=[rel(tmp_path, TS["A"])],
                iso=["2", "8"])  # fmt: skip
    for broken, message in (
        (dict(base, temprature=300), "unknown key"),
        (dict(base, product=[rel(tmp_path, TS["A"])]), "exactly one of"),
        ({k: v for k, v in base.items() if k != "iso"}, "\"iso\""),
        (dict(base, shares=[1, 2]), "belongs to a job with \"channels\""),
        (dict(base, isotopologues=[{"iso": ["2", "8"]}]), "one of them"),
        (dict(base, temperature=-5), "positive"),
        (dict(base, transition_structure=[{"file": "x"}]), "unknown key"),
        ({"channels": [{"reactants": ["a"], "transition_structure": ["b"]}], "iso": [["1", "2"], ["1", "2"]]},
         "one entry"),
        ({"channels": [{"reactants": ["a"]}], "iso": [["1"]]}, "one of \"transition_structure\""),
    ):  # fmt: skip
        with open(os.path.join(str(tmp_path), "bad.json"), "w") as handle:
            json.dump(broken, handle)
        with pytest.raises(KinisotInputError, match=message):
            load_job(os.path.join(str(tmp_path), "bad.json"))
    with open(os.path.join(str(tmp_path), "bad.json"), "w") as handle:
        handle.write("{not json")
    with pytest.raises(KinisotInputError, match="not valid JSON"):
        load_job(os.path.join(str(tmp_path), "bad.json"))
    with pytest.raises(KinisotInputError, match="cannot read"):
        load_job(os.path.join(str(tmp_path), "missing.json"))
    # settings belong in the job file, and an invalid job is a usage error
    monkeypatch.chdir(tmp_path)
    with open("job.json", "w") as handle:
        json.dump(base, handle)
    for extra in (["-t", "300"], ["--iso", "2"], ["--tunneling", "none"], ["--project"]):
        with pytest.raises(SystemExit):
            cli.main(["--job", "job.json", "-q", *extra])
    with pytest.raises(SystemExit):
        cli.main(["--job", "bad.json", "-q"])
    with pytest.raises(SystemExit):
        cli.main(["-q"])  # neither --rct nor --job
    # a chemistry error is reported with exit code 1
    with open("job.json", "w") as handle:
        json.dump(dict(base, iso=["2", "99"]), handle)
    assert quiet(cli.main, ["--job", "job.json", "-q"]) == 1


def test_job_matches_channels_api(tmp_path, monkeypatch):
    job = dict(
        SETTINGS,
        channels=[
            {"name": "A", "reactants": [rel(tmp_path, REACTANT)], "transition_structure": [rel(tmp_path, TS["A"])]},
            {"name": "G", "reactants": [rel(tmp_path, REACTANT)],
             "series": {"steps": [rel(tmp_path, TS["G"]), rel(tmp_path, TS["B"])], "commitment": 2.0}},
        ],
        barriers=[0.0, 1.5],
        iso=[["2", "8"], ["2", "8", "8"]],
    )  # fmt: skip
    code, data = run(tmp_path, monkeypatch, job)
    assert code == 0
    expected = quiet(channels, [
        dict(name="A", rct=[REACTANT], ts=[TS["A"]], iso=["2", "8"]),
        dict(name="G", rct=[REACTANT], ts=Series([TS["G"], TS["B"]], commitment=2.0), iso=["2", "8", "8"]),
    ], barriers=[0.0, 1.5], **OPTIONS)  # fmt: skip
    assert data["kie_tunnel"] == pytest.approx(expected.kie_tunnel, abs=1e-13)
    assert data["shares"] == pytest.approx(list(expected.shares), abs=1e-13)
    assert data["results"][1]["series"] is True


def test_isotopologues_must_be_a_non_empty_list(tmp_path):
    base = dict(SETTINGS, reactants=[rel(tmp_path, REACTANT)], transition_structure=[rel(tmp_path, TS["A"])])
    for value in ([], {}, "2"):
        with open(os.path.join(str(tmp_path), "job.json"), "w") as handle:
            json.dump(dict(base, isotopologues=value), handle)
        with pytest.raises(KinisotInputError, match="non-empty list"):
            load_job(os.path.join(str(tmp_path), "job.json"))
