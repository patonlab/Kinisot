#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The Wittig case (benchmarks/wittig_anisaldehyde): Chen, Nieves-Quinones, Waas, Singleton, JACS 2014.

Gaussian 16 jobs at the SI geometries must reproduce the SI energies, and
Kinisot must reproduce the SI Table S4 KIEs for 4-TS and 6-TS (M06-2X/6-31+G**/PCM,
67 C, scaling factor 0.9614, Bell tunnelling) and, with the two in series, the
weighted predictions of Table 1.
"""

import importlib.util
import os
import warnings

import pytest

from kinisot import KinisotWarning, compute_kie, parse_gaussian
from kinisot.jobs import load_job, run_job

CASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks", "wittig_anisaldehyde")
T, SCALE = 340.15, 0.9614
# SI electronic energies (hartree)
SI = {"anisaldehyde_1": -459.931869643, "ylide_2": -1227.87396966, "ts_4": -1687.79361546, "ts_6": -1687.79472605}
# position: labels (anisaldehyde, ylide, TS; 4-TS and 6-TS share their numbering) and SI Table S4
# (M06-2X/6-31+G**) for 4-TS and 6-TS
TABLE_S4 = {
    "ylide CH carbon": (["0", "1", "4"], 1.022, 0.994),
    "carbonyl CHO": (["9", "0", "2"], 1.043, 1.015),
    "ipso": (["4", "0", "5"], 0.999, 0.999),
    "ortho": (["5", "0", "24"], 1.001, 1.001),
    "ortho'": (["3", "0", "28"], 1.001, 1.001),
    "meta": (["6", "0", "25"], 1.000, 1.000),
    "meta'": (["2", "0", "27"], 1.000, 1.001),
    "para": (["1", "0", "26"], 1.000, 1.001),
    "ketone carbonyl": (["0", "37", "51"], 0.998, 0.999),
    "ketone methyl": (["0", "38", "52"], 1.000, 1.002),
    "methoxy": (["13", "0", "57"], 1.000, 1.000),
    "aldehyde 18O": (["17:18O", "0", "1:18O"], 1.016, 1.045),
}
# Table 1: carbonyl carbon and ylide CH carbon, weighted by the trajectories and by the free energies
TABLE_1 = {"series.json": (1.033, 1.012), "series_free_energies.json": (1.028, 1.008)}


@pytest.fixture(scope="module")
def files():
    return {name: parse_gaussian(os.path.join(CASE, name + ".log")) for name in SI}


def test_jobs_reproduce_the_si(files):
    for name, energy in SI.items():
        assert files[name].energy == pytest.approx(energy, abs=1e-4), name
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)  # the ylide's 9i cm-1 phenyl torsion
        kie = compute_kie(rct=[files["anisaldehyde_1"], files["ylide_2"]], ts=files["ts_4"], iso=["9", "0", "2"],
                          temperature=T, scale=SCALE, project=True)  # fmt: skip
    assert kie.other.light.imaginary == pytest.approx(285.4 * SCALE, abs=0.5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)
        kie = compute_kie(rct=[files["anisaldehyde_1"], files["ylide_2"]], ts=files["ts_6"], iso=["9", "0", "2"],
                          temperature=T, scale=SCALE, project=True)  # fmt: skip
    assert kie.other.light.imaginary == pytest.approx(110.7 * SCALE, abs=0.5)


@pytest.mark.parametrize("ts, column", [("ts_4", 1), ("ts_6", 2)])
def test_transition_structures_reproduce_si_table_s4(files, ts, column):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)  # the ylide's 9i cm-1 phenyl torsion
        for position, row in TABLE_S4.items():
            result = compute_kie(rct=[files["anisaldehyde_1"], files["ylide_2"]], ts=files[ts], iso=row[0],
                                 temperature=T, scale=SCALE, project=True)  # fmt: skip
            # Table S4 prints three decimals; several values lie on a rounding boundary
            assert result.kie_tunnel == pytest.approx(row[column], abs=6e-4), (ts, position)


@pytest.mark.parametrize("job", sorted(TABLE_1))
def test_series_reproduces_table_1(job):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)
        results = dict(run_job(load_job(os.path.join(CASE, job)), T))
    carbonyl, ylide = results["carbonyl carbon"], results["ylide CH carbon"]
    # absolute KIEs against Table 1; its 1.008 carries the rounding of its single-structure inputs (they give 1.0070)
    assert carbonyl.kie_tunnel == pytest.approx(TABLE_1[job][0], abs=6e-4)
    assert ylide.kie_tunnel == pytest.approx(TABLE_1[job][1], abs=1.2e-3)
    # the series formula with the steps' own KIEs: KIE = (KIE_6 + C_f KIE_4) / (1 + C_f)
    kie_4, kie_6 = (step.kie_tunnel for step in carbonyl.steps)
    c_f = carbonyl.commitment
    assert carbonyl.kie_tunnel == pytest.approx((kie_6 + c_f * kie_4) / (1 + c_f), rel=1e-12)
    assert c_f == pytest.approx(128 / 76 if job == "series.json" else 0.822, abs=1e-3)


def test_benchmark_runner_computes_both_weightings():
    spec = importlib.util.spec_from_file_location("bench_run", os.path.join(CASE, "..", "run.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    trajectories, statistical = module.load_cases(["wittig_anisaldehyde", "wittig_anisaldehyde_statistical"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)
        rows = module.run_case(trajectories)
        alternative = module.run_case(statistical)
    # relative to the experiment's standards; the trajectory weighting matches all three measurements
    assert [round(r["computed"], 4) for r in rows] == [1.0326, 1.0111, 0.9979]
    assert max(abs(r["deviation"]) for r in rows) < 1e-3
    assert [round(r["computed"], 4) for r in alternative] == [1.0277, 1.0059, 0.9975]
    assert statistical["alternative_to"] == "wittig_anisaldehyde"
