#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The DYKAT case (benchmarks/dykat_allyl_arylation): van Dijk et al., Nat. Catal. 2021.

Gaussian 16 jobs at the paper's geometries must reproduce SI Table 20, and
Kinisot must reproduce the per-enantiomer KIEs of SI Tables 24 and 25
(wB97X-D/6-31G(d), LANL2DZ(f) on Rh, 313.15 K, unscaled, Bell tunnelling,
relative to the reference carbon) and, through the two channels with
s = 3.3, the combined KIEs of Figure 3d.
"""

import importlib.util
import os

import pytest

from kinisot import compute_kie, parse_gaussian
from kinisot.jobs import load_job, run_job

CASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks", "dykat_allyl_arylation")
T = 313.15
# SI Table 20 electronic energies (hartree)
SI = {
    "allyl_chloride_3": -694.160144,
    "r3_reactant_complex": -3910.736886,
    "s3_reactant_complex": -3910.739063,
    "r3_anti_oa_ts": -3910.712464,
    "s3_anti_oa_ts": -3910.702421,
}
# atom numbers of the allyl carbons, the same in all four Rh structures of a channel
ATOMS = {
    "R": {"C-Cl": 118, "central": 117, "alkene CH": 122, "CH2 next to CH": 121, "CH2 next to C-Cl": 119,
          "reference": 120},
    "S": {"C-Cl": 122, "central": 117, "alkene CH": 118, "CH2 next to CH": 119, "CH2 next to C-Cl": 121,
          "reference": 120},
}  # fmt: skip
# SI Tables 24 (R) and 25 (S): KIE / KIE(ref), semiclassical and with tunnelling
TABLES = {
    "R": {"C-Cl": (1.028, 1.029), "central": (0.998, 0.999), "alkene CH": (1.000, 1.000),
          "CH2 next to CH": (1.002, 1.002), "CH2 next to C-Cl": (1.006, 1.006)},
    "S": {"C-Cl": (1.033, 1.035), "central": (0.994, 0.994), "alkene CH": (0.998, 0.998),
          "CH2 next to CH": (1.003, 1.003), "CH2 next to C-Cl": (1.006, 1.006)},
}  # fmt: skip
# Figure 3d, the paper's combined KIEs with s = 3.3 (case.json)
FIGURE_3D = {"C1": 1.006, "C2": 0.995, "C3": 1.027, "C4": 1.005, "C6": 1.003}
# the tables print three decimals, and the SI gives the reactant complexes' coordinates to three decimals
TOLERANCE = 7e-4


@pytest.fixture(scope="module")
def files():
    return {name: parse_gaussian(os.path.join(CASE, name + ".log")) for name in SI}


def test_jobs_reproduce_si_table_20(files):
    for name, energy in SI.items():
        # the transition structures and free 3 to the printed digits; the reactant complexes' coordinates are rounded
        tolerance = 3e-5 if "complex" in name else 1e-6
        assert files[name].energy == pytest.approx(energy, abs=tolerance), name
    # one imaginary mode each (Gaussian prints 176.3i and 183.3i)
    for channel, imaginary in (("r3", 176.3), ("s3", 183.3)):
        rct, ts = files[channel + "_reactant_complex"], files[channel + "_anti_oa_ts"]
        result = compute_kie(rct=rct, ts=ts, iso="117", temperature=T, scale=1.0)
        assert result.other.light.imaginary == pytest.approx(imaginary, abs=0.1)


@pytest.mark.parametrize("channel", ["R", "S"])
def test_channels_reproduce_si_tables_24_and_25(files, channel):
    rct, ts = files["%s3_reactant_complex" % channel.lower()], files["%s3_anti_oa_ts" % channel.lower()]
    for position, (semiclassical, tunnelling) in TABLES[channel].items():
        result = compute_kie(rct=rct, ts=ts, iso=str(ATOMS[channel][position]),
                             reference=str(ATOMS[channel]["reference"]), temperature=T, scale=1.0)  # fmt: skip
        assert result.kie_relative == pytest.approx(semiclassical, abs=TOLERANCE), (channel, position)
        assert result.kie_tunnel_relative == pytest.approx(tunnelling, abs=TOLERANCE), (channel, position)


def test_job_file_reproduces_figure_3d():
    results = dict(run_job(load_job(os.path.join(CASE, "channels.json")), T))
    assert sorted(results) == sorted(FIGURE_3D)
    for position, value in FIGURE_3D.items():
        result = results[position]
        assert result.selectivity == pytest.approx(3.3)
        assert result.kie_tunnel_relative == pytest.approx(value, abs=TOLERANCE), position
    # C3 is mostly the C-Cl carbon of (S)-3: 1 / KIE = y_R / KIE_R + y_S / KIE_S
    c3 = results["C3"]
    expected = 1 / sum(y / c.kie_tunnel_relative for y, c in zip(c3.shares, c3.channels))
    assert c3.kie_tunnel_relative == pytest.approx(expected, rel=1e-4)


def test_benchmark_runner_computes_the_case_from_its_job_file():
    spec = importlib.util.spec_from_file_location("bench_run", os.path.join(CASE, "..", "run.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    (case,) = module.load_cases(["dykat_allyl_arylation"])
    rows = module.run_case(case)
    assert [r["position"].split()[0] for r in rows] == ["C1", "C2", "C3", "C4", "C6"]
    assert [round(r["computed"], 4) for r in rows] == [1.0055, 0.9955, 1.0264, 1.0053, 1.0035]
    assert all(r["deviation"] is not None for r in rows)
    assert "from the job file `channels.json`" in module.format_case(case, rows)
    # the settings case.json states must agree with the job file's
    case["scale"] = 0.975
    with pytest.raises(ValueError, match="scale"):
        module.run_case(case)
