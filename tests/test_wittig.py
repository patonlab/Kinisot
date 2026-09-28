#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The Wittig case (benchmarks/wittig_anisaldehyde): Chen, Nieves-Quinones, Waas, Singleton, JACS 2014.

Gaussian 16 jobs at the SI geometries must reproduce the SI energies, and
Kinisot must reproduce the SI Table S4 KIEs for 4-TS (M06-2X/6-31+G**/PCM,
67 C, scaling factor 0.9614, Bell tunnelling).
"""

import os
import warnings

import pytest

from kinisot import KinisotWarning, compute_kie, parse_gaussian

CASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks", "wittig_anisaldehyde")
T, SCALE = 340.15, 0.9614
# SI electronic energies (hartree)
SI = {"anisaldehyde_1": -459.931869643, "ylide_2": -1227.87396966, "ts_4": -1687.79361546}
# position: labels (anisaldehyde, ylide, TS) and SI Table S4 (M06-2X/6-31+G**, 4-TS)
TABLE_S4 = {
    "ylide CH carbon": (["0", "1", "4"], 1.022),
    "carbonyl CHO": (["9", "0", "2"], 1.043),
    "ipso": (["4", "0", "5"], 0.999),
    "ortho": (["5", "0", "24"], 1.001),
    "ortho'": (["3", "0", "28"], 1.001),
    "meta": (["6", "0", "25"], 1.000),
    "meta'": (["2", "0", "27"], 1.000),
    "para": (["1", "0", "26"], 1.000),
    "ketone carbonyl": (["0", "37", "51"], 0.998),
    "ketone methyl": (["0", "38", "52"], 1.000),
    "methoxy": (["13", "0", "57"], 1.000),
    "aldehyde 18O": (["17:18O", "0", "1:18O"], 1.016),
}


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


def test_4ts_reproduces_si_table_s4(files):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)  # the ylide's 9i cm-1 phenyl torsion
        for position, (iso, si) in TABLE_S4.items():
            result = compute_kie(rct=[files["anisaldehyde_1"], files["ylide_2"]], ts=files["ts_4"], iso=iso,
                                 temperature=T, scale=SCALE, project=True)  # fmt: skip
            # Table S4 prints three decimals; several values lie on a rounding boundary
            assert result.kie_tunnel == pytest.approx(si, abs=6e-4), position
