#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The nitroarene 18O KIEs (benchmarks/nitroarene_phosphetane).

ORCA 6.1.0 frequency jobs at the SI geometries of Kang & Radosevich
(Tetrahedron 2025) must reproduce ORCA's frequencies once the external modes
are projected out, and the authors' PyQuiver predictions (Table 3; SI Tables
S2 and S3). For the monotopic TS1B, the paper's singly labelled value is the
attacked oxygen alone; averaged over the two equivalent oxygens, as the
experiment measures, it falls within the measured 1.033 +/- 0.003.
"""

import importlib.util
import os
import warnings

import pytest

from kinisot import KinisotWarning, compute_kie, load_hessian

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CASE = os.path.join(ROOT, "benchmarks", "nitroarene_phosphetane")
# transition structure: (its nitrogen, oxygen A, oxygen B); nitrobenzene is N 12, O 13 and 14
ATOMS = {"ts2": (14, 12, 13), "ts1b": (1, 2, 3)}


@pytest.fixture(scope="module")
def files():
    return {name: load_hessian(os.path.join(CASE, name + ".out")) for name in ("phno2", "ts2", "ts1b")}


def kie(files, ts, oxygens, tunneling="none", project=True):
    """KIE of the 18O label(s) relative to [15N]-nitrobenzene; ``oxygens`` indexes (A, B) of the TS."""
    n, *ts_oxygens = ATOMS[ts]
    rct = ",".join(["12:15N"] + ["%d:18O" % (13 + i) for i in oxygens])
    other = ",".join(["%d:15N" % n] + ["%d:18O" % ts_oxygens[i] for i in oxygens])
    result = compute_kie(
        rct=[files["phno2"]], ts=[files[ts]], iso=[rct, other], reference=["12:15N", "%d:15N" % n],
        temperature=393.0, scale=0.9614, tunneling=tunneling, project=project,
    )  # fmt: skip
    return result.kie_relative if tunneling == "none" else result.kie_tunnel_relative


def test_projected_frequencies_match_orca(files):
    # the SI geometries are not exactly stationary at ORCA 6.1.0: unprojected, a 32 cm-1 mode of TS2 is
    # 1.3 cm-1 off ORCA's (projected) value and the self-check warns; projected, every frequency agrees
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        for ts in ATOMS:
            kie(files, ts, (0,))
    with pytest.warns(KinisotWarning, match="differ from the 123"):
        kie(files, "ts2", (0,), project=False)
    # projection moves the KIEs by less than 4e-6
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert kie(files, "ts2", (0, 1), project=False) == pytest.approx(kie(files, "ts2", (0, 1)), abs=4e-6)


def test_ts2_reproduces_the_paper(files):
    one = [kie(files, "ts2", (i,)) for i in (0, 1)]
    assert one == pytest.approx([1.03199, 1.03246], abs=1e-5)
    assert kie(files, "ts2", (0, 1)) == pytest.approx(1.06607, abs=1e-5)
    # PyQuiver (SI Table S2): 1.0321 and 1.0657 uncorrected, 1.0343 and 1.0703 with the inverted parabola
    assert one[0] == pytest.approx(1.0321, abs=2e-4)
    assert kie(files, "ts2", (0, 1)) == pytest.approx(1.0657, abs=6e-4)
    assert kie(files, "ts2", (0,), "bell") == pytest.approx(1.0343, abs=2e-4)
    assert kie(files, "ts2", (0, 1), "bell") == pytest.approx(1.0703, abs=6e-4)


def test_ts1b_singly_labelled_kie_is_an_average_over_both_oxygens(files):
    attacked, spectator = kie(files, "ts1b", (0,)), kie(files, "ts1b", (1,))
    assert (attacked, spectator) == pytest.approx((1.04738, 1.01397), abs=1e-5)
    # the paper's 1.0468 (SI Table S3) is the attacked oxygen alone
    assert attacked == pytest.approx(1.0468, abs=1e-3)
    assert abs(spectator - 1.0468) > 0.03
    # a singly labelled molecule reacts through both isotopomers: harmonic mean of the two
    average = 2 / (1 / attacked + 1 / spectator)
    assert average == pytest.approx(1.0304, abs=1e-4)
    assert abs(average - 1.033) < 0.003  # within the measured value's error
    # and the doubly labelled KIE is close to the product of the two positions, as for TS2
    both = kie(files, "ts1b", (0, 1))
    assert both == pytest.approx(1.06214, abs=1e-5)
    assert both == pytest.approx(attacked * spectator, abs=5e-4)
    assert both == pytest.approx(1.0612, abs=1.1e-3)  # PyQuiver


def test_benchmark_cases():
    spec = importlib.util.spec_from_file_location("bench_run", os.path.join(ROOT, "benchmarks", "run.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    ts2, ts1b = module.load_cases(["nitroarene_phosphetane", "nitroarene_phosphetane_ts1b"])
    rows = module.run_case(ts2)
    assert [round(r["semiclassical"], 4) for r in rows] == [1.0322, 1.0661]
    assert [round(r["computed"], 4) for r in rows] == [1.0346, 1.0709]
    rows = module.run_case(ts1b)
    # the runner averages equivalent positions geometrically (1.0305); the exact average is 1.0304
    assert [round(r["semiclassical"], 4) for r in rows] == [1.0305, 1.0621]
    assert ts1b["alternative_to"] == "nitroarene_phosphetane"
    assert "left out of the overall mean" in module.format_case(ts1b, rows)
