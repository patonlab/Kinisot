#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ORCA backend, program detection, GoodVibes-based scaling and the frequency self-check."""

import numpy as np
import pytest
from conftest import datapath, write_minimum

from kinisot import HessianInput, KinisotInputError, KinisotParseError, KinisotWarning, cli, compute_kie, load_hessian
from kinisot.backends import detect_program
from kinisot.backends.gaussian import is_linear, parse_gaussian
from kinisot.backends.orca import orca_paths, parse_orca
from kinisot.hessian import linear_from_geometry, mass_weight
from kinisot.scaling import find_scaling_factor

GS_G, TS_G = datapath("gaussian/claisen_gs.out"), datapath("gaussian/claisen_ts.out")
GS_O, TS_O = datapath("orca/claisen_gs.out"), datapath("orca/claisen_ts.hess")


def test_detect_program(tmp_path):
    assert detect_program(GS_G) == "Gaussian"
    assert detect_program(GS_O) == "Orca" and detect_program(TS_O) == "Orca"
    other = tmp_path / "other.out"
    other.write_text("some other program\n")
    assert detect_program(str(other)) == "unknown"
    with pytest.raises(KinisotParseError, match="not recognized"):
        load_hessian(str(other))


def test_orca_paths(tmp_path):
    out, hess = orca_paths(GS_O)
    assert out == GS_O and hess.endswith("claisen_gs.hess")
    out, hess = orca_paths(TS_O)
    assert hess == TS_O and out.endswith("claisen_ts.out")
    lonely = tmp_path / "lonely.out"
    lonely.write_text("* O   R   C   A *\n")
    with pytest.raises(KinisotParseError, match="expected .*lonely.hess"):
        parse_orca(str(lonely))


def test_orca_input_matches_gaussian():
    g, o = parse_gaussian(GS_G), parse_orca(GS_O)
    assert o.program == "Orca" and o.level_of_theory == "B3LYP/6-31G(d)"
    assert o.atomic_numbers == g.atomic_numbers and o.symbols[:3] == ("C", "C", "O")
    # ORCA's standard atomic weights are replaced by the pure-isotope masses Gaussian uses
    assert o.masses == pytest.approx(g.masses)
    assert np.abs(o.hessian - g.hessian).max() < 1e-9
    assert np.abs(o.positions - g.positions).max() < 1e-6
    assert o.program_frequencies is not None and len(o.program_frequencies) == 36
    assert not o.linear


def test_orca_path_reproduces_gaussian_golden():
    r = compute_kie(rct=GS_O, ts=TS_O, iso="5", temperature=393.0, scale=0.961)
    assert r.kie == pytest.approx(1.001895135, rel=1e-8) and r.warnings == ()
    mixed = compute_kie(rct=GS_G, ts=TS_O, iso="5", temperature=393.0, scale=0.961)
    assert mixed.kie == pytest.approx(r.kie, rel=1e-8)
    auto = compute_kie(rct=GS_O, ts=TS_O, iso="5", temperature=393.0, scale=None)
    assert auto.scaling.factor == 0.977 and auto.scaling.level == "B3LYP/6-31G(d)"


def test_orca_cli(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert cli.main(["--rct", GS_O, "--ts", TS_O, "--iso", "5", "-t", "393", "-s", "0.961", "-q"]) == 0
    assert "1.001895" in (tmp_path / "Kinisot_output.dat").read_text()


def test_linear_from_geometry():
    co2 = np.array([[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]])
    assert linear_from_geometry(co2, [15.99491, 12.0, 15.99491])
    assert linear_from_geometry(co2 + 3.0, [15.99491, 12.0, 15.99491])  # translation invariant
    water = np.array([[0.0, 0.76, -0.47], [0.0, 0.0, 0.12], [0.0, -0.76, -0.47]])
    assert not linear_from_geometry(water, [1.00783, 15.99491, 1.00783])
    # agrees with the rotational-constant test on every bundled Gaussian output
    for name in ["claisen_gs", "claisen_ts", "DATS", "diene", "dienophile", "tetramethylcyclohexane"]:
        data = parse_gaussian(datapath("gaussian/%s.out" % name))
        assert linear_from_geometry(data.positions, data.masses) == (is_linear(data.source) == "linear") is False


def test_frequency_self_check_warns_on_inconsistent_data():
    g = parse_gaussian(GS_G)
    wrong = HessianInput(
        g.hessian,
        g.masses,
        g.atomic_numbers,
        source="wrong",
        program_frequencies=tuple(f * 1.1 for f in g.program_frequencies),
    )
    with pytest.warns(KinisotWarning, match="differ from the 36 the program printed"):
        r = compute_kie(rct=wrong, ts=TS_G, iso="5", temperature=393.0, scale=0.961)
    assert any("wrong" in w for w in r.warnings)
    # the Gaussian examples pass the check with any cutoff and scaling
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_kie(rct=GS_G, ts=TS_G, iso="5", temperature=393.0, scale=0.961, imag_cutoff=300.0)


def test_scale_types():
    assert find_scaling_factor("RB3LYP/6-31G(d)", "zpe") == (0.977, find_scaling_factor("RB3LYP/6-31G(d)")[1])
    assert find_scaling_factor("RB3LYP/6-31G(d)", "harm")[0] == 0.991
    assert find_scaling_factor("RB3LYP/6-31G(d)", "fund")[0] == 0.952
    with pytest.raises(KinisotInputError, match="unknown scale type"):
        find_scaling_factor("RB3LYP/6-31G(d)", "anharmonic")
    r = compute_kie(rct=GS_G, ts=TS_G, iso="5", temperature=393.0, scale=None, scale_type="harm")
    assert r.scaling.factor == 0.991 and r.scaling.scale_type == "harm" and "harmonic factor" in r.scaling.messages[0]


def test_scale_type_cli(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert cli.main(["--rct", GS_G, "--ts", TS_G, "--iso", "5", "--scale-type", "fund", "-q"]) == 0
    assert "scaling factor 0.952" in (tmp_path / "Kinisot_output.dat").read_text()


def test_goodvibes_hessian_parity():
    from goodvibes.io import parse_hessian

    for name in ["claisen_gs", "claisen_ts", "DATS"]:
        path = datapath("gaussian/%s.out" % name)
        mine, theirs = parse_gaussian(path), parse_hessian(path)
        assert np.abs(mine.hessian - theirs.hessian).max() < 1e-12
        assert mine.masses == pytest.approx(theirs.masses)
        assert np.abs(mass_weight(mine.hessian, mine.masses) - mass_weight(theirs.hessian, theirs.masses)).max() < 1e-12


def test_pyquiver_parity_needs_verbose_output(tmp_path):
    pyquiver = pytest.importorskip("pyquiver")
    from pyquiver.config import Config
    from pyquiver.kie import KIE_Calculation

    config = Config.from_dict(
        isotopologues={"C5": [(5, 5, "13C")]}, temperature=393.0, scaling=0.961, imag_threshold=50
    )
    try:
        calc = KIE_Calculation(config, GS_G, TS_G)
    except ValueError as err:
        pytest.skip("PyQuiver needs #p Gaussian output; the bundled files are not verbose: %s" % str(err)[:80])
    ours = compute_kie(rct=GS_G, ts=TS_G, iso="5", temperature=393.0, scale=0.961)
    assert calc.to_dict()["C5"]["infinite_parabola"] == pytest.approx(ours.kie_tunnel, rel=2e-5)
    assert pyquiver is not None


def test_unsubstitutable_orca_element_keeps_program_mass(tmp_path):
    # elements outside the substitution table keep ORCA's masses (until the Phase 7 isotope table)
    hess = datapath("orca/claisen_gs.hess")
    import re

    text = re.sub(r"^ O  +15\.99900", " N      14.00700", open(hess).read(), count=1, flags=re.M)
    (tmp_path / "n.hess").write_text(text)
    data = parse_orca(str(tmp_path / "n.hess"))
    assert data.symbols[2] == "N" and data.masses[2] == 14.007 and data.level_of_theory is None


def test_synthetic_gaussian_without_frequencies_skips_check(tmp_path):
    data = parse_gaussian(write_minimum(tmp_path / "m.out"))
    assert data.program_frequencies is None
