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
from kinisot.isotopes import light_masses
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
    # ORCA reports standard atomic weights; both isotopologues are built from Kinisot's table anyway
    assert o.masses[0] == 12.011 and g.masses[0] == 12.0
    assert light_masses(o) == light_masses(g)
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


def test_unsubstitutable_orca_element_keeps_program_mass(tmp_path):
    # elements outside the substitution table keep ORCA's masses (until the Phase 7 isotope table)
    hess = datapath("orca/claisen_gs.hess")
    import re

    text = re.sub(r"^ O  +15\.99900", " N      14.00700", open(hess).read(), count=1, flags=re.M)
    (tmp_path / "n.hess").write_text(text)
    data = parse_orca(str(tmp_path / "n.hess"))
    assert data.symbols[2] == "N" and data.masses[2] == 14.007 and data.level_of_theory is None
    assert light_masses(data)[2] == pytest.approx(14.003074, abs=1e-6)


def test_synthetic_gaussian_without_frequencies_skips_check(tmp_path):
    data = parse_gaussian(write_minimum(tmp_path / "m.out"))
    assert data.program_frequencies is None


# Real ORCA 6.1.0 outputs (tests/data/orca/README.md): n-pentane conformers and a hydrogen-atom transfer
PENTANE_TT, PENTANE_GG = datapath("orca/pentane_TT.out"), datapath("orca/pentane_GG.out")
HAT_GS, HAT_TS = datapath("orca/hat_gs_freq.out"), datapath("orca/hat_ts_freq.out")


def test_real_orca_conformer_eqe():
    # TT and GG n-pentane at r2SCAN-3c, which has no scaling factor: EQE for 2H at atom 6
    tt = load_hessian(PENTANE_TT)
    assert tt.program == "Orca" and tt.level_of_theory == "r2SCAN-3c/def2-mTZVPP" and len(tt.atomic_numbers) == 17
    r = compute_kie(rct=PENTANE_TT, prd=PENTANE_GG, iso="6", temperature=298.15, scale=None)
    assert r.scaling.factor == 1.0 and r.warnings == ()
    # the value first pinned for these files, 1.007730117, came from slightly different isotope masses
    assert r.kie == pytest.approx(1.0077302, abs=2e-7)
    via_hess = compute_kie(
        rct=PENTANE_TT.replace(".out", ".hess"), prd=PENTANE_GG.replace(".out", ".hess"), iso="6", temperature=298.15
    )
    assert via_hess.kie == pytest.approx(r.kie, abs=1e-12)


def test_real_orca_deep_tunnelling_hat():
    # Broken-symmetry M06-2X-D3(0)/6-31+G**, SMD(dichloromethane). H26 moves from C13 to C6 in the transition
    # mode (1974.9i cm-1, 1911.7i with the detected factor 0.968), so the crossover temperature is 438 K
    kw = dict(rct=HAT_GS, ts=HAT_TS, temperature=298.15, scale=None)
    primary = compute_kie(iso="26", tunneling="wigner", **kw)
    assert primary.scaling.factor == 0.968 and primary.warnings == ()
    assert primary.other.light.imaginary == pytest.approx(1911.7, abs=0.1)
    assert primary.kie == pytest.approx(3.18884, abs=1e-5) and primary.kie_tunnel == pytest.approx(4.46248, abs=1e-5)
    tritium = compute_kie(iso="26:T", tunneling="none", **kw)
    assert tritium.kie == pytest.approx(5.11494, abs=1e-5) and tritium.kie > primary.kie
    assert compute_kie(iso="24", tunneling="none", **kw).kie == pytest.approx(1.037178, abs=1e-6)  # H24 on C13

    # below the crossover temperature the Bell correction is undefined; the Skodje-Truhlar one is not
    with pytest.raises(KinisotInputError, match="crossover temperature 437.8 K"):
        compute_kie(iso="26", tunneling="bell", **kw)
    assert np.isfinite(compute_kie(iso="26", tunneling="skodje", **kw).kie_tunnel)
    # above it both are defined, and the truncated parabola tunnels less than the infinite one
    hot = dict(kw, temperature=500.0)
    bell = compute_kie(iso="26", tunneling="bell", **hot).kie_tunnel
    skodje = compute_kie(iso="26", tunneling="skodje", **hot).kie_tunnel
    assert bell == pytest.approx(5.67082, abs=1e-5) and skodje == pytest.approx(4.85377, abs=1e-5) and skodje < bell
