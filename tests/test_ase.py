#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ASE backend: VibrationsData JSON input, calculators, finite-difference Hessians, caching, --calc."""

import os

import numpy as np
import pytest
from conftest import datapath

from kinisot import (
    KinisotInputError,
    KinisotParseError,
    build_calculator,
    cli,
    compute_kie,
    hessian_for_geometry,
    hessian_from_calculator,
    load_hessian,
    parse_gaussian,
    save_hessian_json,
)
from kinisot.backends import detect_program
from kinisot.backends.ase import CALCULATORS, is_ase_json, parse_ase_json

ase = pytest.importorskip("ase")

GS_J, TS_J = datapath("ase/claisen_gs.hessian.json"), datapath("ase/claisen_ts.hessian.json")
GS_G, TS_G = datapath("gaussian/claisen_gs.out"), datapath("gaussian/claisen_ts.out")
KW = dict(temperature=393.0, scale=0.961)


def test_json_input_matches_gaussian():
    j, g = parse_ase_json(GS_J), parse_gaussian(GS_G)
    assert detect_program(GS_J) == "ase" and is_ase_json(GS_J) and not is_ase_json(GS_G)
    assert j.program == "Gaussian" and j.level_of_theory == "RB3LYP/6-31G(d)"  # provenance survives the round trip
    assert j.atomic_numbers == g.atomic_numbers and j.masses == pytest.approx(g.masses, abs=1e-5)
    assert np.abs(j.hessian - g.hessian).max() < 1e-10 and np.abs(j.positions - g.positions).max() < 1e-9
    assert j.energy == pytest.approx(g.energy, abs=1e-9)
    assert len(j.program_frequencies) == 36
    assert np.abs(np.sort(j.program_frequencies) - np.sort(g.program_frequencies)).max() < 0.05


def test_json_path_reproduces_golden():
    r = compute_kie(rct=GS_J, ts=TS_J, iso="5", project=False, **KW)
    gaussian = compute_kie(rct=GS_G, ts=TS_G, iso="5", project=False, **KW)
    assert r.kie == pytest.approx(gaussian.kie, rel=1e-10) and r.warnings == ()
    mixed = compute_kie(rct=GS_G, ts=TS_J, iso="5", project=False, **KW)
    assert mixed.kie == pytest.approx(r.kie, rel=1e-10)
    assert compute_kie(rct=GS_J, ts=TS_J, iso="5", scale=None, temperature=393.0).scaling.factor == 0.977


def test_projection_default_follows_the_backend(tmp_path):
    # Gaussian and Gaussian-derived JSON: off; Hessians computed by ASE: on
    assert compute_kie(rct=GS_J, ts=TS_J, iso="5", **KW).project is False
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS

    atoms = molecule("H2O")  # a genuine minimum on the EMT surface (CH4 and C2H6 are not)
    atoms.calc = EMT()
    BFGS(atoms, logfile=None).run(fmax=1e-5)
    data = hessian_from_calculator(atoms, EMT(), delta=0.01, nfree=2)
    assert data.program.startswith("ase") and data.energy is not None and data.positions.shape == (3, 3)
    assert data.masses == pytest.approx((15.99491462, 1.00782503, 1.00782503), abs=1e-7)
    assert len(data.program_frequencies) == 3 and min(data.program_frequencies) > 0
    r = compute_kie(rct=data, prd=data, iso=["2", "3"], temperature=300.0)
    assert r.project is True and r.kind == "EQE"
    assert r.kie == pytest.approx(1.0, abs=1e-6)  # equivalent hydrogens: no isotope effect
    assert max(abs(f) for f in r.reactant.light.species[0].discarded) < 0.5  # projected external modes
    assert compute_kie(rct=data, prd=data, iso=["2", "3"], temperature=300.0, project=False).project is False
    # a linear molecule computed with ASE: five external modes, one vibration
    n2 = molecule("N2")
    n2.calc = EMT()
    BFGS(n2, logfile=None).run(fmax=1e-5)
    linear = hessian_from_calculator(n2, EMT())
    assert linear.linear and len(linear.program_frequencies) == 1
    r = compute_kie(rct=linear, prd=linear, iso=["1", "2"], temperature=300.0)
    assert r.kie == pytest.approx(1.0, abs=1e-9) and len(r.reactant.light.species[0].discarded) == 5


def test_save_and_reload_json(tmp_path):
    g = parse_gaussian(GS_G)
    path = save_hessian_json(g, str(tmp_path / "gs.hessian.json"))
    j = parse_ase_json(path)
    assert np.abs(j.hessian - g.hessian).max() < 1e-10 and j.energy == pytest.approx(g.energy, abs=1e-9)
    assert compute_kie(rct=path, ts=TS_J, iso="5", project=False, **KW).kie == pytest.approx(1.001895135, rel=1e-8)
    bad = tmp_path / "bad.json"
    bad.write_text('{"atoms": 1, "hessian": 2}')
    with pytest.raises(KinisotParseError):
        load_hessian(str(bad))


def test_calc_from_geometry_with_cache(tmp_path, monkeypatch):
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from ase.io import write
    from ase.optimize import BFGS

    atoms = molecule("H2O")
    atoms.calc = EMT()
    BFGS(atoms, logfile=None).run(fmax=1e-5)
    geometry = str(tmp_path / "water.xyz")
    write(geometry, atoms)
    with pytest.raises(KinisotParseError, match="needs --calc"):
        load_hessian(geometry)
    data = hessian_for_geometry(geometry, "emt")
    cached = str(tmp_path / "water.hessian.json")
    assert os.path.exists(cached) and data.source == cached and data.program.startswith("ase")
    stamp = os.path.getmtime(cached)
    again = hessian_for_geometry(geometry, "emt")  # reused, not recomputed
    assert os.path.getmtime(cached) == stamp and np.abs(again.hessian - data.hessian).max() == 0
    other = hessian_for_geometry(geometry, "ase.calculators.emt:EMT")  # a different --calc spec recomputes
    assert np.abs(other.hessian - data.hessian).max() < 1e-8
    # command line: deuterium on one or another of the equivalent hydrogens, the EQE must be 1
    monkeypatch.chdir(tmp_path)
    args = ["--rct", geometry, "--prd", geometry, "--iso", "2", "--iso", "3", "--calc", "emt", "-q", "--json", "e.json"]
    assert cli.main(args) == 0
    text = (tmp_path / "Kinisot_output.dat").read_text()
    assert "external modes projected" in text and "EQE @ 298.15 K" in text
    import json

    assert json.load(open(tmp_path / "e.json"))["kie"] == pytest.approx(1.0, abs=1e-6)


def test_build_calculator_specs():
    from ase.calculators.emt import EMT

    assert isinstance(build_calculator("emt"), EMT)
    assert isinstance(build_calculator("ase.calculators.emt:EMT"), EMT)
    assert isinstance(load_hessian(GS_J, calculator="emt"), type(parse_ase_json(GS_J)))  # JSON wins over --calc
    for name in ("mace_mp", "mace_off", "orb", "sevennet", "aimnet2"):
        assert name in CALCULATORS
    with pytest.raises(KinisotInputError, match="unknown calculator"):
        build_calculator("not-a-calculator")
    with pytest.raises(KinisotInputError, match="install"):
        build_calculator("no_such.module_kinisot:thing")


def test_analytic_hessian_is_used_when_available():
    from ase.build import molecule
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS

    atoms = molecule("H2O")
    atoms.calc = EMT()
    BFGS(atoms, logfile=None).run(fmax=1e-5)
    reference = hessian_from_calculator(atoms, EMT())

    class WithHessian(EMT):
        calls = 0

        def get_hessian(self, atoms):
            WithHessian.calls += 1
            from ase import units

            return reference.hessian * units.Hartree / units.Bohr**2  # eV/A^2, as ASE calculators return

    data = hessian_from_calculator(atoms, WithHessian(), analytic=True)
    assert WithHessian.calls == 1 and np.abs(data.hessian - reference.hessian).max() < 1e-12
    data = hessian_from_calculator(atoms, WithHessian(), analytic=False)
    assert WithHessian.calls == 1


@pytest.mark.skipif(pytest.importorskip("importlib").util.find_spec("mace") is None, reason="mace-torch not installed")
def test_mace_calculator_builds():
    build_calculator("mace_mp:small")
