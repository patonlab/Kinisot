#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Cross-validation against PyQuiver (pip install pyquiver-kie).

Both programs are fed the same Cartesian Hessians (Kinisot's Gaussian
reader, handed to PyQuiver as a ``quiver.System``), so any difference comes
from the physics code alone. PyQuiver uses five-decimal isotope masses and
an atomic mass unit of 1.660468e-27 kg (a transcription error Kinisot
removed in 2.0.3); together these account for the remaining ~2e-6.
"""

import logging

import numpy as np
import pytest
from conftest import datapath

from kinisot import compute_kie, parse_gaussian
from kinisot.thermo import BOHR_TO_ANGSTROM

pyquiver = pytest.importorskip("pyquiver")
from pyquiver import quiver  # noqa: E402
from pyquiver.config import Config  # noqa: E402
from pyquiver.kie import KIE_Calculation  # noqa: E402

TOLERANCE = 5e-6


def as_pyquiver_system(data):
    """A pyquiver.quiver.System built from a Kinisot HessianInput (no file round trip)."""
    system = object.__new__(quiver.System)
    system.filename = data.source
    system.is_linear = True
    system.atomic_numbers = list(data.atomic_numbers)
    system.number_of_atoms = data.natoms
    system.hessian = np.array(data.hessian)
    system.positions = np.array(data.positions)
    system.positions_angstrom = np.array(data.positions) * BOHR_TO_ANGSTROM
    system._detect_linear()
    return system


def pyquiver_kies(gs, ts, isotopologues, temperature, scaling):
    logging.getLogger("pyquiver").setLevel(logging.ERROR)
    config = Config.from_dict(isotopologues=isotopologues, temperature=temperature, scaling=scaling, imag_threshold=50)
    return KIE_Calculation(config, as_pyquiver_system(gs), as_pyquiver_system(ts)).to_dict()


# name -> (Kinisot label, PyQuiver rules)
CLAISEN = {
    "C1": ("1", [(1, 1, "13C")]),
    "C2": ("2", [(2, 2, "13C")]),
    "O3": ("3:18O", [(3, 3, "18O")]),
    "O3_17": ("3:17O", [(3, 3, "17O")]),
    "C4": ("4", [(4, 4, "13C")]),
    "C5": ("5", [(5, 5, "13C")]),
    "C6": ("6", [(6, 6, "13C")]),
    "H7H8": ("7,8", [(7, 7, "2D"), (8, 8, "2D")]),
}


@pytest.fixture(scope="module")
def claisen():
    gs = parse_gaussian(datapath("gaussian/claisen_gs.out"))
    ts = parse_gaussian(datapath("gaussian/claisen_ts.out"))
    rules = {name: spec[1] for name, spec in CLAISEN.items()}
    return gs, ts, pyquiver_kies(gs, ts, rules, 393.0, 0.961)


@pytest.mark.parametrize("name", list(CLAISEN))
def test_claisen_matches_pyquiver(claisen, name):
    gs, ts, reference = claisen
    label = CLAISEN[name][0]
    bell = compute_kie(rct=gs, ts=ts, iso=label, temperature=393.0, scale=0.961)
    wigner = compute_kie(rct=gs, ts=ts, iso=label, temperature=393.0, scale=0.961, tunneling="wigner")
    assert bell.kie == pytest.approx(reference[name]["uncorrected"], rel=TOLERANCE)
    assert bell.kie_tunnel == pytest.approx(reference[name]["infinite_parabola"], rel=TOLERANCE)
    assert wigner.kie_tunnel == pytest.approx(reference[name]["wigner"], rel=TOLERANCE)


def test_diels_alder_matches_pyquiver():
    gs = parse_gaussian(datapath("gaussian/DATS_rct.out"))
    ts = parse_gaussian(datapath("gaussian/DATS.out"))
    reference = pyquiver_kies(gs, ts, {"C15": [(15, 15, "13C")], "C19": [(19, 19, "13C")]}, 298.15, 0.963)
    for name, label in (("C15", "15"), ("C19", "19")):
        r = compute_kie(rct=gs, ts=ts, iso=label, temperature=298.15, scale=0.963)
        assert r.kie == pytest.approx(reference[name]["uncorrected"], rel=TOLERANCE)
        assert r.kie_tunnel == pytest.approx(reference[name]["infinite_parabola"], rel=TOLERANCE)
