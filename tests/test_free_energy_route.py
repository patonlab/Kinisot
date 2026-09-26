#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The numbers quoted in README.md and docs/theory.md (section 6) for Bigeleisen-Mayer versus free-energy KIEs."""

import importlib.util
import os

import numpy as np
import pytest

from kinisot import compute_kie

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "compare_free_energy_route.py")


@pytest.fixture(scope="module")
def route():
    spec = importlib.util.spec_from_file_location("compare_free_energy_route", PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, module.load()


@pytest.mark.parametrize("label", ["1", "3:18O", "7,8"])
def test_bigeleisen_mayer_is_kinisot_and_free_energy_differs_by_the_product_rule(route, label):
    c, (gs, ts) = route
    bm = c.bigeleisen_mayer(gs, ts, label)
    kinisot = compute_kie(rct=gs, ts=ts, iso=label, temperature=393.0, scale=1.0, tunneling="none", project=True)
    assert bm == pytest.approx(kinisot.kie, rel=1e-10)
    # exact identity: ln(KIE_FE / KIE_BM) = -(violation in the reactant - violation in the TS)
    gap = np.log(c.free_energy_kie(gs, ts, label) / bm)
    violation = c.product_rule_violation(gs, label) - c.product_rule_violation(ts, label)
    assert gap == pytest.approx(-violation, abs=1e-9)


def test_violation_on_converged_gaussian_files(route):
    c, (gs, ts) = route
    assert abs(c.product_rule_violation(gs, "7,8")) == pytest.approx(1.0e-3, abs=0.1e-3)
    assert max(abs(c.product_rule_violation(d, label)) for d in (gs, ts) for label in ("1", "4", "6")) < 4e-4


def test_sensitivity_to_one_soft_mode(route):
    # +0.1 cm-1 on the heavy isotopologue's 70.5 cm-1 reactant torsion
    c, (gs, ts) = route
    assert c.modes(gs, "")[1][0] == pytest.approx(70.5, abs=0.1)
    bm, fe = c.bigeleisen_mayer(gs, ts, "1"), c.free_energy_kie(gs, ts, "1")
    assert abs(c.bigeleisen_mayer(gs, ts, "1", shift=(gs, 0.1)) - bm) < 1e-5
    assert fe - c.free_energy_kie(gs, ts, "1", shift=(gs, 0.1)) == pytest.approx(1.45e-3, abs=0.05e-3)


def test_quasi_harmonic_free_energies_are_not_for_isotope_effects(route):
    c, (gs, ts) = route
    for label in ("1", "3:18O", "7,8"):
        bm = c.bigeleisen_mayer(gs, ts, label)
        assert abs(c.bigeleisen_mayer(gs, ts, label, quasi_harmonic=100.0) - bm) < 2e-4
    fe = c.free_energy_kie(gs, ts, "3:18O")
    assert fe - c.free_energy_kie(gs, ts, "3:18O", quasi_harmonic=100.0) == pytest.approx(0.0189, abs=0.0005)


def test_printed_precision(route):
    c, (gs, ts) = route
    kt = c.BOLTZMANN_CONSTANT * 393.0 / c.ENERGY_AU
    assert np.expm1(1e-6 / kt) == pytest.approx(8.0e-4, abs=0.05e-4)
    rounded, exact = c.free_energy_kie(gs, ts, "6", decimals=6), c.free_energy_kie(gs, ts, "6")
    assert abs(rounded - exact) == pytest.approx(6.5e-4, abs=0.1e-4)
