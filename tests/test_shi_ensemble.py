#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""The 18 Shi epoxidation transition structures (benchmarks/shi_epoxidation/ensemble).

Gaussian frequency jobs at the SI geometries of Singleton & Wang (JACS 2005)
must reproduce the SI, Kinisot must reproduce the authors' per-structure
QUIVER predictions (SI Table 1), and the ensemble prototype must stay close
to TS 10, which carries most of the rate.
"""

import importlib.util
import json
import os

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ENSEMBLE = os.path.join(ROOT, "benchmarks", "shi_epoxidation", "ensemble")


@pytest.fixture(scope="module")
def shi():
    pytest.importorskip("goodvibes")
    spec = importlib.util.spec_from_file_location("shi_ensemble", os.path.join(ENSEMBLE, "analyze.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows, ensembles = module.analyze()
    return module, rows, ensembles


def test_frequency_jobs_reproduce_the_si(shi):
    module, rows, _ = shi
    with open(os.path.join(ENSEMBLE, "si_atom_order.json")) as handle:
        si = json.load(handle)
    si["A"] = {"si_energy": -1343.49880272, "si_zpe": 0.466153}
    assert set(rows) == set(module.LABELS) and len(rows) == 18
    for label, row in rows.items():
        # TS AB includes the repair of the SI's misprinted coordinates (ensemble/README.md)
        assert row["energy"] == pytest.approx(si[label]["si_energy"], abs=2e-7), label
        assert row["zpe"] == pytest.approx(si[label]["si_zpe"], abs=1.5e-5), label
        assert row["n_imaginary"] == 1, label


def test_per_structure_kies_reproduce_si_table_1(shi):
    module, rows, _ = shi
    exact = 0
    for label, row in rows.items():
        for site, predicted in zip(module.SITES, module.SI_PREDICTED[label]):
            assert row["absolute"][site] == pytest.approx(predicted, abs=1e-3), (label, site)
            exact += round(row["absolute"][site], 3) == predicted
    assert exact >= 102  # of 108; the other six are within 0.0008


def test_ensemble_is_dominated_by_ts_10(shi):
    module, rows, ensembles = shi
    for name, ensemble in ensembles.items():
        assert 0.84 < ensemble["share"]["A"] < 0.91 and 0.05 < ensemble["share"]["B"] < 0.14, name
        for site in module.SITES:
            assert ensemble["relative"][site] == pytest.approx(rows["A"]["relative"][site], abs=5e-4), (name, site)
