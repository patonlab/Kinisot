#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Conformer ensembles (kinisot/ensemble.py; IMPLEMENTATION_PLAN.md, Phase 10).

The Shi epoxidation transition structures (benchmarks/shi_epoxidation) are
the test set: 18 conformers of one transition structure, one reactant.
"""

import dataclasses
import json
import os
import warnings

import numpy as np
import pytest

from kinisot import (
    ConformerResult,
    Conformers,
    EnsembleIsotopeEffect,
    IsotopeEffect,
    KinisotInputError,
    KinisotWarning,
    compute_kie,
    equivalent_positions,
    load_hessian,
)
from kinisot.ensemble import GAS_CONSTANT_KCAL
from kinisot.thermo import HARTREE_TO_KCAL_PER_MOL, skodje_truhlar_kappa, tunneling_kappa

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SHI = os.path.join(ROOT, "benchmarks", "shi_epoxidation")
T, SCALE = 273.15, 0.9614
BETA = ["2", "8"]  # C-beta: reactant atom 2, transition-structure atom 8


def ts_path(label):
    return os.path.join(SHI, "ts10.log") if label == "A" else os.path.join(SHI, "ensemble", "ts_%s.log" % label)


@pytest.fixture(scope="module")
def shi():
    reactant = load_hessian(os.path.join(SHI, "methylstyrene.log"))
    return reactant, {label: load_hessian(ts_path(label)) for label in ("A", "B", "D", "G")}


def kie(rct, ts, iso=BETA, **kwargs):
    kwargs.setdefault("temperature", T)
    kwargs.setdefault("scale", SCALE)
    kwargs.setdefault("project", True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)
        return compute_kie(rct=rct, ts=ts, iso=iso, **kwargs)


def test_one_conformer_is_the_ordinary_calculation(shi):
    reactant, ts = shi
    plain = kie(reactant, ts["A"])
    for nested in ([[ts["A"]]], Conformers([ts["A"]]), [Conformers([ts["A"]])]):
        result = kie(reactant, nested)
        assert isinstance(result, IsotopeEffect)
        assert result.kie_tunnel == plain.kie_tunnel and result.kie == plain.kie


def test_the_same_file_twice_changes_nothing_but_warns(shi):
    reactant, ts = shi
    plain = kie(reactant, ts["A"])
    with pytest.warns(KinisotWarning, match="same structure"):
        twice = compute_kie(rct=reactant, ts=[[ts["A"], ts["A"]]], iso=BETA, temperature=T, scale=SCALE, project=True)
    assert isinstance(twice, EnsembleIsotopeEffect)
    assert twice.kie_tunnel == pytest.approx(plain.kie_tunnel, abs=1e-12)
    assert twice.kie == pytest.approx(plain.kie, abs=1e-12)
    assert twice.n_effective == pytest.approx(2.0)
    both = kie([[reactant, reactant]], [[ts["A"], ts["A"]]])
    assert both.kie_tunnel == pytest.approx(plain.kie_tunnel, abs=1e-12)


def test_degeneracy_counts_like_a_duplicate(shi):
    reactant, ts = shi
    weighted = kie(reactant, Conformers([ts["A"], ts["B"]], degeneracy=[2, 1]))
    duplicated = kie(reactant, [[ts["A"], ts["A"], ts["B"]]])
    assert weighted.kie_tunnel == pytest.approx(duplicated.kie_tunnel, abs=1e-13)


def test_the_ensemble_formula(shi):
    """KIE = rho_R / sum_j y_j rho'_j with y_j proportional to kappa_L,j exp(-G_j / RT)."""
    reactant, ts = shi
    energies = [0.0, 0.3, 1.1]  # kcal/mol
    labels = ["A", "B", "G"]
    pairs = [kie(reactant, ts[k]) for k in labels]
    for tunneling in ("bell", "none"):
        result = kie(reactant, Conformers([ts[k] for k in labels], free_energies=energies), weights="user",
                     tunneling=tunneling)  # fmt: skip
        kappa = [tunneling_kappa(tunneling, p.other.light.imaginary, T) for p in pairs]
        y = np.array(kappa) * np.exp(-np.array(energies) / (GAS_CONSTANT_KCAL * T))
        y /= y.sum()
        rho_ts = [p.other.rpfr / (p.imag_ratio * (p.tunnel_corr if tunneling == "bell" else 1.0)) for p in pairs]
        assert result.kie_tunnel == pytest.approx(pairs[0].reactant.rpfr / np.dot(y, rho_ts), abs=1e-13)
        rows = result.rows("transition structure")
        assert [r.population for r in rows] == pytest.approx(list(y), abs=1e-13)
        # each conformer's own KIE against the single reactant is the pairwise one
        expected = [p.kie_tunnel if tunneling == "bell" else p.kie for p in pairs]
        assert [r.kie_tunnel for r in rows] == pytest.approx(expected, abs=1e-12)
        assert [r.free_energy for r in rows] == pytest.approx(energies)
    # kJ/mol and hartree give the same weights
    kcal = kie(reactant, Conformers([ts[k] for k in labels], free_energies=energies), weights="user")
    for unit, factor in (("kJ/mol", 4.184), ("hartree", 1 / 627.5094740631)):
        other = kie(reactant, Conformers([ts[k] for k in labels], [e * factor for e in energies], energy_unit=unit),
                    weights="user")  # fmt: skip
        assert other.kie_tunnel == pytest.approx(kcal.kie_tunnel, abs=1e-12)


def test_weighting_schemes(shi):
    reactant, ts = shi
    labels = ["A", "B", "D", "G"]
    files = [ts[k] for k in labels]
    pairs = {k: kie(reactant, ts[k]) for k in labels}
    qrrho = kie(reactant, [files])
    assert qrrho.weights == "qrrho"
    energies = [r.free_energy for r in qrrho.rows("transition structure")]
    assert min(energies) == 0.0 and energies.index(0.0) == 0  # TS A (10) is the lowest
    assert 1.0 < qrrho.n_effective < 2.0
    assert qrrho.kie_tunnel_range[0] < qrrho.kie_tunnel < qrrho.kie_tunnel_range[1]
    assert qrrho.kie_tunnel_lowest == pytest.approx(pairs["A"].kie_tunnel, abs=1e-12)
    lowest = kie(reactant, [files], weights="lowest")
    assert lowest.kie_tunnel == pytest.approx(pairs["A"].kie_tunnel, abs=1e-12)
    equal = kie(reactant, [files], weights="equal", tunneling="none")
    rho = [p.other.rpfr / p.imag_ratio for p in pairs.values()]
    assert equal.kie == pytest.approx(pairs["A"].reactant.rpfr / np.mean(rho), abs=1e-13)
    assert all(r.free_energy is None for r in equal.rows("transition structure"))
    rrho = kie(reactant, [files], weights="rrho")
    assert rrho.kie_tunnel != qrrho.kie_tunnel
    assert rrho.kie_tunnel == pytest.approx(qrrho.kie_tunnel, abs=5e-4)
    # given free energies are used whenever present
    given = kie(reactant, Conformers(files, free_energies=[0, 5, 5, 5]))
    assert given.kie_tunnel == pytest.approx(pairs["A"].kie_tunnel, abs=2e-5)
    assert {r.weight_source for r in given.rows("transition structure")} == {"user"}


def test_qrrho_free_energies_match_goodvibes(shi):
    goodvibes = pytest.importorskip("goodvibes.api")
    reactant, ts = shi
    labels = ["A", "B", "D", "G"]
    result = kie(reactant, [[ts[k] for k in labels]])
    g = [
        goodvibes.compute_thermo(ts_path(k), temperature=T, freq_scale_factor=SCALE,
                                 zpe_scale_factor=SCALE).qh_gibbs_free_energy * 627.5094740631
        for k in labels
    ]  # fmt: skip
    expected = np.array(g) - min(g)
    computed = [r.free_energy for r in result.rows("transition structure")]
    # GoodVibes reads Gaussian's projected frequencies; Kinisot projects its own
    assert computed == pytest.approx(list(expected), abs=0.006)


def test_reference_and_json(shi):
    reactant, ts = shi
    files = [[ts["A"], ts["B"]]]
    result = kie(reactant, files, reference=["11", "4"])
    meta = kie(reactant, files, iso=["11", "4"])
    assert result.kie_tunnel_relative == pytest.approx(result.kie_tunnel / meta.kie_tunnel, abs=1e-14)
    data = json.loads(result.to_json())
    assert data["ensemble"] is True and data["weights"] == "qrrho"
    assert len(data["conformers"]) == 3 and data["reference"]["kie_tunnel"] == pytest.approx(meta.kie_tunnel)
    assert isinstance(result.conformers[0], ConformerResult) and result.conformers[0].role == "reactant"
    row = result.summary_row()
    assert row["conformers"] == 3 and row["labels"] == "2;8" and row["reference_labels"] == "11;4"


def test_equivalent_positions_average_rho_on_each_side(shi):
    reactant, ts = shi
    # the two meta carbons: one NMR signal, two positions in the transition structure
    pairs = [kie(reactant, ts["A"], iso=iso) for iso in (["11", "4"], ["13", "6"])]
    exact = equivalent_positions(pairs)
    expected = np.mean([p.reactant.rpfr for p in pairs]) / np.mean([p.other.rpfr / p.imag_ratio for p in pairs])
    assert exact[0] == pytest.approx(expected, abs=1e-14)
    assert exact[1] != exact[0]
    # an ensemble of one placement each averages in the same way
    ensembles = [kie(reactant, [[ts["A"], ts["B"]]], iso=iso) for iso in (["11", "4"], ["13", "6"])]
    both = equivalent_positions(ensembles)
    assert min(e.kie_tunnel for e in ensembles) <= both[1] <= max(e.kie_tunnel for e in ensembles)
    with pytest.raises(KinisotInputError):
        equivalent_positions([])


def test_skodje_barriers_are_per_conformer(shi):
    reactant, ts = shi
    # the dioxirane is left out of the reactant files, so place the reactant 1 kcal/mol below TS A: a barrier
    # low enough for the Skodje-Truhlar factor to depend on it
    reactant = dataclasses.replace(reactant, energy=ts["A"].energy - 1.0 / HARTREE_TO_KCAL_PER_MOL)
    plain = kie(reactant, ts["A"], tunneling="skodje")
    assert plain.barrier == pytest.approx(1.0)
    twice = kie(reactant, [[ts["A"], ts["A"]]], tunneling="skodje")
    assert twice.kie_tunnel == pytest.approx(plain.kie_tunnel, abs=1e-12)
    pair = kie(reactant, [[ts["A"], ts["G"]]], tunneling="skodje")
    g = pair.rows("transition structure")[1]
    barrier = 1.0 + (ts["G"].energy - ts["A"].energy) * HARTREE_TO_KCAL_PER_MOL
    assert g.kappa_light == pytest.approx(skodje_truhlar_kappa(g.imaginary_light, T, barrier), rel=1e-12)
    fixed = kie(reactant, [[ts["A"], ts["G"]]], tunneling="skodje", barrier=1.0)
    assert fixed.barrier == 1.0 and fixed.rows("transition structure")[1].kappa_light != g.kappa_light


def test_input_errors(shi):
    reactant, ts = shi
    with pytest.raises(KinisotInputError, match="same atoms"):
        kie(reactant, [[ts["A"], reactant]])
    with pytest.raises(KinisotInputError, match="unknown weights"):
        kie(reactant, [[ts["A"], ts["B"]]], weights="boltzmann")
    with pytest.raises(KinisotInputError, match="needs a free energy"):
        kie(reactant, [[ts["A"], ts["B"]]], weights="user")
    with pytest.raises(KinisotInputError, match="2 conformers but 3 free energies"):
        Conformers([ts["A"], ts["B"]], free_energies=[0, 1, 2])
    with pytest.raises(KinisotInputError, match="unknown energy unit"):
        Conformers([ts["A"]], energy_unit="eV")
    with pytest.raises(KinisotInputError, match="positive"):
        Conformers([ts["A"]], degeneracy=[0])
    with pytest.raises(KinisotInputError, match="no electronic energy"):
        kie(reactant, [[ts["A"], dataclasses.replace(ts["B"], energy=None)]])
    with pytest.raises(KinisotInputError, match="imaginary frequency.*reactant conformer"):
        kie([[ts["A"], ts["B"]]], ts["G"], iso=["8", "8"])


def test_renumbered_conformer_warns(shi):
    reactant, ts = shi
    # swap the coordinates of a methyl hydrogen and a ring hydrogen: the bonds no longer match
    positions = np.array(ts["B"].positions)
    hydrogens = [i for i, z in enumerate(ts["B"].atomic_numbers) if z == 1]
    distance = np.linalg.norm(positions[hydrogens][:, None] - positions[hydrogens][None], axis=-1)
    i, j = np.unravel_index(np.argmax(distance), distance.shape)
    positions[[hydrogens[i], hydrogens[j]]] = positions[[hydrogens[j], hydrogens[i]]]
    renumbered = dataclasses.replace(ts["B"], positions=positions)
    with pytest.warns(KinisotWarning, match="differ in bonding"):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*(differ from|mixes).*")
            compute_kie(rct=reactant, ts=[[ts["A"], renumbered]], iso=BETA, temperature=T, scale=SCALE,
                        weights="equal")  # fmt: skip
