#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The Python API: compute_kie(), the result object, JSON/CSV output, tunnelling models."""

import csv
import json
import warnings

import numpy as np
import pytest
from conftest import datapath, write_minimum, write_ts

from kinisot import (
    HessianInput,
    IsotopeEffect,
    KinisotInputError,
    KinisotWarning,
    ScalingChoice,
    choose_scaling_factor,
    cli,
    compute_kie,
    load_hessian,
    parse_gaussian,
)
from kinisot.thermo import bell_correction, crossover_temperature, tunneling_correction, wigner_correction

GS = datapath("gaussian/claisen_gs.out")
TS = datapath("gaussian/claisen_ts.out")
TMCH = datapath("gaussian/tetramethylcyclohexane.out")


@pytest.fixture(scope="module")
def claisen():
    return compute_kie(rct=GS, ts=TS, iso="5", temperature=393.0, scale=0.961)


def test_result_object(claisen):
    r = claisen
    assert isinstance(r, IsotopeEffect)
    assert r.kind == "KIE" and r.tunneling == "bell"
    assert r.kie_tunnel == pytest.approx(1.001939694, rel=1e-6)
    assert r.kie == pytest.approx(r.imag_ratio * r.zpe * r.exc * r.trpf)
    assert r.kie_tunnel == pytest.approx(r.kie * r.tunnel_corr)
    # the final factors are ratios of the per-side factors
    assert r.zpe == pytest.approx(r.reactant.zpe_factor / r.other.zpe_factor)
    assert r.trpf == pytest.approx(r.reactant.trpf_factor / r.other.trpf_factor)
    assert r.transition_structure is r.other and r.product is None
    assert r.reactant.rpfr == pytest.approx(r.reactant.zpe_factor * r.reactant.exc_factor * r.reactant.trpf_factor)
    assert r.scaling == ScalingChoice(0.961, "user")
    assert r.scale_factor == 0.961
    assert len(r.species) == 4 and r.species[0] is r.reactant.light


def test_species_details(claisen):
    light, heavy = claisen.other.light.species[0], claisen.other.heavy.species[0]
    assert light.name == "claisen_ts" and light.label == "0" and heavy.label == "5"
    assert light.imaginary == pytest.approx(463.9, abs=0.05)
    assert len(light.frequencies) == 35 and len(light.discarded) == 6
    assert [s.atom for s in heavy.substitutions] == [5] and heavy.masses[4] == pytest.approx(13.00335)
    assert claisen.reactant.light.species[0].imaginary is None
    assert claisen.reactant.heavy.labels == ("5",)
    assert claisen.other.name == "claisen_ts"


def test_to_dict_and_json_round_trip(claisen):
    d = claisen.to_dict()
    text = claisen.to_json()
    assert json.loads(text) == d
    assert d["kie_tunnel"] == claisen.kie_tunnel and d["kind"] == "KIE"
    assert d["transition_structure"]["heavy"]["species"][0]["substitutions"] == [
        {"atom": 5, "element": "C", "light_mass": 12.0, "heavy_mass": 13.00335}
    ]
    assert d["reactant"]["light"]["species"][0]["imaginary_frequency"] is None
    assert d["scale_source"] == "user" and d["warnings"] == []
    row = claisen.summary_row()
    assert row["kie_tunnel"] == claisen.kie_tunnel and row["labels"] == "5;5"


def test_eqe_result():
    r = compute_kie(rct=TMCH, prd=TMCH, iso=["24,25,26", "28,29,30"], temperature=290.0, scale=1.0)
    assert r.kind == "EQE" and r.product is r.other and r.transition_structure is None
    assert r.imag_ratio == 1.0 and r.tunnel_corr == 1.0 and r.tunneling == "none"
    assert r.kie == pytest.approx(1.038453539, rel=1e-6)
    assert "product" in r.to_dict() and "transition_structure" not in r.to_dict()


def test_inputs_may_be_hessian_objects(claisen):
    gs, ts = parse_gaussian(GS), parse_gaussian(TS)
    assert isinstance(gs, HessianInput) and gs.program == "Gaussian" and gs.symbols[2] == "O"
    assert gs.positions is not None and gs.positions.shape == (14, 3)
    assert load_hessian(gs) is gs
    r = compute_kie(rct=[gs], ts=[ts], iso=["5", "5"], temperature=393.0, scale=0.961)
    assert r.kie_tunnel == claisen.kie_tunnel
    # a hand-made HessianInput works too (masses, atomic numbers, hessian)
    made = HessianInput(gs.hessian, gs.masses, gs.atomic_numbers, source="made", level_of_theory=gs.level_of_theory)
    assert compute_kie(rct=made, ts=ts, iso="5", temperature=393.0, scale=0.961).kie == claisen.kie
    with pytest.raises(ValueError, match="shape"):
        HessianInput(np.zeros((3, 3)), (1.0, 12.0), (1, 6))
    with pytest.raises(KinisotInputError):
        load_hessian(42)


def test_auto_scaling_records_provenance():
    r = compute_kie(rct=GS, ts=TS, iso="5", temperature=393.0, scale=None)
    assert r.scaling.source == "truhlar" and r.scaling.factor == 0.977
    assert r.scaling.level == "RB3LYP/6-31G(d)" and "Truhlar" in r.scaling.reference
    assert any("0.977" in m for m in r.scaling.messages)
    assert r.to_dict()["scale_source"] == "truhlar"


def test_choose_scaling_factor_mismatch(tmp_path):
    a = parse_gaussian(write_minimum(tmp_path / "a.out", level="RM062X"))
    b = parse_gaussian(write_ts(tmp_path / "b.out"))
    choice = choose_scaling_factor([a, b])
    assert choice.factor == 1.0 and choice.source == "default"
    assert any("not computed at the same level" in m for m in choice.messages)
    assert choose_scaling_factor([a], 0.9) == ScalingChoice(0.9, "user")


def test_label_forms(claisen):
    assert compute_kie(rct=GS, ts=TS, iso=5, temperature=393.0, scale=0.961).kie == claisen.kie
    assert compute_kie(rct=GS, ts=TS, iso=["5", "5"], temperature=393.0, scale=0.961).kie == claisen.kie
    with pytest.raises(KinisotInputError, match="iso is required"):
        compute_kie(rct=GS, ts=TS)
    with pytest.raises(KinisotInputError, match="at least one"):
        compute_kie(rct=[], ts=TS, iso="5")


def test_tunneling_models(claisen):
    im_light, im_heavy = claisen.other.light.imaginary, claisen.other.heavy.imaginary
    assert claisen.tunnel_corr == pytest.approx(bell_correction(im_light, im_heavy, 393.0))
    none = compute_kie(rct=GS, ts=TS, iso="5", temperature=393.0, scale=0.961, tunneling="none")
    assert none.tunnel_corr == 1.0 and none.kie_tunnel == none.kie == claisen.kie
    wigner = compute_kie(rct=GS, ts=TS, iso="5", temperature=393.0, scale=0.961, tunneling="wigner")
    assert wigner.tunnel_corr == pytest.approx(wigner_correction(im_light, im_heavy, 393.0))
    # Bell reduces to Wigner to first order in u^2 (the quartic term is ~1e-5 here)
    assert wigner.tunnel_corr == pytest.approx(claisen.tunnel_corr, rel=2e-5)
    assert tunneling_correction("none", 1.0, 2.0, 300.0) == 1.0
    with pytest.raises(KinisotInputError, match="unknown tunnelling model"):
        compute_kie(rct=GS, ts=TS, iso="5", tunneling="eckart")


def test_bell_below_crossover_temperature():
    t_c = crossover_temperature(1000.0)
    assert t_c == pytest.approx(228.9, abs=0.1)
    assert bell_correction(1000.0, 990.0, t_c + 1.0) > 1.0
    with pytest.raises(KinisotInputError, match="crossover temperature"):
        bell_correction(1000.0, 990.0, t_c - 1.0)
    assert wigner_correction(1000.0, 990.0, t_c - 1.0) > 1.0  # Wigner still defined


def test_warnings_are_emitted_and_recorded(tmp_path):
    from conftest import MASSES_CH2, Z_CH2, synthetic_hessian, write_gaussian_like

    rct = write_minimum(tmp_path / "rct.out")
    ts = write_gaussian_like(
        tmp_path / "ts.out", Z_CH2, MASSES_CH2,
        synthetic_hessian(MASSES_CH2, [-500.0, -80.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3000.0], seed=2), nimag=2,
    )  # fmt: skip
    with pytest.warns(KinisotWarning) as record:
        r = compute_kie(rct=rct, ts=ts, iso="1")
    assert len(record) == 2 and len(r.warnings) == 2
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_kie(rct=GS, ts=TS, iso="5")  # clean inputs raise nothing


def test_cli_json_and_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = [
        "--rct",
        GS,
        "--ts",
        TS,
        "--iso",
        "5",
        "-t",
        "393",
        "-s",
        "0.961",
        "-q",
        "--json",
        "run.json",
        "--csv",
        "runs.csv",
    ]
    assert cli.main(args) == 0
    assert cli.main(args[:-4] + ["--csv", "runs.csv", "--iso", "4", "--tunneling", "wigner"]) == 0
    with open("run.json") as handle:
        data = json.load(handle)
    assert data["kie_tunnel"] == pytest.approx(1.001939694, rel=1e-6)
    assert data["transition_structure"]["light"]["imaginary_frequency"] == pytest.approx(463.9, abs=0.05)
    with open("runs.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2 and rows[0]["labels"] == "5;5" and rows[1]["tunneling"] == "wigner"
    assert float(rows[0]["kie_tunnel"]) == pytest.approx(1.001939694, rel=1e-6)
    with open("Kinisot_output.dat") as handle:
        assert "tunnelling: wigner" in handle.read()
