#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Transition structures in series and parallel channels (kinisot/pathways.py; IMPLEMENTATION_PLAN.md, Phase 10).

The formulas are checked against the published combinations of Chen et al.
(Wittig, JACS 2014; series) and van Dijk et al. (DyKAT arylation, Nat.
Catal. 2021; channels), against a kinetic scheme solved exactly, and, with
real files, against the conformer-ensemble code they must agree with.
"""

import json
import os
import warnings

import numpy as np
import pytest

from kinisot import (
    ChannelIsotopeEffect,
    Conformers,
    KinisotInputError,
    KinisotWarning,
    Series,
    SeriesIsotopeEffect,
    channel_kie,
    channels,
    compute_kie,
    load_hessian,
    series_kie,
)
from kinisot.ensemble import GAS_CONSTANT_KCAL
from kinisot.thermo import tunneling_kappa

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SHI = os.path.join(ROOT, "benchmarks", "shi_epoxidation")
T, SCALE = 273.15, 0.9614
OPTIONS = dict(temperature=T, scale=SCALE, project=True)


@pytest.fixture(scope="module")
def shi():
    """Styrene and three of its epoxidation transition structures, used as steps and channels."""
    reactant = load_hessian(os.path.join(SHI, "methylstyrene.log"))
    ts = {k: load_hessian(os.path.join(SHI, "ensemble", "ts_%s.log" % k)) for k in ("B", "G")}
    ts["A"] = load_hessian(os.path.join(SHI, "ts10.log"))
    return reactant, ts


def quiet(function, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)
        return function(*args, **kwargs)


# ------------------------------------------------------------------ the formulas against the papers


def test_series_formula_reproduces_the_wittig_predictions():
    """Chen et al. 2014, Table 1: 4-TS and 6-TS alone, weighted by free energies and by trajectories."""
    for alone, by_energy, by_trajectories in (((1.043, 1.015), 1.0280, 1.0326), ((1.022, 0.994), 1.0070, 1.0116)):
        # Figure 2: 25.9 and 26.0 kcal/mol at 340.15 K, so C_f = exp(-0.1 / RT) = 0.862
        assert series_kie(alone, free_energies=[25.9, 26.0], temperature=340.15) == pytest.approx(by_energy, abs=5e-5)
        # 128 trajectories go on to product, 76 go back
        assert series_kie(alone, commitment=128 / 76) == pytest.approx(by_trajectories, abs=5e-5)
    # the printed 1.028 / 1.033 and 1.008 / 1.012 are within the rounding of the inputs
    assert series_kie([1.043, 1.015], free_energies=[25.9, 26.0], temperature=340.15) == pytest.approx(1.028, abs=5e-4)
    assert series_kie([1.022, 0.994], commitment=128 / 76) == pytest.approx(1.012, abs=5e-4)
    # the free energies in kJ/mol give the same
    kj = series_kie([1.043, 1.015], free_energies=[25.9 * 4.184, 26.0 * 4.184], temperature=340.15,
                    energy_unit="kJ/mol")  # fmt: skip
    assert kj == pytest.approx(series_kie([1.043, 1.015], free_energies=[25.9, 26.0], temperature=340.15))


def test_series_limits():
    kies = [1.043, 1.015]
    assert series_kie(kies, commitment=1e9) == pytest.approx(1.043, abs=1e-8)  # the betaine always goes on
    assert series_kie(kies, commitment=1e-9) == pytest.approx(1.015, abs=1e-8)  # it always returns
    assert series_kie(kies, free_energies=[20.0, 30.0], temperature=300.0) == pytest.approx(1.015, abs=1e-6)
    assert series_kie([1.03], free_energies=[5.0]) == pytest.approx(1.03)
    with pytest.raises(KinisotInputError):
        series_kie(kies)
    with pytest.raises(KinisotInputError):
        series_kie(kies, free_energies=[0, 1], commitment=1)
    with pytest.raises(KinisotInputError):
        series_kie([1, 1, 1], commitment=1)


def slowest_rate(energies, rt):
    """R <=> I1 <=> I2 -> P: the slowest eigenvalue of the rate matrix, with transition-state theory rates."""

    def k(ts, start):
        return np.exp(-(energies[ts] - energies[start]) / rt)

    m = np.array([
        [-k("TS1", "R"), k("TS1", "I1"), 0.0],
        [k("TS1", "R"), -k("TS1", "I1") - k("TS2", "I1"), k("TS2", "I2")],
        [0.0, k("TS2", "I1"), -k("TS2", "I2") - k("TS3", "I2")],
    ])  # fmt: skip
    return -max(np.linalg.eigvals(m).real)


def test_series_formula_matches_the_full_kinetic_scheme():
    """Three steps for both isotopologues, solved exactly, against the steady-state formula."""
    rng = np.random.default_rng(7)
    temperature = 298.15
    rt = GAS_CONSTANT_KCAL * temperature
    for _ in range(5):
        # light free energies (kcal/mol); the intermediates are high enough for the steady state to hold
        light = {"R": 0.0, "I1": rng.uniform(9, 11), "I2": rng.uniform(9, 11)}
        light.update({"TS1": rng.uniform(17, 20), "TS2": rng.uniform(17, 20), "TS3": rng.uniform(17, 20)})
        # the heavy isotopologue: every state shifted a little (its free energy of substitution)
        shift = {state: rng.uniform(-0.05, 0.05) for state in light}
        heavy = {state: light[state] + shift[state] for state in light}
        exact = slowest_rate(light, rt) / slowest_rate(heavy, rt)
        # KIE_n from R to TS n, weighted by the light free energies of the transition structures
        steps = ("TS1", "TS2", "TS3")
        kies = [np.exp((shift[ts] - shift["R"]) / rt) for ts in steps]
        combined = series_kie(kies, free_energies=[light[ts] for ts in steps], temperature=temperature)
        assert combined == pytest.approx(exact, abs=1e-6)


def test_channel_formula_reproduces_the_dykat_predictions():
    """van Dijk et al. 2021, Figure 3d: (S)-3 and (R)-3 channels with s = k_S / k_R = 3.3."""
    # per-channel KIEs (Supplementary Tables 24 and 25) of the product carbons: (S)-3, (R)-3
    channel_kies = {
        "C1": (0.998, 1.029), "C2": (0.994, 0.999), "C3": (1.035, 1.000), "C4": (1.006, 1.002), "C6": (1.003, 1.006),
    }  # fmt: skip
    combined = {"C1": 1.0050, "C2": 0.9952, "C3": 1.0266, "C4": 1.0051, "C6": 1.0037}
    figure = {"C1": 1.006, "C2": 0.995, "C3": 1.027, "C4": 1.005, "C6": 1.003}
    s = 3.3
    for carbon, (k_s, k_r) in channel_kies.items():
        value = channel_kie([k_s, k_r], [s, 1.0])
        assert value == pytest.approx(combined[carbon], abs=5e-5), carbon
        assert value == pytest.approx(figure[carbon], abs=1.1e-3), carbon
        # Dale et al. eq S111
        assert value == pytest.approx((1 + s) * k_r * k_s / (k_s + s * k_r), abs=1e-12)
    assert channel_kie([1.02, 1.04], [1, 0]) == pytest.approx(1.02)
    assert channel_kie([1.02, 1.04], [0, 1]) == pytest.approx(1.04)
    with pytest.raises(KinisotInputError):
        channel_kie([1.02, 1.04], [1])
    with pytest.raises(KinisotInputError):
        channel_kie([1.02, 1.04], [0, 0])


# ------------------------------------------------------------------ series from files


def pair(reactant, ts, iso=("2", "8"), **kwargs):
    return quiet(compute_kie, rct=reactant, ts=ts, iso=list(iso), **{**OPTIONS, **kwargs})


def test_series_from_files(shi):
    reactant, ts = shi
    a, g = pair(reactant, ts["A"]), pair(reactant, ts["G"])
    kappa = [tunneling_kappa("bell", p.other.light.imaginary, T) for p in (a, g)]
    rt = GAS_CONSTANT_KCAL * T

    # given free energies: w_n proportional to exp(G_n / RT) / kappa_n
    result = quiet(compute_kie, rct=reactant, ts=Series([ts["A"], ts["G"]], free_energies=[0.0, 1.0]),
                   iso=["2", "8", "8"], **OPTIONS)  # fmt: skip
    assert isinstance(result, SeriesIsotopeEffect)
    w = np.array([1.0 / kappa[0], np.exp(1.0 / rt) / kappa[1]])
    assert result.kie_tunnel == pytest.approx(np.dot(w / w.sum(), [a.kie_tunnel, g.kie_tunnel]), abs=1e-13)
    w = np.array([1.0, np.exp(1.0 / rt)])
    assert result.kie == pytest.approx(np.dot(w / w.sum(), [a.kie, g.kie]), abs=1e-13)
    assert result.steps[0].kie_tunnel == a.kie_tunnel and result.weight_source == "user"
    assert result.kie_tunnel_range[0] < result.kie_tunnel < result.kie_tunnel_range[1]

    # computed free energies agree with the conformer code's
    computed = quiet(compute_kie, rct=reactant, ts=Series([ts["A"], ts["G"]]), iso=["2", "8", "8"], **OPTIONS)
    ensemble = quiet(compute_kie, rct=reactant, ts=[[ts["A"], ts["G"]]], iso=["2", "8"], **OPTIONS)
    gap = ensemble.rows("transition structure")[1].free_energy
    w = np.array([1.0 / kappa[0], np.exp(gap / rt) / kappa[1]])
    assert computed.kie_tunnel == pytest.approx(np.dot(w / w.sum(), [a.kie_tunnel, g.kie_tunnel]), abs=1e-12)
    assert computed.weight_source == "qrrho"

    # a commitment factor, its inverse and the limits
    committed = quiet(compute_kie, rct=reactant, ts=Series([ts["A"], ts["G"]], commitment=2.0), iso=["2", "8", "8"],
                      **OPTIONS)  # fmt: skip
    assert committed.kie_tunnel == pytest.approx((g.kie_tunnel + 2 * a.kie_tunnel) / 3, abs=1e-13)
    assert committed.kie == pytest.approx((g.kie + 2 * a.kie) / 3, abs=1e-13)
    assert committed.commitment == pytest.approx(2.0)
    assert committed.commitment_for(committed.kie_at(0.5)) == pytest.approx(0.5)
    assert committed.kie_at(1e12) == pytest.approx(a.kie_tunnel)
    with pytest.raises(KinisotInputError, match="outside"):
        committed.commitment_for(1.5)

    # one step is the ordinary calculation
    single = quiet(compute_kie, rct=reactant, ts=Series([ts["A"]]), iso=["2", "8"], **OPTIONS)
    assert single.kie_tunnel == pytest.approx(a.kie_tunnel, abs=1e-14)


def test_series_reference_and_ensemble_steps(shi):
    reactant, ts = shi
    series = Series([[ts["A"], ts["B"]], ts["G"]], free_energies=[0.0, 0.5])
    result = quiet(compute_kie, rct=reactant, ts=series, iso=["2", "8", "8"], reference=["11", "4", "4"], **OPTIONS)
    meta = quiet(compute_kie, rct=reactant, ts=series, iso=["11", "4", "4"], **OPTIONS)
    assert result.kie_tunnel_relative == pytest.approx(result.kie_tunnel / meta.kie_tunnel, abs=1e-13)
    assert result.shares == pytest.approx(meta.shares, abs=1e-14)
    # the first step is a conformer ensemble: its weight carries its conformer sum
    step = result.steps[0]
    rows = step.rows("transition structure")
    rt = GAS_CONSTANT_KCAL * T
    total = sum(r.kappa_light * np.exp(-r.free_energy / rt) for r in rows)
    kappa_g = tunneling_kappa("bell", result.steps[1].other.light.imaginary, T)
    w = np.array([1.0 / total, np.exp(0.5 / rt) / kappa_g])
    assert result.shares == pytest.approx(list(w / w.sum()), abs=1e-13)
    # relative commitment
    two = quiet(compute_kie, rct=reactant, ts=Series([ts["A"], ts["G"]], commitment=1.0), iso=["2", "8", "8"],
                reference=["11", "4", "4"], **OPTIONS)  # fmt: skip
    assert two.commitment_for(two.kie_tunnel_relative, relative=True) == pytest.approx(1.0)
    data = json.loads(two.to_json())
    assert data["series"] and data["commitment"] == pytest.approx(1.0) and len(data["steps"]) == 2
    assert data["reference"]["kie_tunnel_range"] is None
    assert two.summary_row()["kind"] == "KIE (series)"


def test_series_errors(shi):
    reactant, ts = shi
    diene = load_hessian(os.path.join(ROOT, "tests", "data", "gaussian", "DATS.out"))
    with pytest.raises(KinisotInputError, match="different atoms"):
        quiet(compute_kie, rct=reactant, ts=Series([ts["A"], diene]), iso=["2", "8", "10"], **OPTIONS)
    with pytest.raises(KinisotInputError, match="exactly two steps"):
        Series([ts["A"], ts["B"], ts["G"]], commitment=1.0)
    with pytest.raises(KinisotInputError, match="not both"):
        Series([ts["A"], ts["B"]], free_energies=[0, 1], commitment=1.0)
    with pytest.raises(KinisotInputError, match="2 steps but 3"):
        Series([ts["A"], ts["B"]], free_energies=[0, 1, 2])
    with pytest.raises(KinisotInputError, match="alone"):
        compute_kie(rct=reactant, ts=[Series([ts["A"]]), ts["B"]], iso="2")
    with pytest.raises(KinisotInputError, match="free energies"):
        quiet(compute_kie, rct=reactant, ts=Series([ts["A"], ts["B"]]), iso=["2", "8", "8"], weights="user", **OPTIONS)
    with pytest.raises(KinisotInputError, match="isotope labels"):
        quiet(compute_kie, rct=reactant, ts=Series([ts["A"], ts["B"]]), iso=["2", "8"], **OPTIONS)


# ------------------------------------------------------------------ channels from files


def test_channels_from_files(shi):
    reactant, ts = shi
    a, g = pair(reactant, ts["A"]), pair(reactant, ts["G"])
    jobs = [dict(rct=reactant, ts=ts["A"], iso=["2", "8"]), dict(rct=reactant, ts=ts["G"], iso=["2", "8"])]

    # given shares: a measured selectivity, and its inverse
    given = quiet(channels, jobs, shares=[1, 3], **OPTIONS)
    assert isinstance(given, ChannelIsotopeEffect) and given.share_source == "given"
    assert given.kie_tunnel == pytest.approx(channel_kie([a.kie_tunnel, g.kie_tunnel], [1, 3]), abs=1e-13)
    assert given.kie == pytest.approx(channel_kie([a.kie, g.kie], [1, 3]), abs=1e-13)
    assert given.selectivity == pytest.approx(3.0)
    assert given.selectivity_for(given.kie_at([2, 5])) == pytest.approx(2.5)
    assert given.kie_tunnel_range == (given.kie_tunnel, given.kie_tunnel)

    # computed shares from one reactant: the transition structures as a conformer ensemble, exactly
    computed = quiet(channels, jobs, **OPTIONS)
    ensemble = quiet(compute_kie, rct=reactant, ts=[[ts["A"], ts["G"]]], iso=["2", "8"], **OPTIONS)
    assert computed.kie_tunnel == pytest.approx(ensemble.kie_tunnel, abs=1e-13)
    assert computed.kie == pytest.approx(ensemble.kie, abs=1e-13)
    assert computed.shares == pytest.approx([r.population for r in ensemble.rows("transition structure")], abs=1e-12)
    three = quiet(channels, [dict(rct=reactant, ts=[[ts["A"], ts["B"]]], iso=["2", "8"]), jobs[1]], **OPTIONS)
    everything = quiet(compute_kie, rct=reactant, ts=[[ts["A"], ts["B"], ts["G"]]], iso=["2", "8"], **OPTIONS)
    assert three.kie_tunnel == pytest.approx(everything.kie_tunnel, abs=1e-13)

    # barriers equal to the computed free energies give the computed shares; amounts scale them
    gap = ensemble.rows("transition structure")[1].free_energy
    by_barrier = quiet(channels, jobs, barriers=[0.0, gap], **OPTIONS)
    assert by_barrier.kie_tunnel == pytest.approx(computed.kie_tunnel, abs=1e-13)
    assert by_barrier.share_source == "barriers"
    doubled = quiet(channels, jobs, barriers=[0.0, gap], amounts=[1, 2], **OPTIONS)
    assert doubled.shares[1] / doubled.shares[0] == pytest.approx(2 * computed.shares[1] / computed.shares[0])


def test_channels_with_their_own_labels():
    """Kang and Radosevich's TS1B attacks one of two equivalent oxygens: two channels of equal share."""
    case = os.path.join(ROOT, "benchmarks", "nitroarene_phosphetane")
    phno2, ts1b = (load_hessian(os.path.join(case, name + ".out")) for name in ("phno2", "ts1b"))
    options = dict(temperature=393.0, scale=0.9614, tunneling="none", project=True)
    jobs = [
        dict(name="attacked", rct=phno2, ts=ts1b, iso=["12:15N,13:18O", "1:15N,2:18O"], reference=["12:15N", "1:15N"]),
        dict(name="spectator", rct=phno2, ts=ts1b, iso=["12:15N,14:18O", "1:15N,3:18O"], reference=["12:15N", "1:15N"]),
    ]  # fmt: skip
    result = quiet(channels, jobs, shares=[1, 1], **options)
    attacked, spectator = (r.kie_relative for r in result.channels)
    assert result.kie_relative == pytest.approx(2 / (1 / attacked + 1 / spectator), abs=1e-12)
    assert result.kie_relative == pytest.approx(1.0304, abs=1e-4)
    assert result.names == ("attacked", "spectator")
    # computed shares are equal, since both channels are the same structures
    assert quiet(channels, jobs, **options).shares == pytest.approx([0.5, 0.5], abs=1e-14)
    assert result.selectivity_for(result.kie_relative, tunnel=False, relative=True) == pytest.approx(1.0)


def test_channel_errors(shi):
    reactant, ts = shi
    jobs = [dict(rct=reactant, ts=ts["A"], iso=["2", "8"]), dict(rct=reactant, ts=ts["G"], iso=["2", "8"])]
    with pytest.raises(KinisotInputError, match="not both"):
        channels(jobs, shares=[1, 1], barriers=[0, 0])
    with pytest.raises(KinisotInputError, match="needs rct, ts and iso"):
        channels([dict(rct=reactant, ts=ts["A"])])
    with pytest.raises(KinisotInputError, match="for every channel or for none"):
        channels([dict(jobs[0], reference=["11", "4"]), jobs[1]], **OPTIONS)
    with pytest.raises(KinisotInputError, match="2 channels but 1 barriers"):
        channels(jobs, barriers=[0])
    diene, dats = (load_hessian(os.path.join(ROOT, "tests", "data", "gaussian", f)) for f in ("diene.out", "DATS.out"))
    mixed = [jobs[0], dict(rct=diene, ts=dats, iso=["6", "15"])]
    with pytest.raises(KinisotInputError, match="same molecules"):
        quiet(channels, mixed, **OPTIONS)
    assert quiet(channels, mixed, shares=[1, 1], **OPTIONS).kie_tunnel > 0
    with pytest.raises(KinisotInputError, match="shares or barriers"):
        quiet(channels, jobs, weights="user", **OPTIONS)


def test_a_channel_can_be_a_series(shi):
    reactant, ts = shi
    series = Series([ts["A"], ts["G"]], commitment=1.0)
    jobs = [dict(rct=reactant, ts=series, iso=["2", "8", "8"]), dict(rct=reactant, ts=ts["B"], iso=["2", "8"])]
    result = quiet(channels, jobs, shares=[1, 1], **OPTIONS)
    assert isinstance(result.channels[0], SeriesIsotopeEffect)
    expected = channel_kie([result.channels[0].kie_tunnel, pair(reactant, ts["B"]).kie_tunnel], [1, 1])
    assert result.kie_tunnel == pytest.approx(expected, abs=1e-13)
    with pytest.raises(KinisotInputError, match="series"):
        quiet(channels, jobs, **OPTIONS)
    by_barrier = quiet(channels, jobs, barriers=[10.0, 10.0], **OPTIONS)
    assert 0 < by_barrier.shares[0] < 1
    assert quiet(channels, [dict(rct=reactant, ts=Conformers([ts["A"]]), iso=["2", "8"])], **OPTIONS).kie_tunnel == \
        pytest.approx(pair(reactant, ts["A"]).kie_tunnel)  # fmt: skip
