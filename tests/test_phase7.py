#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Isotope table and syntax, Eckart projection, Skodje-Truhlar, reference isotopologue, temperature scans."""

import csv
import json
import warnings

import numpy as np
import pytest
from conftest import FREQS_MINIMUM, datapath, synthetic_hessian, write_gaussian_like

from kinisot import HessianInput, KinisotInputError, KinisotWarning, cli, compute_kie, parse_gaussian
from kinisot.hessian import mass_weight
from kinisot.isotope_data import ISOTOPE_MASSES, MOST_ABUNDANT
from kinisot.isotopes import isotope_mass, light_mass, light_masses, parse_label, substitute
from kinisot.projection import external_mode_vectors, project_external_modes
from kinisot.thermo import (
    HARTREE_TO_KCAL_PER_MOL,
    bell_correction,
    harmonic_frequencies,
    skodje_truhlar_correction,
    skodje_truhlar_kappa,
    tunneling_correction,
)

GS, TS = datapath("gaussian/claisen_gs.out"), datapath("gaussian/claisen_ts.out")
KW = dict(temperature=393.0, scale=0.961)


# --- isotope table -----------------------------------------------------------


def test_isotope_table_values():
    assert light_mass("H") == pytest.approx(1.00782503, abs=1e-8) and MOST_ABUNDANT["H"] == 1
    assert light_mass("C") == 12.0 and light_mass("O") == pytest.approx(15.99491462, abs=1e-8)
    assert isotope_mass("C", 13) == pytest.approx(13.00335484, abs=1e-8)
    assert isotope_mass("O", 18) == pytest.approx(17.99915961, abs=1e-8)
    assert isotope_mass("H", 3) == pytest.approx(3.01604928, abs=1e-8)  # tritium is included
    assert isotope_mass("Cl", 37) == pytest.approx(36.96590, abs=1e-5) and MOST_ABUNDANT["Sn"] == 120
    assert len(ISOTOPE_MASSES) >= 80 and all(MOST_ABUNDANT[s] in ISOTOPE_MASSES[s] for s in MOST_ABUNDANT)
    with pytest.raises(KinisotInputError, match="no mass for 19O"):
        isotope_mass("O", 19)
    with pytest.raises(KinisotInputError, match="no isotope data"):
        light_mass("Xx")


def test_parse_label_explicit_forms():
    assert parse_label("5:13C", 14, "f") == [(4, "C", 13, None)]
    assert parse_label("7:D,8:T", 14, "f") == [(6, "H", 2, None), (7, "H", 3, None)]
    assert parse_label("3:18O 5", 14, "f") == [(2, "O", 18, None), (4, None, None, None)]
    assert parse_label("5:13.5", 14, "f") == [(4, None, None, 13.5)]
    with pytest.raises(KinisotInputError, match="cannot read the isotope"):
        parse_label("5:heavy", 14, "f")
    with pytest.raises(KinisotInputError, match="unknown element"):
        parse_label("5:13Xx", 14, "f")


def test_substitute_explicit_and_default():
    data = parse_gaussian(GS)
    masses, applied = substitute(data, "3:17O,5:13C,7:D")
    assert masses[2] == pytest.approx(isotope_mass("O", 17)) and applied[0].isotope == "17O"
    assert masses[4] == pytest.approx(isotope_mass("C", 13)) and masses[6] == pytest.approx(isotope_mass("H", 2))
    assert applied[2].isotope == "2H" and "as 2H" in str(applied[2])
    notes = []
    masses, applied = substitute(data, "3", notes)  # bare oxygen index: 18O, with the one-release note
    assert masses[2] == pytest.approx(isotope_mass("O", 18)) and applied[0].isotope == "18O"
    assert len(notes) == 1 and "18O" in notes[0] and "3:17O" in notes[0]
    assert substitute(data, "5", notes)[1][0].isotope == "13C" and len(notes) == 1
    with pytest.raises(KinisotInputError, match="atom 5 is C but the label asks for 18O"):
        substitute(data, "5:18O")
    with pytest.raises(KinisotInputError, match="already 12C"):
        substitute(data, "5:12C")
    masses, applied = substitute(data, "5:13.5")
    assert masses[4] == 13.5 and applied[0].isotope == "m=13.5"


def test_light_masses_from_table_and_program_check():
    data = parse_gaussian(GS)
    assert light_masses(data) == [light_mass(s) for s in data.symbols]
    std = HessianInput(
        data.hessian, [12.011 if z == 6 else m for z, m in zip(data.atomic_numbers, data.masses)], data.atomic_numbers
    )
    assert light_masses(std) == light_masses(data)  # standard atomic weights pass the consistency check
    with pytest.raises(KinisotInputError, match="do not substitute isotopes in the program"):
        light_masses(HessianInput(data.hessian, [13.0] + list(data.masses[1:]), data.atomic_numbers))


def test_bare_oxygen_warning_reaches_the_result_and_the_cli(tmp_path, monkeypatch):
    with pytest.warns(KinisotWarning, match="bare atom number now means 18O"):
        r = compute_kie(rct=GS, ts=TS, iso="3", **KW)
    explicit = compute_kie(rct=GS, ts=TS, iso="3:18O", **KW)
    assert r.kie == explicit.kie and len(r.warnings) == 2  # one note per file
    assert compute_kie(rct=GS, ts=TS, iso="3:17O", **KW).kie == pytest.approx(1.018834, abs=1e-6)
    monkeypatch.chdir(tmp_path)
    assert cli.main(["--rct", GS, "--ts", TS, "--iso", "3", "-t", "393", "-s", "0.961", "-q"]) == 0
    text = (tmp_path / "Kinisot_output.dat").read_text()
    assert text.count("WARNING") == 2 and "1.039564" in text


# --- projection --------------------------------------------------------------


@pytest.mark.parametrize("name", ["claisen_gs", "claisen_ts", "DATS", "tetramethylcyclohexane"])
def test_projection_zeroes_external_modes_and_matches_gaussian(name):
    data = parse_gaussian(datapath("gaussian/%s.out" % name))
    projected, n_external = project_external_modes(mass_weight(data.hessian, data.masses), data.positions, data.masses)
    assert n_external == 6
    freqs = harmonic_frequencies(projected)
    order = np.argsort(np.abs(freqs))
    assert np.abs(freqs[order[:6]]).max() < 0.01  # external modes are numerically zero
    kept = np.sort(freqs[order[6:]])
    assert np.abs(kept - np.sort(data.program_frequencies)).max() < 0.01  # Gaussian projects too


def test_projection_of_linear_molecule():
    positions = np.array([[0.0, 0.0, -2.2], [0.0, 0.0, 0.0], [0.0, 0.0, 2.2]])  # CO2-like, Bohr
    assert external_mode_vectors(positions, [15.99491, 12.0, 15.99491]).shape == (9, 5)
    assert external_mode_vectors(
        np.array([[0, 0.76, -0.47], [0, 0, 0.12], [0, -0.76, -0.47]]), [1.0, 16.0, 1.0]
    ).shape == (9, 6)


def test_projected_kie_matches_unprojected_on_converged_geometries():
    for iso in ("4", "7,8"):
        plain = compute_kie(rct=GS, ts=TS, iso=iso, **KW)
        proj = compute_kie(rct=GS, ts=TS, iso=iso, project=True, **KW)
        assert proj.project and proj.other.light.species[0].projected
        assert proj.kie_tunnel == pytest.approx(plain.kie_tunnel, rel=1e-6)
        assert max(abs(f) for f in proj.other.light.species[0].discarded) < 0.01
        assert proj.other.light.imaginary == pytest.approx(plain.other.light.imaginary, abs=0.05)
        assert len(proj.other.light.frequencies) == len(plain.other.light.frequencies)
        assert proj.to_dict()["project"] is True and proj.summary_row()["project"] is True


def test_projection_needs_positions(tmp_path):
    data = parse_gaussian(GS)
    no_geometry = HessianInput(data.hessian, data.masses, data.atomic_numbers, source="nogeom")
    with pytest.raises(KinisotInputError, match="needs the geometry"):
        compute_kie(rct=no_geometry, ts=TS, iso="5", project=True, **KW)


def test_projection_rescues_a_noisy_hessian():
    # Add spurious curvature along the six external directions of the real Claisen Hessians, as
    # finite-difference (MLIP) Hessians have. The lowest-six rule then discards a genuine vibration
    # and keeps a spurious one; projection removes exactly the noise and recovers the clean KIE.
    from kinisot.thermo import HESSIAN_TO_WAVENUMBER_SQ

    def noisy(path, spurious_wn):
        data = parse_gaussian(path)
        masses = np.asarray(light_masses(data))
        q = external_mode_vectors(data.positions, masses)
        lam = np.sign(spurious_wn) * np.asarray(spurious_wn, dtype=float) ** 2 / HESSIAN_TO_WAVENUMBER_SQ
        mw = mass_weight(data.hessian, masses) + q @ np.diag(lam) @ q.T
        sqrt_m = np.repeat(np.sqrt(masses), 3)
        hessian = mw * sqrt_m[:, None] * sqrt_m[None, :]
        return HessianInput(hessian, data.masses, data.atomic_numbers, source=path + "+noise", positions=data.positions)

    gs = noisy(GS, [-45.0, -30.0, -20.0, 90.0, 120.0, 200.0])  # below the 50 cm-1 cutoff, above the 70 cm-1 mode
    ts = noisy(TS, [-45.0, -30.0, -20.0, 90.0, 120.0, 200.0])
    clean = compute_kie(rct=GS, ts=TS, iso="4", **KW)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", KinisotWarning)  # the self-check rightly complains about the noisy inputs
        plain = compute_kie(rct=gs, ts=ts, iso="4", **KW)
        projected = compute_kie(rct=gs, ts=ts, iso="4", project=True, **KW)
    # The noise is exactly external only for the light masses; for the heavy isotopologue a small part
    # leaks into the vibrations, so projection recovers the clean KIE to ~1e-5 rather than exactly,
    # still far better than the lowest-six rule.
    assert abs(plain.kie_tunnel - clean.kie_tunnel) > 1e-4  # the lowest-six rule is fooled
    assert projected.kie_tunnel == pytest.approx(clean.kie_tunnel, rel=2e-5)
    assert abs(projected.kie_tunnel - clean.kie_tunnel) < 0.05 * abs(plain.kie_tunnel - clean.kie_tunnel)
    assert max(abs(f) for f in projected.reactant.light.species[0].discarded) < 0.01


# --- Skodje-Truhlar ----------------------------------------------------------


def test_skodje_truhlar_limits_and_branches():
    # high barrier: the correction reduces to Bell
    assert skodje_truhlar_correction(463.9, 458.1, 393.0, 30.0) == pytest.approx(
        bell_correction(463.9, 458.1, 393.0), rel=1e-9
    )
    # low barrier reduces tunnelling; still above 1
    low = skodje_truhlar_correction(1000.0, 990.0, 300.0, 3.0)
    assert 1.0 < low < bell_correction(1000.0, 990.0, 300.0)
    # below the crossover temperature the deep-tunnelling branch is used and stays finite
    assert skodje_truhlar_kappa(1500.0, 200.0, 10.0) > 1.0
    with pytest.raises(KinisotInputError, match="positive barrier"):
        skodje_truhlar_kappa(1000.0, 300.0, -1.0)
    with pytest.raises(KinisotInputError, match="needs the barrier height"):
        tunneling_correction("skodje", 1000.0, 990.0, 300.0)


def test_skodje_in_api_and_cli(tmp_path, monkeypatch):
    gs, ts = parse_gaussian(GS), parse_gaussian(TS)
    assert gs.energy == pytest.approx(-270.5072318) and ts.energy == pytest.approx(-270.461128)
    barrier = (ts.energy - gs.energy) * HARTREE_TO_KCAL_PER_MOL
    r = compute_kie(rct=GS, ts=TS, iso="4", tunneling="skodje", **KW)
    assert r.barrier == pytest.approx(barrier) and r.tunneling == "skodje"
    assert r.tunnel_corr == pytest.approx(
        skodje_truhlar_correction(r.other.light.imaginary, r.other.heavy.imaginary, 393.0, barrier)
    )
    given = compute_kie(rct=GS, ts=TS, iso="4", tunneling="skodje", barrier=5.0, **KW)
    assert given.barrier == 5.0 and given.tunnel_corr < r.tunnel_corr
    assert compute_kie(rct=GS, ts=TS, iso="4", **KW).barrier is None
    no_energy = HessianInput(gs.hessian, gs.masses, gs.atomic_numbers, source="noE", positions=gs.positions)
    with pytest.raises(KinisotInputError, match="needs the barrier height"):
        compute_kie(rct=no_energy, ts=TS, iso="4", tunneling="skodje", **KW)
    monkeypatch.chdir(tmp_path)
    assert (
        cli.main(
            [
                "--rct",
                GS,
                "--ts",
                TS,
                "--iso",
                "4",
                "-t",
                "393",
                "-s",
                "0.961",
                "--tunneling",
                "skodje",
                "--barrier",
                "5",
                "-q",
            ]
        )
        == 0
    )
    assert "tunnelling: skodje (barrier 5.00 kcal/mol)" in (tmp_path / "Kinisot_output.dat").read_text()


# --- reference isotopologue --------------------------------------------------


def test_reference_isotopologue(tmp_path, monkeypatch):
    r = compute_kie(rct=GS, ts=TS, iso="4", reference="5", **KW)
    ref = compute_kie(rct=GS, ts=TS, iso="5", **KW)
    assert r.reference.kie == ref.kie and r.kie_relative == pytest.approx(r.kie / ref.kie)
    assert r.kie_tunnel_relative == pytest.approx(r.kie_tunnel / ref.kie_tunnel)
    assert r.to_dict()["reference"]["kie"] == ref.kie and r.summary_row()["reference_labels"] == "5;5"
    assert compute_kie(rct=GS, ts=TS, iso="4", **KW).kie_relative is None
    monkeypatch.chdir(tmp_path)
    assert (
        cli.main(
            [
                "--rct",
                GS,
                "--ts",
                TS,
                "--iso",
                "4",
                "-t",
                "393",
                "-s",
                "0.961",
                "--reference",
                "5",
                "--csv",
                "r.csv",
                "-q",
            ]
        )
        == 0
    )
    text = (tmp_path / "Kinisot_output.dat").read_text()
    assert "relative to iso @ 5 / 5:" in text and "%.6f" % r.kie_relative in text
    with open(tmp_path / "r.csv", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert float(row["kie_relative"]) == pytest.approx(r.kie_relative)
    with pytest.raises(SystemExit):
        cli.main(["--rct", GS, "--ts", TS, "--iso", "4", "--reference", "5", "--reference", "5", "--reference", "5"])


# --- temperature scans -------------------------------------------------------


def test_parse_temperatures():
    assert cli.parse_temperatures("393") == [393.0]
    assert cli.parse_temperatures("273,298.15, 323") == [273.0, 298.15, 323.0]
    assert cli.parse_temperatures("250:350:50") == [250.0, 300.0, 350.0]
    assert cli.parse_temperatures("300:310:2.5") == [300.0, 302.5, 305.0, 307.5, 310.0]
    for bad in ("", "300:200:10", "300:400", "300:400:0", "abc"):
        with pytest.raises(ValueError):
            cli.parse_temperatures(bad)


def test_temperature_scan_cli(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = [
        "--rct",
        GS,
        "--ts",
        TS,
        "--iso",
        "4",
        "-s",
        "0.961",
        "-t",
        "300:400:50",
        "-q",
        "--json",
        "scan.json",
        "--csv",
        "scan.csv",
    ]
    assert cli.main(args) == 0
    text = (tmp_path / "Kinisot_output.dat").read_text()
    assert text.count("KIE @") == 3 and "Temp = 300.0, 350.0, 400.0K" in text
    assert text.count("Vibrational modes") == 1
    with open(tmp_path / "scan.json") as handle:
        data = json.load(handle)
    assert [d["temperature"] for d in data] == [300.0, 350.0, 400.0]
    with open(tmp_path / "scan.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [float(r["temperature"]) for r in rows] == [300.0, 350.0, 400.0]
    assert float(rows[0]["kie"]) > float(rows[2]["kie"])  # KIEs shrink with temperature
    with pytest.raises(SystemExit):
        cli.main(["--rct", GS, "--ts", TS, "--iso", "4", "-t", "300:200:10"])


def test_clean_inputs_emit_no_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compute_kie(rct=GS, ts=TS, iso="4:13C", project=True, reference="5:13C", tunneling="skodje", **KW)


def test_synthetic_file_masses_are_accepted(tmp_path):
    masses = [12.0, 1.00783, 1.00783]
    path = write_gaussian_like(tmp_path / "m.out", [6, 1, 1], masses, synthetic_hessian(masses, FREQS_MINIMUM))
    assert light_masses(parse_gaussian(path)) == pytest.approx([12.0, 1.00782503, 1.00782503], abs=1e-7)
