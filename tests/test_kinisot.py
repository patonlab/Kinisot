#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Characterization tests pinning the numerical behavior of Kinisot.

Golden values were generated from v2.0.2 and cross-checked against the
reference outputs shipped with that release, which they matched to all
printed digits. They were regenerated after the v2.0.3 physical constants
fix (ATOMIC_MASS_UNIT typo), which moved KIE values by <3e-6 relative, and
again for v2.4.0, when the isotope masses moved from the five-decimal values
of Kinisot 1.x/2.x to AME 2020 (largest change 1.2e-6 relative).
"""

import warnings

import pytest
from conftest import datapath

from kinisot import compute_isotope_effect, compute_kie, find_scaling_factor
from kinisot.backends.gaussian import is_linear

REL = 1e-6


def run_kie(reactants, ts, prd, iso, temperature, scaling, freq_cutoff=50.0):
    rct = [datapath(p) for p in reactants]
    ts = [datapath(p) for p in ts] if ts else None
    prd = [datapath(p) for p in prd] if prd else None
    return compute_kie(rct, ts, prd, iso=iso, temperature=temperature, scale=scaling, imag_cutoff=freq_cutoff)


# Columns: name, reactants, ts, prd, iso labels (one per file), T/K, scale factor,
#          V-ratio, ZPE, EXC, TRPF, KIE, 1D-tunn, corr-KIE
CASES = [
    # Claisen rearrangement, 13C/2H KIEs at 393 K with 0.961 scaling
    ("claisen_C5_393K", ["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None,
     ["5", "5"], 393.0, 0.961,
     1.000175990, 0.999077117, 1.000680386, 1.001962452, 1.001895143, 1.000044475, 1.001939702),
    ("claisen_C4_393K", ["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None,
     ["4", "4"], 393.0, 0.961,
     1.012716276, 1.036593905, 1.001969546, 0.978988650, 1.029742457, 1.003156985, 1.032993339),
    ("claisen_H78_393K", ["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None,
     ["7,8", "7,8"], 393.0, 0.961,
     1.007258283, 0.885606164, 1.063923796, 1.006138202, 0.954881853, 1.001815899, 0.956615822),
    # Same substitution, different temperature and no scaling
    ("claisen_C5_298K_unscaled", ["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None,
     ["5", "5"], 298.15, 1.0,
     1.000175990, 0.998734371, 1.001118946, 1.001962452, 1.001990371, 1.000087838, 1.002078384),
    # Diels-Alder with two separate reactant files
    ("DA_multi_rct_C19", ["gaussian/dienophile.out", "gaussian/diene.out"], ["gaussian/DATS.out"], None,
     ["0", "10", "19"], 298.15, 0.963,
     1.000072559, 0.995775433, 0.981489432, 1.023884989, 1.000759502, 1.000027858, 1.000787381),
    ("DA_multi_rct_C15", ["gaussian/dienophile.out", "gaussian/diene.out"], ["gaussian/DATS.out"], None,
     ["0", "6", "15"], 298.15, 0.963,
     1.009797786, 1.001688125, 0.986483642, 1.020845549, 1.018630948, 1.003711415, 1.022411511),
    # Same TS carbon via the pre-formed reactant complex (single reactant file)
    ("DA_single_rct_C15", ["gaussian/DATS_rct.out"], ["gaussian/DATS.out"], None,
     ["15", "15"], 298.15, 0.963,
     1.009797786, 1.003033133, 1.012514270, 0.992401147, 1.017742952, 1.003711415, 1.021520218),
    # Equilibrium isotope effect (--prd path): CD3 axial/equatorial preference
    ("EQE_tmch_290K", ["gaussian/tetramethylcyclohexane.out"], None, ["gaussian/tetramethylcyclohexane.out"],
     ["24,25,26", "28,29,30"], 290.0, 1.0,
     1.000000000, 1.027782603, 1.025380992, 0.985373090, 1.038453903, 1.000000000, 1.038453903),
    ("EQE_tmch_300K", ["gaussian/tetramethylcyclohexane.out"], None, ["gaussian/tetramethylcyclohexane.out"],
     ["24,25,26", "28,29,30"], 300.0, 1.0,
     1.000000000, 1.026844198, 1.025140087, 0.985373090, 1.037262000, 1.000000000, 1.037262000),
]  # fmt: skip


@pytest.mark.parametrize(
    "name, reactants, ts, prd, iso, temperature, scaling, vratio, ZPE, EXC, TRPF, KIE, tunn, corrKIE",
    CASES,
    ids=[c[0] for c in CASES],
)
def test_isotope_effect(
    name, reactants, ts, prd, iso, temperature, scaling, vratio, ZPE, EXC, TRPF, KIE, tunn, corrKIE
):
    r = run_kie(reactants, ts, prd, iso, temperature, scaling)
    assert r.imag_ratio == pytest.approx(vratio, rel=REL)
    assert r.zpe == pytest.approx(ZPE, rel=REL)
    assert r.exc == pytest.approx(EXC, rel=REL)
    assert r.trpf == pytest.approx(TRPF, rel=REL)
    assert r.kie == pytest.approx(KIE, rel=REL)
    assert r.tunnel_corr == pytest.approx(tunn, rel=REL)
    assert r.kie_tunnel == pytest.approx(corrKIE, rel=REL)
    assert r.kind == ("KIE" if ts else "EQE")


def test_deprecated_compute_isotope_effect_matches_compute_kie():
    args = (["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None, ["5", "5"], 393.0, 0.961)
    new = run_kie(*args)
    with pytest.warns(DeprecationWarning, match="compute_isotope_effect is deprecated"):
        species, zpe, exc, trpf, kie, kie_tunnel, tunnel_corr, freq_fac = compute_isotope_effect(
            [datapath(args[0][0])], [datapath(args[1][0])], None, ["5", "5"], 393.0, 0.961
        )
    assert (zpe, exc, trpf, kie, kie_tunnel, tunnel_corr, freq_fac) == (
        new.zpe,
        new.exc,
        new.trpf,
        new.kie,
        new.kie_tunnel,
        new.tunnel_corr,
        new.imag_ratio,
    )
    assert species[2].im_frequency_wn == pytest.approx(463.9, abs=0.05)
    assert not hasattr(species[0], "im_frequency_wn")  # ground state has none
    assert species[1].PF == new.reactant.heavy.log_pf


@pytest.mark.parametrize(
    "level, expected",
    [
        ("RB3LYP/6-31G(d)", 0.977),  # plain lookup, R prefix stripped
        ("RM062X/MG3S", 0.970),  # Gaussian writes M06-2X without hyphen
        ("RCAM-B3LYP/ma-TZVP", 0.976),  # matches the CAM-B3LYP entry...
        ("RCAM-B3LYP/6-31G(d)", None),  # ...but must NOT fall back to B3LYP/6-31G(d)
        ("UB3LYP/6-31G(d)", 0.977),  # U prefix stripped
        ("RHF/3-21G", 0.919),
        ("M06/maug-cc-pVTZ", 0.982),  # exact match, not shadowed by later rows
        ("RM062X/maug-cc-pVTZ", None),  # in Truhlar v3b2 (0.971) but not in v5, which lists maug-cc-pV(T+D)Z
        ("MN15-L/MG3S", 0.977),  # needs the MN15-L alias Kinisot adds to GoodVibes' table
        ("PBE1PBE/MG3S", 0.975),  # Gaussian's name for PBE0, resolved by GoodVibes
        ("B3LYP/6-31G*", 0.977),  # star shorthand resolved by GoodVibes
        ("B3LYP/STO-3G", None),  # basis set not in database
        ("MADEUP/nonsense", None),
    ],
    ids=lambda v: str(v),
)
def test_find_scaling_factor(level, expected):
    factor, ref = find_scaling_factor(level)
    if expected is None:
        assert factor is None and ref is None
    else:
        assert factor == pytest.approx(expected, rel=1e-6)
        assert ref  # a literature reference is returned alongside


def test_is_linear(tmp_path):
    # Linear molecule (CO2-like): first rotational constant is zero
    co2 = tmp_path / "co2.out"
    co2.write_text(" Rotational constants (GHZ):      0.0000000     11.6919157     11.6919157\n")
    assert is_linear(str(co2)) == "linear"
    # Non-linear prolate symmetric top (CH3Cl-like): all three constants nonzero.
    # The pre-v2.0.3 string heuristic misclassified these as linear.
    ch3cl = tmp_path / "ch3cl.out"
    ch3cl.write_text(
        " This molecule is a prolate symmetric top.\n"
        " Rotational constants (GHZ):    152.8000000     13.2900000     13.2900000\n"
    )
    assert is_linear(str(ch3cl)) == "none"


def test_examples_are_nonlinear():
    for f in ["gaussian/claisen_gs.out", "gaussian/DATS.out", "gaussian/tetramethylcyclohexane.out"]:
        assert is_linear(datapath(f)) == "none"


def test_ts_without_imaginary_frequency_raises():
    # A ground-state file passed as the TS must raise a clear error,
    # not crash with NameError (bug fixed in v2.0.3)
    with pytest.raises(ValueError, match="imaginary frequency"):
        run_kie(["gaussian/dienophile.out"], ["gaussian/diene.out"], None, ["0", "0"], 298.15, 1.0)


def test_ts_imaginary_frequency_detected():
    # The uphill TS mode of the Claisen TS at 0.961 scaling
    r = run_kie(["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None, ["5", "5"], 393.0, 0.961)
    assert r.other.light.imaginary == pytest.approx(463.9, abs=0.05)
    assert r.reactant.light.imaginary is None  # ground state has none
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        run_kie(["gaussian/claisen_gs.out"], ["gaussian/claisen_ts.out"], None, ["5", "5"], 393.0, 0.961)
