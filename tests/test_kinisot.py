#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Characterization tests pinning the numerical behavior of Kinisot.

Golden values were generated from v2.0.2 and cross-checked against the
reference outputs shipped in kinisot/examples/gaussian/ (claisen_kinisot.dat,
tetramethylcyclohexane_kinisot.dat, comparison.dat), which they matched to
all printed digits. They were then regenerated after the v2.0.3 physical
constants fix (ATOMIC_MASS_UNIT typo), which moved KIE values by <3e-6
relative; the shipped .dat files still agree to their printed 6 decimals.
"""

import pytest

import kinisot
from conftest import datapath

REL = 1e-6

# ORCA end-to-end golden: n-pentane TT vs GG conformer EQE, D at atom 6,
# 298.15 K, unscaled. Pinned from the Kinisot+GoodVibes computation with
# principal-isotope mass normalization.
EQE_PENTANE_TT_GG_D6 = 1.007730117


def run_kie(reactants, ts, prd, iso, temperature, scaling, freq_cutoff=50.0):
    rct = [datapath(p) for p in reactants]
    ts = [datapath(p) for p in ts] if ts else None
    prd = [datapath(p) for p in prd] if prd else None
    return kinisot.compute_isotope_effect(
        rct, ts, prd, iso, temperature, scaling, freq_cutoff
    )


# Columns: name, reactants, ts, prd, iso labels (one per file), T/K, scale factor,
#          V-ratio, ZPE, EXC, TRPF, KIE, 1D-tunn, corr-KIE
CASES = [
    # Claisen rearrangement, 13C/2H KIEs at 393 K with 0.961 scaling
    (
        "claisen_C5_393K",
        ["gaussian/claisen_gs.out"],
        ["gaussian/claisen_ts.out"],
        None,
        ["5", "5"],
        393.0,
        0.961,
        1.000175989,
        0.999077127,
        1.000680383,
        1.001962438,
        1.001895135,
        1.000044475,
        1.001939694,
    ),
    (
        "claisen_C4_393K",
        ["gaussian/claisen_gs.out"],
        ["gaussian/claisen_ts.out"],
        None,
        ["4", "4"],
        393.0,
        0.961,
        1.012716215,
        1.036593731,
        1.001969536,
        0.978988748,
        1.029742315,
        1.003156971,
        1.032993182,
    ),
    (
        "claisen_H78_393K",
        ["gaussian/claisen_gs.out"],
        ["gaussian/claisen_ts.out"],
        None,
        ["7,8", "7,8"],
        393.0,
        0.961,
        1.007258229,
        0.885607175,
        1.063923244,
        1.006138148,
        0.954882344,
        1.001815886,
        0.956616301,
    ),
    # Same substitution, different temperature and no scaling
    (
        "claisen_C5_298K_unscaled",
        ["gaussian/claisen_gs.out"],
        ["gaussian/claisen_ts.out"],
        None,
        ["5", "5"],
        298.15,
        1.0,
        1.000175989,
        0.998734385,
        1.001118940,
        1.001962438,
        1.001990364,
        1.000087837,
        1.002078376,
    ),
    # Diels-Alder with two separate reactant files
    (
        "DA_multi_rct_C19",
        ["gaussian/dienophile.out", "gaussian/diene.out"],
        ["gaussian/DATS.out"],
        None,
        ["0", "10", "19"],
        298.15,
        0.963,
        1.000072559,
        0.995775463,
        0.981489550,
        1.023884832,
        1.000759499,
        1.000027858,
        1.000787377,
    ),
    (
        "DA_multi_rct_C15",
        ["gaussian/dienophile.out", "gaussian/diene.out"],
        ["gaussian/DATS.out"],
        None,
        ["0", "6", "15"],
        298.15,
        0.963,
        1.009797738,
        1.001688142,
        0.986483717,
        1.020845417,
        1.018630864,
        1.003711397,
        1.022411407,
    ),
    # Same TS carbon via the pre-formed reactant complex (single reactant file)
    (
        "DA_single_rct_C15",
        ["gaussian/DATS_rct.out"],
        ["gaussian/DATS.out"],
        None,
        ["15", "15"],
        298.15,
        0.963,
        1.009797738,
        1.003033140,
        1.012514194,
        0.992401183,
        1.017742871,
        1.003711397,
        1.021520119,
    ),
    # Equilibrium isotope effect (--prd path): CD3 axial/equatorial preference
    (
        "EQE_tmch_290K",
        ["gaussian/tetramethylcyclohexane.out"],
        None,
        ["gaussian/tetramethylcyclohexane.out"],
        ["24,25,26", "28,29,30"],
        290.0,
        1.0,
        1.000000000,
        1.027782351,
        1.025380776,
        0.985373194,
        1.038453539,
        1.000000000,
        1.038453539,
    ),
    (
        "EQE_tmch_300K",
        ["gaussian/tetramethylcyclohexane.out"],
        None,
        ["gaussian/tetramethylcyclohexane.out"],
        ["24,25,26", "28,29,30"],
        300.0,
        1.0,
        1.000000000,
        1.026843954,
        1.025139873,
        0.985373194,
        1.037261647,
        1.000000000,
        1.037261647,
    ),
]


@pytest.mark.parametrize(
    "name, reactants, ts, prd, iso, temperature, scaling, vratio, ZPE, EXC, TRPF, KIE, tunn, corrKIE",
    CASES,
    ids=[c[0] for c in CASES],
)
def test_isotope_effect(
    name,
    reactants,
    ts,
    prd,
    iso,
    temperature,
    scaling,
    vratio,
    ZPE,
    EXC,
    TRPF,
    KIE,
    tunn,
    corrKIE,
):
    r = run_kie(reactants, ts, prd, iso, temperature, scaling)
    assert r.v_ratio == pytest.approx(vratio, rel=REL)
    assert r.zpe == pytest.approx(ZPE, rel=REL)
    assert r.exc == pytest.approx(EXC, rel=REL)
    assert r.trpf == pytest.approx(TRPF, rel=REL)
    assert r.kie == pytest.approx(KIE, rel=REL)
    assert r.bell_correction == pytest.approx(tunn, rel=REL)
    assert r.kie_bell == pytest.approx(corrKIE, rel=REL)


@pytest.mark.parametrize(
    "level, expected",
    [
        ("RB3LYP/6-31G(d)", 0.977),  # plain lookup, R prefix stripped
        ("RM062X/MG3S", 0.970),  # Gaussian writes M06-2X without hyphen
        ("M06-2X/MG3S", 0.970),  # hyphenated (ORCA-style) alias
        ("RCAM-B3LYP/ma-TZVP", 0.976),  # matches the CAM-B3LYP entry...
        ("RCAM-B3LYP/6-31G(d)", None),  # ...but must NOT fall back to B3LYP/6-31G(d)
        ("UB3LYP/6-31G(d)", 0.977),  # U prefix stripped
        ("RHF/3-21G", 0.919),
        ("M06/maug-cc-pVTZ", 0.982),
        ("B3LYP/STO-3G", None),  # basis set not in database
        ("MADEUP/nonsense", None),
    ],
    ids=lambda v: str(v),
)
def test_find_scaling_factor(level, expected):
    factor, ref = kinisot.find_scaling_factor(level)
    if expected is None:
        assert factor is None and ref is None
    else:
        assert factor == pytest.approx(expected, rel=1e-6)
        assert ref  # a literature reference is returned alongside


def test_is_linear(tmp_path):
    from kinisot.Hess_to_Freq import is_linear

    # Linearity comes from the point group detected by GoodVibes' parsers.
    # Linear molecule (CO2-like): D*H point group
    co2 = tmp_path / "co2.log"
    co2.write_text(
        " Gaussian 16\n"
        " Full point group                 D*H\n"
        " Rotational constants (GHZ):      0.0000000     11.6919157     11.6919157\n"
    )
    assert is_linear(str(co2)) == "linear"
    # Non-linear prolate symmetric top (CH3Cl-like, C3V): the pre-v2.0.3
    # string heuristic misclassified these as linear.
    ch3cl = tmp_path / "ch3cl.log"
    ch3cl.write_text(
        " Gaussian 16\n"
        " Full point group                 C3V\n"
        " This molecule is a prolate symmetric top.\n"
        " Rotational constants (GHZ):    152.8000000     13.2900000     13.2900000\n"
    )
    assert is_linear(str(ch3cl)) == "none"


def test_examples_are_nonlinear():
    from kinisot.Hess_to_Freq import is_linear

    for f in [
        "gaussian/claisen_gs.out",
        "gaussian/DATS.out",
        "gaussian/tetramethylcyclohexane.out",
    ]:
        assert is_linear(datapath(f)) == "none"


def test_ts_without_imaginary_frequency_raises():
    # A ground-state file passed as the TS must raise a clear error,
    # not crash with NameError (bug fixed in v2.0.3)
    with pytest.raises(ValueError, match="imaginary frequency"):
        run_kie(
            ["gaussian/dienophile.out"],
            ["gaussian/diene.out"],
            None,
            ["0", "0"],
            298.15,
            1.0,
        )


def test_orca_eqe_pentane_conformers():
    # End-to-end ORCA support via GoodVibes parse_hessian: conformational
    # equilibrium isotope effect for n-pentane TT vs GG, with one methyl H
    # replaced by D in each. Exercises the ORCA .hess path (companion file
    # next to the .out) including ORCA's average-atomic-mass convention.
    r = run_kie(
        ["orca/pentane_TT.out"], None, ["orca/pentane_GG.out"], ["6", "6"], 298.15, 1.0
    )
    # An H/D conformer EQE must be a small effect close to unity
    assert 0.9 < r.kie < 1.1
    assert r.is_eqe
    assert r.v_ratio == 1.0 and r.bell_correction == 1.0  # no TS: no tunneling terms
    assert r.wigner_correction == 1.0
    assert r.kie == r.kie_bell == r.kie_wigner
    # Pinned at first computation (Kinisot + GoodVibes parse_hessian)
    assert r.kie == pytest.approx(EQE_PENTANE_TT_GG_D6, rel=1e-6)


def test_orca_hess_direct_path():
    # Passing the .hess file itself must give identical numbers to the .out
    via_out = run_kie(
        ["orca/pentane_TT.out"], None, ["orca/pentane_GG.out"], ["6", "6"], 298.15, 1.0
    )
    via_hess = run_kie(
        ["orca/pentane_TT.hess"],
        None,
        ["orca/pentane_GG.hess"],
        ["6", "6"],
        298.15,
        1.0,
    )
    assert via_hess.kie == pytest.approx(via_out.kie, rel=1e-12)


def test_ts_imaginary_frequency_detected():
    # The uphill TS mode of the Claisen TS at 0.961 scaling, from claisen_kinisot.dat
    r = run_kie(
        ["gaussian/claisen_gs.out"],
        ["gaussian/claisen_ts.out"],
        None,
        ["5", "5"],
        393.0,
        0.961,
    )
    assert r.rpfrs[2].im_frequency_wn == pytest.approx(463.9, abs=0.05)
    assert r.rpfrs[0].im_frequency_wn is None  # ground state has none


def test_wigner_between_uncorrected_and_bell():
    # For a normal KIE the Wigner correction is smaller than Bell's
    # infinite parabola but still > 1
    r = run_kie(
        ["gaussian/claisen_gs.out"],
        ["gaussian/claisen_ts.out"],
        None,
        ["4", "4"],
        393.0,
        0.961,
    )
    assert 1.0 < r.wigner_correction < r.bell_correction
    assert r.kie < r.kie_wigner < r.kie_bell


def test_iso_out_of_range_raises():
    # The Claisen system has 14 atoms; atom 99 must raise, naming the file
    with pytest.raises(ValueError, match="out of range.*claisen_gs"):
        run_kie(
            ["gaussian/claisen_gs.out"],
            ["gaussian/claisen_ts.out"],
            None,
            ["99", "99"],
            298.15,
            1.0,
        )


def test_iso_not_a_number_raises():
    with pytest.raises(ValueError, match="Invalid atom number"):
        run_kie(
            ["gaussian/claisen_gs.out"],
            ["gaussian/claisen_ts.out"],
            None,
            ["C5", "C5"],
            298.15,
            1.0,
        )


def test_neither_ts_nor_prd_raises():
    with pytest.raises(ValueError, match="TS .*or a product"):
        kinisot.compute_isotope_effect(
            [datapath("gaussian/claisen_gs.out")],
            None,
            None,
            ["5", "5"],
            298.15,
            1.0,
            50.0,
        )


def test_multi_isotopologue_matches_single():
    # compute_isotope_effects shares the light-species RPFRs; results must be
    # identical to one-at-a-time calls
    rct = [datapath("gaussian/claisen_gs.out")]
    ts = [datapath("gaussian/claisen_ts.out")]
    labels = [["5", "5"], ["4", "4"], ["7,8", "7,8"]]
    together = kinisot.compute_isotope_effects(rct, ts, None, labels, 393.0, 0.961)
    for label, effect in zip(labels, together):
        single = kinisot.compute_isotope_effect(rct, ts, None, label, 393.0, 0.961)
        assert effect.kie == pytest.approx(single.kie, rel=1e-12)
        assert effect.kie_bell == pytest.approx(single.kie_bell, rel=1e-12)


def test_explicit_isotope_syntax():
    # 5:13C is the same as the default 5; 5:14C differs; 3:2D is rejected (atom 3 is C)
    rct = [datapath("gaussian/claisen_gs.out")]
    ts = [datapath("gaussian/claisen_ts.out")]
    default = kinisot.compute_isotope_effect(rct, ts, None, ["5", "5"], 393.0, 0.961)
    explicit = kinisot.compute_isotope_effect(
        rct, ts, None, ["5:13C", "5:13C"], 393.0, 0.961
    )
    assert explicit.kie == pytest.approx(default.kie, rel=1e-12)
    c14 = kinisot.compute_isotope_effect(
        rct, ts, None, ["5:14C", "5:14C"], 393.0, 0.961
    )
    assert c14.kie != pytest.approx(default.kie, rel=1e-9)
    with pytest.raises(ValueError, match="cannot replace"):
        kinisot.compute_isotope_effect(rct, ts, None, ["5:2D", "5:2D"], 393.0, 0.961)
    with pytest.raises(ValueError, match="Unknown isotope"):
        kinisot.compute_isotope_effect(rct, ts, None, ["5:99X", "5:99X"], 393.0, 0.961)


def test_legacy_kinisot_module_shim():
    # The historical kinisot.Kinisot tuple API still works, with a DeprecationWarning
    from kinisot import Kinisot

    args = (
        [datapath("gaussian/claisen_gs.out")],
        [datapath("gaussian/claisen_ts.out")],
        None,
        ["5", "5"],
        393.0,
        0.961,
        50.0,
    )
    new = kinisot.compute_isotope_effect(*args)
    with pytest.warns(DeprecationWarning):
        legacy = Kinisot.compute_isotope_effect(*args)
    rpfrs, zpe, exc, trpf, kie, kie_bell, bell, v_ratio = legacy
    assert kie == pytest.approx(new.kie, rel=1e-12)
    assert kie_bell == pytest.approx(new.kie_bell, rel=1e-12)
    assert v_ratio == pytest.approx(new.v_ratio, rel=1e-12)
    # the legacy calc_rpfr class only has im_frequency_wn when a mode was found
    with pytest.warns(DeprecationWarning):
        gs_rpfr = Kinisot.calc_rpfr(
            [datapath("gaussian/claisen_gs.out")], ["0"], 393.0, 0.961
        )
    assert not hasattr(gs_rpfr, "im_frequency_wn")


def test_normalize_principal_masses():
    from kinisot.Hess_to_Freq import normalize_principal_masses

    # ORCA-style abundance-averaged masses map onto principal isotope masses
    assert normalize_principal_masses([12.011, 1.008, 15.999]) == pytest.approx(
        [12.00000, 1.00783, 15.99491]
    )
    # Gaussian-style principal masses are unchanged
    assert normalize_principal_masses([12.00000, 1.00783, 15.99491]) == pytest.approx(
        [12.00000, 1.00783, 15.99491]
    )
    # nitrogen normalizes too; elements without isotope data keep the program's value
    assert normalize_principal_masses([35.453, 14.0067]) == pytest.approx(
        [35.453, 14.00307]
    )


def test_substitute_isotopes_windows():
    from kinisot.Hess_to_Freq import substitute_isotopes

    # Gaussian isotopic masses and ORCA average atomic masses both substitute
    gaussian_masses = [12.00000, 1.00783, 15.99491]
    orca_masses = [12.011, 1.008, 15.999]
    for masses in (gaussian_masses, orca_masses):
        out = substitute_isotopes(masses, "1,2,3")
        assert out == pytest.approx([13.00335, 2.0141, 16.9991])
    # '0' leaves everything alone
    assert substitute_isotopes(orca_masses, "0") == orca_masses
    # nitrogen now substitutes to 15N by default
    assert substitute_isotopes([14.00307], "1") == pytest.approx([15.00011])
    # elements without isotope data raise instead of being silently skipped
    with pytest.raises(ValueError, match="No isotope substitution"):
        substitute_isotopes([35.453], "1")
