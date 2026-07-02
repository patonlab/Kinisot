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
from kinisot import Kinisot
from conftest import datapath

REL = 1e-6


def run_kie(reactants, ts, prd, iso, temperature, scaling, freq_cutoff=50.0):
    rct = [datapath(p) for p in reactants]
    ts = [datapath(p) for p in ts] if ts else None
    prd = [datapath(p) for p in prd] if prd else None
    return Kinisot.compute_isotope_effect(rct, ts, prd, iso, temperature, scaling, freq_cutoff)


# Columns: name, reactants, ts, prd, iso labels (one per file), T/K, scale factor,
#          V-ratio, ZPE, EXC, TRPF, KIE, 1D-tunn, corr-KIE
CASES = [
    # Claisen rearrangement, 13C/2H KIEs at 393 K with 0.961 scaling
    ("claisen_C5_393K", ['gaussian/claisen_gs.out'], ['gaussian/claisen_ts.out'], None,
     ['5', '5'], 393.0, 0.961,
     1.000175989, 0.999077127, 1.000680383, 1.001962438, 1.001895135, 1.000044475, 1.001939694),
    ("claisen_C4_393K", ['gaussian/claisen_gs.out'], ['gaussian/claisen_ts.out'], None,
     ['4', '4'], 393.0, 0.961,
     1.012716215, 1.036593731, 1.001969536, 0.978988748, 1.029742315, 1.003156971, 1.032993182),
    ("claisen_H78_393K", ['gaussian/claisen_gs.out'], ['gaussian/claisen_ts.out'], None,
     ['7,8', '7,8'], 393.0, 0.961,
     1.007258229, 0.885607175, 1.063923244, 1.006138148, 0.954882344, 1.001815886, 0.956616301),
    # Same substitution, different temperature and no scaling
    ("claisen_C5_298K_unscaled", ['gaussian/claisen_gs.out'], ['gaussian/claisen_ts.out'], None,
     ['5', '5'], 298.15, 1.0,
     1.000175989, 0.998734385, 1.001118940, 1.001962438, 1.001990364, 1.000087837, 1.002078376),
    # Diels-Alder with two separate reactant files
    ("DA_multi_rct_C19", ['gaussian/dienophile.out', 'gaussian/diene.out'], ['gaussian/DATS.out'], None,
     ['0', '10', '19'], 298.15, 0.963,
     1.000072559, 0.995775463, 0.981489550, 1.023884832, 1.000759499, 1.000027858, 1.000787377),
    ("DA_multi_rct_C15", ['gaussian/dienophile.out', 'gaussian/diene.out'], ['gaussian/DATS.out'], None,
     ['0', '6', '15'], 298.15, 0.963,
     1.009797738, 1.001688142, 0.986483717, 1.020845417, 1.018630864, 1.003711397, 1.022411407),
    # Same TS carbon via the pre-formed reactant complex (single reactant file)
    ("DA_single_rct_C15", ['gaussian/DATS_rct.out'], ['gaussian/DATS.out'], None,
     ['15', '15'], 298.15, 0.963,
     1.009797738, 1.003033140, 1.012514194, 0.992401183, 1.017742871, 1.003711397, 1.021520119),
    # Equilibrium isotope effect (--prd path): CD3 axial/equatorial preference
    ("EQE_tmch_290K", ['gaussian/tetramethylcyclohexane.out'], None, ['gaussian/tetramethylcyclohexane.out'],
     ['24,25,26', '28,29,30'], 290.0, 1.0,
     1.000000000, 1.027782351, 1.025380776, 0.985373194, 1.038453539, 1.000000000, 1.038453539),
    ("EQE_tmch_300K", ['gaussian/tetramethylcyclohexane.out'], None, ['gaussian/tetramethylcyclohexane.out'],
     ['24,25,26', '28,29,30'], 300.0, 1.0,
     1.000000000, 1.026843954, 1.025139873, 0.985373194, 1.037261647, 1.000000000, 1.037261647),
]


@pytest.mark.parametrize(
    "name, reactants, ts, prd, iso, temperature, scaling, vratio, ZPE, EXC, TRPF, KIE, tunn, corrKIE",
    CASES, ids=[c[0] for c in CASES])
def test_isotope_effect(name, reactants, ts, prd, iso, temperature, scaling,
                        vratio, ZPE, EXC, TRPF, KIE, tunn, corrKIE):
    (rpfrs, zpe_val, exc_val, trpf_val, kie_val, corr_kie_val,
     tunn_corr_val, freq_fac) = run_kie(reactants, ts, prd, iso, temperature, scaling)
    assert freq_fac == pytest.approx(vratio, rel=REL)
    assert zpe_val == pytest.approx(ZPE, rel=REL)
    assert exc_val == pytest.approx(EXC, rel=REL)
    assert trpf_val == pytest.approx(TRPF, rel=REL)
    assert kie_val == pytest.approx(KIE, rel=REL)
    assert tunn_corr_val == pytest.approx(tunn, rel=REL)
    assert corr_kie_val == pytest.approx(corrKIE, rel=REL)


@pytest.mark.parametrize("level, expected", [
    ("RB3LYP/6-31G(d)", 0.977),          # plain lookup, R prefix stripped
    ("RM062X/MG3S", 0.970),              # Gaussian writes M06-2X without hyphen
    ("RCAM-B3LYP/ma-TZVP", 0.976),       # matches the CAM-B3LYP entry...
    ("RCAM-B3LYP/6-31G(d)", None),       # ...but must NOT fall back to B3LYP/6-31G(d)
    ("UB3LYP/6-31G(d)", 0.977),          # U prefix stripped
    ("RHF/3-21G", 0.919),
    ("M06/maug-cc-pVTZ", 0.982),         # exact match, not shadowed by later rows
    ("RM062X/maug-cc-pVTZ", 0.971),      # row was mislabeled M06/maug-cc-pVTZ before v2.0.3
    ("B3LYP/STO-3G", None),              # basis set not in database
    ("MADEUP/nonsense", None),
], ids=lambda v: str(v))
def test_find_scaling_factor(level, expected):
    factor, ref = Kinisot.find_scaling_factor(level)
    if expected is None:
        assert factor is None and ref is None
    else:
        assert factor == pytest.approx(expected, rel=1e-6)
        assert ref  # a literature reference is returned alongside


def test_ts_without_imaginary_frequency_raises():
    # A ground-state file passed as the TS must raise a clear error,
    # not crash with NameError (bug fixed in v2.0.3)
    with pytest.raises(ValueError, match="imaginary frequency"):
        run_kie(['gaussian/dienophile.out'], ['gaussian/diene.out'], None,
                ['0', '0'], 298.15, 1.0)


def test_ts_imaginary_frequency_detected():
    # The uphill TS mode of the Claisen TS at 0.961 scaling, from claisen_kinisot.dat
    rpfrs, *_ = run_kie(['gaussian/claisen_gs.out'], ['gaussian/claisen_ts.out'], None,
                        ['5', '5'], 393.0, 0.961)
    assert rpfrs[2].im_frequency_wn == pytest.approx(463.9, abs=0.05)
    assert not hasattr(rpfrs[0], "im_frequency_wn")  # ground state has none
