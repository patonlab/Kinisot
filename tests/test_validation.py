#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Error paths: unusable files, bad labels, inconsistent structures."""

import numpy as np
import pytest

from kinisot import Kinisot
from kinisot.exceptions import KinisotInputError, KinisotParseError, KinisotWarning
from kinisot.Hess_to_Freq import (is_linear, level_of_theory, mass_weight, parse_gaussian,
                                  parse_label, read_hess, substitute)
from conftest import (FREQS_MINIMUM, MASSES_CH2, Z_CH2, datapath, synthetic_hessian,
                      write_gaussian_like, write_minimum, write_ts)

GS = datapath('gaussian/claisen_gs.out')
TS = datapath('gaussian/claisen_ts.out')


def kie(rct, ts, labels, **kw):
    return Kinisot.compute_isotope_effect(rct, ts, None, labels, kw.get('T', 298.15), kw.get('s', 1.0), kw.get('cutoff', 50.0))


# --- parsing -----------------------------------------------------------------

def test_missing_file():
    with pytest.raises(KinisotParseError, match="cannot read"):
        parse_gaussian('/no/such/file.out')


def test_file_without_archive(tmp_path):
    path = write_minimum(tmp_path / 'sp.out', archive=False)
    with pytest.raises(KinisotParseError, match="archive"):
        parse_gaussian(path)
    assert level_of_theory(path) is None


def test_file_without_natoms(tmp_path):
    path = write_minimum(tmp_path / 'x.out', natoms_line=False)
    with pytest.raises(KinisotParseError, match="NAtoms"):
        parse_gaussian(path)


def test_wrong_number_of_force_constants(tmp_path):
    path = write_minimum(tmp_path / 'x.out')
    text = open(path).read().replace('NAtoms=      3', 'NAtoms=      4')  # 3-atom Hessian, 4 atoms claimed
    open(path, 'w').write(text)
    with pytest.raises(KinisotParseError, match="force constants"):
        parse_gaussian(path)


def test_truncated_file(tmp_path):
    with open(GS) as handle:
        lines = handle.readlines()
    path = tmp_path / 'truncated.out'
    path.write_text(''.join(lines[:len(lines) // 2]))
    with pytest.raises(KinisotParseError):
        parse_gaussian(str(path))


def test_windows_archive_separator_and_wrapping(tmp_path):
    path = write_minimum(tmp_path / 'win.out', separator='|', wrap=7)
    data = parse_gaussian(path)
    assert data.level_of_theory == 'RB3LYP/6-31G(d)'
    freqs = Kinisot.harmonic_frequencies(mass_weight(data.hessian, data.masses))
    assert freqs[6:] == pytest.approx([900.0, 1400.0, 3000.0], abs=1e-3)


def test_synthetic_frequencies_round_trip(tmp_path):
    path = write_ts(tmp_path / 'ts.out')
    data = parse_gaussian(path)
    freqs = Kinisot.harmonic_frequencies(mass_weight(data.hessian, data.masses))
    assert freqs == pytest.approx([-500.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 1400.0, 3000.0], abs=1e-3)
    assert data.atomic_numbers == (6, 1, 1)
    assert data.natoms == 3 and not data.linear


def test_parse_gaussian_reads_real_file():
    data = parse_gaussian(GS)
    assert data.natoms == 14
    assert data.level_of_theory == 'RB3LYP/6-31G(d)'
    assert data.masses[0] == pytest.approx(12.0) and data.atomic_numbers[2] == 8
    assert np.allclose(data.hessian, data.hessian.T)
    assert np.allclose(read_hess(GS, '0'), mass_weight(data.hessian, data.masses))


def test_linear_molecule_drops_five_modes(tmp_path):
    masses = [15.99491, 12.0, 15.99491]
    freqs = [0.5, 1.0, 1.5, 2.0, 2.5, 600.0, 600.0, 1300.0, 2300.0]
    path = write_gaussian_like(tmp_path / 'co2.out', [8, 6, 8], masses, synthetic_hessian(masses, freqs),
                               rotational="0.0000000 11.6919157 11.6919157")
    assert parse_gaussian(path).linear and is_linear(path) == 'linear'
    rpfr = Kinisot.calc_rpfr([path], ['0'])
    assert len(rpfr.frequency_wn) == 4
    assert rpfr.discarded_wn[path] == pytest.approx([0.5, 1.0, 1.5, 2.0, 2.5], abs=1e-3)


# --- labels ------------------------------------------------------------------

@pytest.mark.parametrize("label, expected", [
    ('0', []), ('5', [4]), ('7,8', [6, 7]), ('7 8', [6, 7]), (' 1, 3 ', [0, 2]),
])
def test_parse_label(label, expected):
    assert parse_label(label, 14, 'f.out') == expected


@pytest.mark.parametrize("label, message", [
    ('99', "out of range"), ('0,5', "cannot be combined"), ('5,5', "listed twice"),
    ('C5', "not an atom number"), ('-1', "out of range"),
])
def test_bad_labels(label, message):
    with pytest.raises(KinisotInputError, match=message):
        parse_label(label, 14, 'f.out')


def test_out_of_range_label_is_an_error_not_a_silent_noop():
    with pytest.raises(KinisotInputError, match="atom 99 is out of range"):
        read_hess(GS, '99')


def test_substitution_records():
    data = parse_gaussian(GS)
    masses, applied = substitute(data, '5,7')
    assert masses[4] == pytest.approx(13.00335) and masses[6] == pytest.approx(2.0141)
    assert [(s.atom, s.symbol) for s in applied] == [(5, 'C'), (7, 'H')]
    assert str(applied[0]).startswith('C (')


def test_unsupported_element(tmp_path):
    masses = [14.00307, 1.00783, 1.00783]
    path = write_gaussian_like(tmp_path / 'nh2.out', [7, 1, 1], masses, synthetic_hessian(masses, FREQS_MINIMUM))
    with pytest.raises(KinisotInputError, match=r"atom 1 is N.*supported"):
        substitute(parse_gaussian(path), '1')


def test_already_substituted_atom(tmp_path):
    masses = [12.0, 2.0141, 1.00783]
    path = write_gaussian_like(tmp_path / 'chd.out', Z_CH2, masses, synthetic_hessian(masses, FREQS_MINIMUM))
    with pytest.raises(KinisotInputError, match="has mass 2.01410, not the 1H mass"):
        substitute(parse_gaussian(path), '2')


# --- consistency between species ---------------------------------------------

def test_no_substitution_anywhere():
    with pytest.raises(KinisotInputError, match="no isotopic substitution"):
        kie([GS], [TS], ['0', '0'])


def test_different_elements_on_the_two_sides():
    # C5 in the reactant but H7 in the TS: not the same isotopologue
    with pytest.raises(KinisotInputError, match="substitutions differ"):
        kie([GS], [TS], ['5', '7'])


def test_reactant_with_imaginary_frequency():
    with pytest.raises(KinisotInputError, match="imaginary frequency.*given as a reactant"):
        kie([TS], [GS], ['5', '5'])


def test_product_with_imaginary_frequency():
    with pytest.raises(KinisotInputError, match="given as a product"):
        Kinisot.compute_isotope_effect([GS], None, [TS], ['5', '5'])


def test_ts_and_prd_together():
    with pytest.raises(KinisotInputError, match="not both"):
        Kinisot.compute_isotope_effect([GS], [TS], [TS], ['5', '5', '5'])


def test_label_count_mismatch():
    with pytest.raises(KinisotInputError, match="isotope labels"):
        kie([GS], [TS], ['5'])


@pytest.mark.parametrize("kw, message", [({'T': -5.0}, "temperature"), ({'s': 0.0}, "scaling factor")])
def test_nonpositive_parameters(kw, message):
    with pytest.raises(KinisotInputError, match=message):
        kie([GS], [TS], ['5', '5'], **kw)


def test_two_files_with_imaginary_modes_on_one_side(tmp_path):
    ts1, ts2 = write_ts(tmp_path / 'ts1.out'), write_ts(tmp_path / 'ts2.out')
    with pytest.raises(KinisotInputError, match="more than one file"):
        Kinisot.calc_rpfr([ts1, ts2], ['0', '0'])


def test_ts_with_two_imaginary_modes_warns(tmp_path):
    rct = write_minimum(tmp_path / 'rct.out')
    # two imaginary modes, six external modes and one vibration
    freqs = [-500.0, -80.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3000.0]
    ts = write_gaussian_like(tmp_path / 'ts2i.out', Z_CH2, MASSES_CH2, synthetic_hessian(MASSES_CH2, freqs, seed=2), nimag=2)
    with pytest.warns(KinisotWarning, match="2 imaginary frequencies"):
        species, *_, freq_fac = kie([rct], [ts], ['1', '1'])
    assert species[2].im_frequency_wn == pytest.approx(500.0, abs=1e-3)
    # the second imaginary mode was discarded with the external modes, so one external mode (3.0) leaks
    # into the vibrational product: exactly what the warning tells the user
    assert species[2].frequency_wn == pytest.approx([3.0, 3000.0], abs=1e-3)


def test_negative_frequency_never_reaches_the_partition_function(tmp_path):
    freqs = [-500.0, -400.0, -300.0, -200.0, -100.0, -90.0, -80.0, -70.0, 1500.0]
    path = write_gaussian_like(tmp_path / 'bad.out', Z_CH2, MASSES_CH2, synthetic_hessian(MASSES_CH2, freqs, seed=3), nimag=8)
    with pytest.warns(KinisotWarning):
        with pytest.raises(KinisotInputError, match="non-positive frequencies remain"):
            Kinisot.calc_rpfr([path], ['0'])


def test_heavy_ts_mode_below_cutoff_is_reported():
    with pytest.raises(KinisotInputError, match="imaginary frequency beyond the 470.0 cm-1 cutoff"):
        kie([GS], [TS], ['5', '5'], s=0.961, cutoff=470.0)  # light 463.9i at 0.961 scaling


def test_errors_are_value_errors_for_old_callers():
    with pytest.raises(ValueError):
        kie([GS], [TS], ['5', '7'])
    with pytest.raises(ValueError):
        parse_gaussian('/no/such/file.out')
