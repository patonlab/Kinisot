#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Kinisot's frequencies must match the ones Gaussian prints.

For every bundled output, the modes Kinisot keeps must agree with the
``Frequencies --`` lines and the modes it discards as external must agree
with Gaussian's (unprojected) ``Low frequencies ---`` values. This guards
the mass weighting, the unit conversion and the mode-dropping rule for any
future backend.
"""

import glob
import os

import numpy as np
import pytest

from kinisot import Kinisot
from kinisot.Hess_to_Freq import mass_weight, parse_gaussian
from conftest import datapath

EXAMPLES = sorted(glob.glob(datapath('gaussian/*.out')))


def gaussian_printed(path):
    freqs, low = [], []
    with open(path) as handle:
        for line in handle:
            s = line.strip()
            if s.startswith('Frequencies --'):
                freqs += [float(x) for x in s.split()[2:]]
            elif s.startswith('Low frequencies ---'):
                low += [float(x) for x in s.split()[3:]]
    return np.array(freqs), np.sort(np.array(low))


@pytest.mark.parametrize('path', EXAMPLES, ids=os.path.basename)
def test_modes_match_gaussian(path):
    data = parse_gaussian(path)
    freqs = Kinisot.harmonic_frequencies(mass_weight(data.hessian, data.masses))
    printed, low = gaussian_printed(path)
    assert len(printed) > 0

    n_imag = int(np.sum(freqs < -50.0))
    assert n_imag in (0, 1)
    n_external = (5 if data.linear else 6) + n_imag

    # Gaussian lists the imaginary mode first among the vibrational frequencies
    kept = np.concatenate([freqs[:n_imag], freqs[n_external:]])
    assert len(kept) == len(printed)
    assert np.abs(kept - printed).max() < 0.05

    # The lowest 6 (+1) unprojected modes are what Gaussian prints as low frequencies
    assert np.abs(np.sort(freqs[:n_external]) - low[:n_external]).max() < 0.5


def test_unit_conversion_constant():
    # 1 Hartree/(amu Bohr^2) corresponds to a wavenumber of about 5140 cm-1
    assert np.sqrt(Kinisot.HESSIAN_TO_WAVENUMBER_SQ) == pytest.approx(5140.49, abs=0.05)
