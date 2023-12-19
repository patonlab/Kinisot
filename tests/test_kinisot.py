#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import math
from kinisot import Kinisot
from conftest import datapath

@pytest.mark.parametrize("reactants, ts, iso, scaling, KIE, ZPE, EXC, TRPF, KIE_no_tunnel, KIE_tunnel, parabolic_tunn_corr, freq_fac", [
    # Diels Alder, 298.15K with 0.963 scaling factor
    (['gaussian/dienophile.out', 'gaussian/diene.out'], ['gaussian/DATS.out'], ['0','10', '10'], 0.963, 1.000073, 0.995775, 0.981490, 1.023885, 1.000760, 1.000028, 1.000787)
])
def test_DA(reactants, ts, iso, scaling):
    path = datapath(path)
    temperature = 298.15
    freq_cutoff = 50.0
    kie_val, zpe_val, exc_val, trpf_val, no_tunn_val, tunn_val, tunn_corr_val, freq_val = compute_isotope_effect(rct, ts, prd, iso, temperature, scaling, freq_cutoff)
    precision = 6 
    assert KIE == round(kie_val, precision)
    assert ZPE == round(zpe_val, precision)
    assert EXC == round(exc_val, precision)
    assert TRPF == round(trpf_val, precision)
    assert KIE_no_tunnel == round(no_tunn_val, precision)
    assert KIE_tunnel == round(tunn_val, precision)
    assert parabolic_tunn_corr == round(tunn_corr_val, precision)
    assert freq_fac == round(freq_val, precision)

    
