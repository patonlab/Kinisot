#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import math
from kinisot import Kinisot
from conftest import datapath

@pytest.mark.parametrize("reactants, ts, iso, scaling, ZPE, EXC, TRPF, KIE, tunn, corrKIE", [
    # Diels Alder, 298.15K with 0.963 scaling factor
    (['gaussian/dienophile.out', 'gaussian/diene.out'], ['gaussian/DATS.out'], ['0','10', '19'], 0.963, 0.995775, 0.981490, 1.023885, 1.000760, 1.000028, 1.000787)
])
def test_DA(reactants, ts, iso, scaling, ZPE, EXC, TRPF, KIE, tunn, corrKIE):
    rct = [datapath(path) for path in reactants]
    ts = [datapath(path) for path in ts]
    temperature = 298.15
    freq_cutoff = 50.0
    kie_val, zpe_val, exc_val, trpf_val, no_tunn_val, tunn_val, tunn_corr_val, freq_val = Kinisot.compute_isotope_effect(rct, ts, [], iso, temperature, scaling, freq_cutoff)
    precision = 6 
    assert ZPE == round(zpe_val, precision)
    assert EXC == round(exc_val, precision)
    assert TRPF == round(trpf_val, precision)
    assert KIE == round(no_tunn_val, precision)
    assert tunn == round(tunn_corr_val, precision)
    assert corrKIE == round(tunn_val, precision)
    

    
