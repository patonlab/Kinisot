#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Cross-validation against PyQuiverHS on the SI data of Grazioli et al. (ChemRxiv 2026).

Each benchmark case built from that SI keeps the PyQuiverHS input
(``pyquiverhs/<variant>.config``) and its outputs for the same Gaussian files
(``<variant>_range.csv`` from 10 to 1000 K, ``<variant>_<T>K.csv`` at the
experimental temperature). Kinisot is run with the same labels, scaling
factor and imaginary-frequency threshold and must reproduce every
Bigeleisen-Mayer term and total at every temperature. The residual grows as
1/T through the zero-point term, as expected from the two programs' slightly
different isotope masses and constants (see tests/test_pyquiver.py).
"""

import csv
import glob
import os
import re
import warnings

import pytest

from kinisot import KinisotInputError, compute_kie, load_hessian, parse_gaussian

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks")
# case directory: (reactant file, transition structure or product file, "KIE" or "EQE")
CASES = {
    "dihydrophenanthrene": ("gs.log", "ts.log", "KIE"),
    "biaryl_diketone": ("gs.log", "ts.log", "KIE"),
    "metaparacyclophane": ("gs.log", "ts.log", "KIE"),
    "sn2_chloride_methyl_bromide": ("methyl_bromide.log", "ts.log", "KIE"),
    "tetramethylcyclohexane_eie": ("tetramethylcyclohexane.log", "tetramethylcyclohexane.log", "EQE"),
}
# older PyQuiverHS output files use other column names
ALIASES = {
    "EXC": "BM EXC",
    "ZPE": "BM ZPE",
    "Approx. MMI": "BM MMI",
    "Uncorrected (BM) KIE": "BM (total)",
    "Wigner (BM) KIE": "BM (total-Wigner)",
    "Bell (BM) KIE": "BM (total-Bell)",
}


def variants():
    for case in CASES:
        for config in sorted(glob.glob(os.path.join(ROOT, case, "pyquiverhs", "*.config"))):
            yield case, os.path.basename(config)[: -len(".config")]


def read_config(path):
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    scale = float(re.search(r"^scaling\s+(\S+)", text, re.M).group(1))
    threshold = float(re.search(r"^imag_threshold\s+(\S+)", text, re.M).group(1))
    pairs = re.findall(r"^isotopomer\s+\S+\s+(\d+)\s+(\d+)\s+2D", text, re.M)
    return scale, threshold, ",".join(a for a, _ in pairs), ",".join(b for _, b in pairs)


def read_outputs(case, variant):
    rows = []
    for path in sorted(glob.glob(os.path.join(ROOT, case, "pyquiverhs", variant + "_*.csv"))):
        with open(path, encoding="utf-8") as handle:
            rows += [{ALIASES.get(k, k): v for k, v in row.items()} for row in csv.DictReader(handle)]
    return rows


def value(row, key):
    text = (row.get(key) or "").strip()
    return float(text) if text else None  # PyQuiverHS leaves a cell empty when the number overflows


@pytest.mark.parametrize("case,variant", list(variants()))
def test_bigeleisen_mayer_terms_match_pyquiverhs(case, variant):
    reactant_file, other_file, kind = CASES[case]
    config = os.path.join(ROOT, case, "pyquiverhs", variant + ".config")
    scale, threshold, reactant_atoms, other_atoms = read_config(config)
    reactant = load_hessian(os.path.join(ROOT, case, reactant_file))
    other = load_hessian(os.path.join(ROOT, case, other_file))
    rows = read_outputs(case, variant)
    assert len(rows) >= 100
    checked = 0
    for row in rows:
        temperature = float(row["Temperature"])
        for model, key in (("none", "BM (total)"), ("wigner", "BM (total-Wigner)"), ("bell", "BM (total-Bell)")):
            if kind == "EQE" and model != "none":
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                side = {"ts": [other]} if kind == "KIE" else {"prd": [other]}
                try:
                    result = compute_kie(
                        rct=[reactant], iso=[reactant_atoms, other_atoms], temperature=temperature, scale=scale,
                        imag_cutoff=threshold, tunneling=model, project=False, **side,
                    )  # fmt: skip
                except KinisotInputError as error:
                    # Bell's formula has no meaning below the crossover temperature; Kinisot refuses
                    # where PyQuiverHS prints a number (negative for the 45i cm-1 biaryl mode at 10 K)
                    assert model == "bell" and "crossover" in str(error)
                    continue
            # the zero-point term turns small differences in masses and constants into a relative
            # difference that grows as 1/T: largest observed 8.3e-3 K / T (SN2, Bell, 100 K)
            tolerance = 1e-2 / temperature
            expected = value(row, key)
            if expected is not None:
                assert result.kie_tunnel == pytest.approx(expected, rel=tolerance), (variant, temperature, key)
                checked += 1
            if model == "none":
                terms = (("BM ZPE", result.zpe), ("BM EXC", result.exc), ("BM MMI", result.trpf * result.imag_ratio))
                for term, mine in terms:
                    expected = value(row, term)
                    if expected is not None:
                        assert mine == pytest.approx(expected, rel=tolerance), (variant, temperature, term)
    assert checked >= 100


def test_opt_calcall_log_self_check():
    # opt=(calcall,ts) prints a frequency block during the optimization; only the last block is the
    # frequency job's (a 3N-5 slice used to pick up one stray value and trip the self-check)
    data = parse_gaussian(os.path.join(ROOT, "dihydrophenanthrene", "ts.log"))
    assert len(data.program_frequencies) == 3 * len(data.atomic_numbers) - 6
    gs = load_hessian(os.path.join(ROOT, "dihydrophenanthrene", "gs.log"))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = compute_kie(rct=[gs], ts=[data], iso="20,21,23,24", temperature=315.0, scale=0.97)
    assert result.kie == pytest.approx(0.95186027, rel=1e-5)  # PyQuiverHS, SI Table 3
