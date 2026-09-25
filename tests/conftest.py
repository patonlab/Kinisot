#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

try:
    import kinisot
    BASEPATH = os.path.join(kinisot.__path__[0])
except ImportError:
    here = os.path.dirname(os.path.abspath(__file__))
    BASEPATH = os.path.normpath(os.path.join(here, '..', 'kinisot'))


def datapath(path):
    return os.path.join(BASEPATH, 'examples', path)


def synthetic_hessian(masses, freqs_cm, seed=0):
    """Cartesian Hessian (Hartree/Bohr^2) whose mass-weighted eigenvalues give
    exactly ``freqs_cm`` (cm-1, negative for imaginary) with masses ``masses``.

    Used to build small Gaussian-like files with chosen frequencies so that
    mode-counting logic can be tested without real quantum chemistry output.
    """
    from kinisot.Kinisot import HESSIAN_TO_WAVENUMBER_SQ
    n = 3 * len(masses)
    freqs = np.asarray(freqs_cm, dtype=float)
    assert len(freqs) == n, "need 3N frequencies"
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    eigenvalues = np.sign(freqs) * freqs ** 2 / HESSIAN_TO_WAVENUMBER_SQ
    mw = q @ np.diag(eigenvalues) @ q.T
    w = np.repeat(np.sqrt(np.asarray(masses, dtype=float)), 3)
    return mw * w[:, None] * w[None, :]


def write_gaussian_like(path, atomic_numbers, masses, hessian, rotational="12.49 1.87 1.64",
                        level="RB3LYP", basis="6-31G(d)", nimag=0, separator="\\",
                        archive=True, natoms_line=True, wrap=70):
    """Write a minimal file with the parts of a Gaussian freq output Kinisot reads."""
    n = len(masses)
    dof = 3 * n
    lines = [" Entering Gaussian System, Link 0=g16\n"]
    if natoms_line:
        lines.append(" NAtoms=%7d NActive=%7d NUniq=%7d SFac= 1.00D+00 NAtFMM=   60 NAOKFM=F Big=F\n" % (n, n, n))
    for i, (z, m) in enumerate(zip(atomic_numbers, masses)):
        lines.append(" Atom %5d has atomic number %2d and mass %9.5f\n" % (i + 1, z, m))
    lines.append(" Rotational constants (GHZ):     %s\n" % rotational)
    if archive:
        tri = np.asarray(hessian)[np.tril_indices(dof)]
        body = ("1\\1\\GINC-TEST\\Freq\\%s\\%s\\C1H1\\USER\\01-Jan-2026\\0\\\\# freq\\\\title\\\\0,1\\C,0.,0.,0."
                "\\\\Version=ES64L-G16RevC.01\\NImag=%d\\\\" % (level, basis, nimag)
                + ",".join("%.12e" % x for x in tri) + "\\\\0.,0.,0.\\\\@")
        body = body.replace("\\", separator)
        for k in range(0, len(body), wrap):
            lines.append(" " + body[k:k + wrap] + "\n")
    lines.append(" Normal termination of Gaussian 16 at Thu Jan  1 00:00:00 2026.\n")
    with open(path, "w") as handle:
        handle.writelines(lines)
    return str(path)


# Three-atom, non-linear test molecules: six external modes near zero plus
# three vibrations. Atom 1 is 12C, atoms 2-3 are 1H.
MASSES_CH2 = [12.0, 1.00783, 1.00783]
Z_CH2 = [6, 1, 1]
FREQS_MINIMUM = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 900.0, 1400.0, 3000.0]
FREQS_TS = [-500.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 1400.0, 3000.0]


def write_minimum(path, **kwargs):
    return write_gaussian_like(path, Z_CH2, MASSES_CH2, synthetic_hessian(MASSES_CH2, FREQS_MINIMUM), **kwargs)


def write_ts(path, **kwargs):
    kwargs.setdefault("nimag", 1)
    return write_gaussian_like(path, Z_CH2, MASSES_CH2, synthetic_hessian(MASSES_CH2, FREQS_TS, seed=1), **kwargs)
