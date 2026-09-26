"""Isotope effects from the Bigeleisen-Mayer equation versus from free-energy differences.

    python scripts/compare_free_energy_route.py

Computes the Claisen KIEs in tests/data/gaussian (B3LYP/6-31G(d), 393 K,
unscaled, no tunnelling) two ways from the same projected harmonic
frequencies:

  BM  the Bigeleisen-Mayer equation as Kinisot evaluates it: the
      translational and rotational partition-function ratios are replaced by
      the product of vibrational frequency ratios (Teller-Redlich), so every
      vibration enters only through ratios of isotopologue frequencies;
  FE  exp(ddG/RT) from rigid-rotor/harmonic-oscillator free energies, as
      when a KIE is taken from the "Sum of electronic and thermal Free
      Energies" a quantum chemistry program prints for each isotopologue:
      translation and rotation from masses and moments of inertia.

With exact harmonic frequencies at an exact stationary point the two are
identical. The table shows what separates them in practice: the six
decimals Gaussian prints, a 0.1 cm-1 error in the isotope shift of one soft
mode, and quasi-harmonic free energies (every mode below 100 cm-1 raised to
100 cm-1). See docs/theory.md, section 6.
"""

import os

import numpy as np

from kinisot import compute_kie, parse_gaussian
from kinisot.hessian import mass_weight
from kinisot.isotopes import substitute
from kinisot.projection import project_external_modes
from kinisot.thermo import (
    BOHR_TO_ANGSTROM,
    BOLTZMANN_CONSTANT,
    ENERGY_AU,
    WAVENUMBER_TO_KELVIN,
    harmonic_frequencies,
)

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
TEMPERATURE = 393.0
SITES = (("C1", "1"), ("C2", "2"), ("O3 (18O)", "3:18O"), ("C4", "4"), ("C5", "5"), ("C6", "6"), ("H7,H8", "7,8"))


def load():
    folder = os.path.join(ROOT, "tests", "data", "gaussian")
    return tuple(parse_gaussian(os.path.join(folder, "claisen_%s.out" % name)) for name in ("gs", "ts"))


def modes(data, label):
    """Masses, real frequencies (cm-1, ascending) and imaginary frequency (or None) of one isotopologue."""
    masses, _ = substitute(data, label, [])
    projected, n_external = project_external_modes(mass_weight(data.hessian, masses), data.positions, masses)
    frequencies = harmonic_frequencies(projected)
    frequencies = np.sort(frequencies[np.argsort(np.abs(frequencies))[n_external:]])
    imaginary = -frequencies[0] if frequencies[0] < -50.0 else None
    return np.asarray(masses), (frequencies[1:] if imaginary else frequencies), imaginary


def _ln_translation_rotation(data, masses):
    """Mass-dependent part of ln(q_trans q_rot): 3/2 ln M + 1/2 ln det(I)."""
    x = np.asarray(data.positions) * BOHR_TO_ANGSTROM
    r = x - masses @ x / masses.sum()
    inertia = sum(m * (v @ v * np.eye(3) - np.outer(v, v)) for m, v in zip(masses, r))
    return 1.5 * np.log(masses.sum()) + 0.5 * np.log(np.linalg.det(inertia))


def _raise(frequencies, floor):
    return np.where(frequencies < floor, floor, frequencies) if floor else frequencies


def _shifted(data, label, shift):
    """Modes of an isotopologue, with ``shift`` = (data, cm-1) added to the heavy one's lowest real mode."""
    masses, real, imaginary = modes(data, label)
    if shift and shift[0] is data and label:
        real = real.copy()
        real[0] += shift[1]
    return masses, real, imaginary


def bigeleisen_mayer(gs, ts, label, temperature=TEMPERATURE, quasi_harmonic=None, shift=None):
    """Semiclassical KIE (with the imaginary-frequency ratio) from the Bigeleisen-Mayer equation."""
    ln_kie = 0.0
    for data, sign in ((gs, 1.0), (ts, -1.0)):
        _, light, imaginary_light = modes(data, "")
        _, heavy, imaginary_heavy = _shifted(data, label, shift)
        ul = WAVENUMBER_TO_KELVIN * _raise(light, quasi_harmonic) / temperature
        uh = WAVENUMBER_TO_KELVIN * _raise(heavy, quasi_harmonic) / temperature
        ln_kie += sign * np.sum(np.log(uh / ul) + (ul - uh) / 2 + np.log1p(-np.exp(-ul)) - np.log1p(-np.exp(-uh)))
        if imaginary_light:
            ln_kie += np.log(imaginary_light / imaginary_heavy)
    return float(np.exp(ln_kie))


def free_energy(data, label, temperature=TEMPERATURE, quasi_harmonic=None, shift=None):
    """E + G_corr in Hartree (mass-independent constants omitted), rigid rotor / harmonic oscillator."""
    masses, real, _ = _shifted(data, label, shift)
    u = WAVENUMBER_TO_KELVIN * _raise(real, quasi_harmonic) / temperature
    ln_q = _ln_translation_rotation(data, masses) + np.sum(-u / 2 - np.log1p(-np.exp(-u)))
    return data.energy - BOLTZMANN_CONSTANT * temperature / ENERGY_AU * ln_q


def free_energy_kie(gs, ts, label, temperature=TEMPERATURE, quasi_harmonic=None, shift=None, decimals=None):
    """exp(ddG/RT) from the four free energies, optionally rounded to ``decimals`` (Gaussian prints 6)."""

    def g(data, iso):
        value = free_energy(data, iso, temperature, quasi_harmonic, shift)
        return round(value, decimals) if decimals is not None else value

    ddg = (g(ts, label) - g(gs, label)) - (g(ts, "") - g(gs, ""))
    return float(np.exp(ddg / (BOLTZMANN_CONSTANT * temperature / ENERGY_AU)))


def product_rule_violation(data, label):
    """ln prod(nu_H/nu_L) over all 3N-6 modes minus its Teller-Redlich value from masses and moments of inertia."""
    ml, light, imaginary_light = modes(data, "")
    mh, heavy, imaginary_heavy = modes(data, label)
    frequencies = np.sum(np.log(heavy / light)) + (np.log(imaginary_heavy / imaginary_light) if imaginary_light else 0)
    masses = _ln_translation_rotation(data, mh) - _ln_translation_rotation(data, ml) - 1.5 * np.sum(np.log(mh / ml))
    return float(frequencies - masses)


def main():
    gs, ts = load()
    kt = BOLTZMANN_CONSTANT * TEMPERATURE / ENERGY_AU
    print("T = %.1f K: 1e-6 Hartree in one free energy changes a KIE by %.1e" % (TEMPERATURE, np.expm1(1e-6 / kt)))
    print("shift: +0.1 cm-1 on the heavy isotopologue's lowest reactant mode (%.1f cm-1)" % modes(gs, "")[1][0])
    header = ("site", "Kinisot", "FE", "FE 6 dp", "BM shift", "FE shift", "BM qh", "FE qh")
    print(("%-9s" + " %9s" * 7) % header)
    for name, label in SITES:
        kinisot = compute_kie(
            rct=gs, ts=ts, iso=label, temperature=TEMPERATURE, scale=1.0, tunneling="none", project=True
        )
        row = (
            kinisot.kie,
            free_energy_kie(gs, ts, label),
            free_energy_kie(gs, ts, label, decimals=6),
            bigeleisen_mayer(gs, ts, label, shift=(gs, 0.1)),
            free_energy_kie(gs, ts, label, shift=(gs, 0.1)),
            bigeleisen_mayer(gs, ts, label, quasi_harmonic=100.0),
            free_energy_kie(gs, ts, label, quasi_harmonic=100.0),
        )
        print(("%-9s" + " %9.5f" * 7) % ((name,) + row))
    print("product-rule violation (ln): " + ", ".join(
        "%s R %+.1e / TS %+.1e" % (name, product_rule_violation(gs, label), product_rule_violation(ts, label))
        for name, label in SITES
    ))  # fmt: skip


if __name__ == "__main__":
    main()
