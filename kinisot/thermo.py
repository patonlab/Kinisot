"""Physical constants and the Bigeleisen-Mayer / tunnelling equations.

Everything here works on arrays of wavenumbers (cm-1) and is independent of
where the Hessian came from. The equations are written out in
docs/theory.md.
"""

import numpy as np

from .exceptions import KinisotInputError

# PHYSICAL CONSTANTS (CODATA 2010; SI apart from the speed of light in cm/s)
PLANCK_CONSTANT = 6.62606957e-34  # J s
BOLTZMANN_CONSTANT = 1.3806488e-23  # J / K
SPEED_OF_LIGHT = 2.99792458e10  # cm / s
ENERGY_AU = 4.35974434e-18  # J
BOHR_RADIUS = 5.2917721092e-11  # m
ATOMIC_MASS_UNIT = 1.660538921e-27  # kg
BOHR_TO_ANGSTROM = BOHR_RADIUS * 1e10

# Multiply a mass-weighted Hessian in Hartree/(amu Bohr^2) by this to get eigenvalues in cm^-2
HESSIAN_TO_WAVENUMBER_SQ = ENERGY_AU / (BOHR_RADIUS**2 * ATOMIC_MASS_UNIT) / ((SPEED_OF_LIGHT * 2 * np.pi) ** 2)

# h c / k: multiply a wavenumber by this and divide by T to get u = h c nu / k T
WAVENUMBER_TO_KELVIN = PLANCK_CONSTANT * SPEED_OF_LIGHT / BOLTZMANN_CONSTANT


def harmonic_frequencies(mw_hessian):
    """Harmonic frequencies in cm-1 (negative for imaginary modes), ascending."""
    eigenvalues = np.linalg.eigvalsh(np.asarray(mw_hessian) * HESSIAN_TO_WAVENUMBER_SQ)
    return np.copysign(np.sqrt(np.abs(eigenvalues)), eigenvalues)


def reduced_energies(frequency_wn, temperature):
    """u = h c nu / k T for an array of wavenumbers."""
    return WAVENUMBER_TO_KELVIN * np.asarray(frequency_wn, dtype=float) / temperature


def log_product_factor(frequency_wn):
    """ln of the product of vibrational temperatures h c nu / k: the Teller-Redlich
    product term of the Bigeleisen-Mayer equation. Temperature independent."""
    return float(np.sum(np.log(WAVENUMBER_TO_KELVIN * np.asarray(frequency_wn, dtype=float))))


def log_zpe_factor(frequency_wn, temperature):
    """ln of the zero-point term: sum(u / 2)."""
    return float(0.5 * np.sum(reduced_energies(frequency_wn, temperature)))


def log_excitation_factor(frequency_wn, temperature):
    """ln of the excitation term: sum(ln(1 - exp(-u)))."""
    u = reduced_energies(frequency_wn, temperature)
    return float(np.sum(np.log1p(-np.exp(-u))))


def crossover_temperature(imaginary_wn):
    """Temperature h c |nu| / 2 pi k below which the Bell correction is not defined."""
    return WAVENUMBER_TO_KELVIN * abs(imaginary_wn) / (2.0 * np.pi)


def bell_correction(imaginary_light, imaginary_heavy, temperature):
    """Ratio Q_t(light) / Q_t(heavy) of Bell infinite-parabola tunnelling factors.

    Q_t = (u/2) / sin(u/2) with u = h c |nu| / k T, so the ratio equals
    (nu_L / nu_H) sin(u_H / 2) / sin(u_L / 2). Raises KinisotInputError below
    the crossover temperature (u >= 2 pi), where the model diverges.
    """
    u_light = reduced_energies(imaginary_light, temperature)
    u_heavy = reduced_energies(imaginary_heavy, temperature)
    if u_light >= 2.0 * np.pi or u_heavy >= 2.0 * np.pi:
        raise KinisotInputError(
            "the Bell tunnelling correction is not defined at %.1f K for an imaginary frequency of "
            "%.1fi cm-1 (crossover temperature %.1f K); use --tunneling wigner or --tunneling none"
            % (temperature, max(imaginary_light, imaginary_heavy), crossover_temperature(imaginary_light))
        )
    return float((imaginary_light / imaginary_heavy) * np.sin(0.5 * u_heavy) / np.sin(0.5 * u_light))


def wigner_correction(imaginary_light, imaginary_heavy, temperature):
    """Ratio of Wigner tunnelling factors Q_t = 1 + u^2 / 24."""
    u_light = reduced_energies(imaginary_light, temperature)
    u_heavy = reduced_energies(imaginary_heavy, temperature)
    return float((1.0 + u_light**2 / 24.0) / (1.0 + u_heavy**2 / 24.0))


TUNNELING_MODELS = {
    "none": lambda light, heavy, temperature: 1.0,
    "bell": bell_correction,
    "wigner": wigner_correction,
}


def tunneling_correction(model, imaginary_light, imaginary_heavy, temperature):
    """Tunnelling correction factor for the KIE, by model name ('none', 'bell', 'wigner')."""
    try:
        function = TUNNELING_MODELS[model]
    except KeyError:
        raise KinisotInputError(
            "unknown tunnelling model %r (choose from %s)" % (model, ", ".join(TUNNELING_MODELS))
        ) from None
    return function(imaginary_light, imaginary_heavy, temperature)
