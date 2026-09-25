"""Deprecated module kept for one minor release (removed in 3.0).

Everything here forwards to the new modules: :mod:`kinisot.api`
(``compute_kie``), :mod:`kinisot.thermo` (constants and Bigeleisen-Mayer
terms), :mod:`kinisot.scaling` and :mod:`kinisot.cli`. The 2.x functions
keep their signatures and return values.
"""

import warnings

from . import __version__
from .api import compute_kie as _compute_kie
from .api import evaluate_isotopologue
from .backends import load_hessian
from .cli import Logger, build_parser, main, write_results  # noqa: F401
from .scaling import choose_scaling_factor, find_scaling_factor  # noqa: F401
from .thermo import (  # noqa: F401
    ATOMIC_MASS_UNIT,
    BOHR_RADIUS,
    BOLTZMANN_CONSTANT,
    ENERGY_AU,
    HESSIAN_TO_WAVENUMBER_SQ,
    PLANCK_CONSTANT,
    SPEED_OF_LIGHT,
    harmonic_frequencies,
    log_excitation_factor,
    log_product_factor,
    log_zpe_factor,
)

__all__ = [
    "compute_isotope_effect",
    "calc_rpfr",
    "calc_product_factor",
    "calc_zpe_factor",
    "calc_excitation_factor",
    "harmonic_frequencies",
    "find_scaling_factor",
    "get_frequency_scaling",
    "Logger",
    "build_parser",
    "main",
    "__version__",
]

# print formatting (names used by 2.x scripts)
space = "   "
dash = "--"
dash_line = space * 17 + " " + dash * 37


def _deprecated(old, new):
    warnings.warn(
        "kinisot.Kinisot.%s is deprecated and will be removed in Kinisot 3.0; use %s" % (old, new),
        DeprecationWarning,
        stacklevel=3,
    )


def calc_product_factor(frequency_wn):
    """Deprecated alias of kinisot.thermo.log_product_factor."""
    return log_product_factor(frequency_wn)


def calc_zpe_factor(frequency_wn, temperature):
    """Deprecated alias of kinisot.thermo.log_zpe_factor."""
    return log_zpe_factor(frequency_wn, temperature)


def calc_excitation_factor(frequency_wn, temperature):
    """Deprecated alias of kinisot.thermo.log_excitation_factor."""
    return log_excitation_factor(frequency_wn, temperature)


class calc_rpfr:
    """Deprecated: the 2.x per-side result object, now built from IsotopologueResult.

    Attributes: ``PF``, ``ZPE``, ``EXC`` (logarithmic sums), ``frequency_wn``
    (kept modes), ``discarded_wn`` and ``im_frequencies`` (per file),
    ``im_frequency_wn`` (set only when the side has a reaction coordinate),
    ``substituted``, ``files`` and ``isomer``.
    """

    def __init__(self, files, isomer, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0, _result=None):
        if _result is None:
            _deprecated("calc_rpfr", "kinisot.compute_kie")
            inputs = [load_hessian(f) for f in files]
            _result = evaluate_isotopologue(inputs, list(isomer), temperature, freq_scale_factor, freq_cutoff, [])
        self.result = _result
        self.files = [s.source for s in _result.species]
        self.isomer = [s.label for s in _result.species]
        self.PF, self.ZPE, self.EXC = _result.log_pf, _result.log_zpe, _result.log_exc
        self.frequency_wn = list(_result.frequencies)
        self.kept_wn = {s.source: list(s.frequencies) for s in _result.species}
        self.discarded_wn = {s.source: list(s.discarded) for s in _result.species}
        self.im_frequencies = {s.source: s.imaginary for s in _result.species if s.imaginary is not None}
        self.substituted = list(_result.substitutions)
        if _result.imaginary is not None:
            self.im_frequency_wn = _result.imaginary


def compute_isotope_effect(rct, ts, prd, label, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):
    """Deprecated: use kinisot.compute_kie(), which returns an IsotopeEffect.

    Returns the 2.x tuple ``(species, ZPE, EXC, TRPF, KIE, KIE_tunnel,
    tunnel_corr, freq_fac)`` where ``species`` lists four calc_rpfr objects
    (reactant light/heavy, TS-or-product light/heavy).
    """
    _deprecated("compute_isotope_effect", "kinisot.compute_kie")
    r = _compute_kie(
        rct, ts, prd, iso=list(label), temperature=temperature, scale=freq_scale_factor, imag_cutoff=freq_cutoff
    )
    species = [calc_rpfr(None, None, _result=iso) for iso in r.species]
    return species, r.zpe, r.exc, r.trpf, r.kie, r.kie_tunnel, r.tunnel_corr, r.imag_ratio


def get_frequency_scaling(files, log):
    """Deprecated: use kinisot.scaling.choose_scaling_factor(). Writes its messages to ``log``."""
    _deprecated("get_frequency_scaling", "kinisot.scaling.choose_scaling_factor")
    choice = choose_scaling_factor([load_hessian(f) for f in files])
    for message in choice.messages:
        log.Write("\n  " + message)
    return choice.factor
