"""Deprecated module kept for one minor release (removed in 3.0).

The Gaussian parser lives in :mod:`kinisot.backends.gaussian`, the
interchange type in :mod:`kinisot.hessian` and the isotope handling in
:mod:`kinisot.isotopes`. ``FrequencyData`` is the 2.1 name of
``HessianInput``.
"""

from .backends.gaussian import is_linear, level_of_theory, parse_gaussian
from .hessian import HessianInput, mass_weight
from .isotopes import Substitution, parse_label, substitute

FrequencyData = HessianInput

__all__ = [
    "parse_gaussian",
    "level_of_theory",
    "is_linear",
    "read_hess",
    "substitute",
    "parse_label",
    "mass_weight",
    "FrequencyData",
    "HessianInput",
    "Substitution",
]


def read_hess(file, iso):
    """Mass-weighted Hessian of ``file`` with the substitutions in ``iso`` applied."""
    data = parse_gaussian(file)
    masses, _ = substitute(data, iso)
    return mass_weight(data.hessian, masses)
