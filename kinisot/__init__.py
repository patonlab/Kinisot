"""Kinisot: kinetic and equilibrium isotope effects from computed Hessians.

Command line::

    kinisot --rct reactant.out --ts ts.out --iso 5 -t 393

Python::

    from kinisot import compute_isotope_effect
    species, zpe, exc, trpf, kie, kie_tunnel, tunnel_corr, freq_ratio = compute_isotope_effect(
        ["reactant.out"], ["ts.out"], None, ["5", "5"], temperature=393.0, freq_scale_factor=0.961)
"""

__version__ = "2.1.0.dev0"

from .exceptions import KinisotError, KinisotInputError, KinisotParseError, KinisotWarning
from .Hess_to_Freq import FrequencyData, Substitution, mass_weight, parse_gaussian, read_hess, substitute
from .Kinisot import calc_rpfr, compute_isotope_effect, find_scaling_factor, harmonic_frequencies

__all__ = [
    "__version__",
    "compute_isotope_effect",
    "calc_rpfr",
    "harmonic_frequencies",
    "find_scaling_factor",
    "parse_gaussian",
    "substitute",
    "mass_weight",
    "read_hess",
    "FrequencyData",
    "Substitution",
    "KinisotError",
    "KinisotInputError",
    "KinisotParseError",
    "KinisotWarning",
]
