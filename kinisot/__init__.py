"""Kinisot: kinetic and equilibrium isotope effects from computed Hessians.

Command line::

    kinisot --rct reactant.out --ts ts.out --iso 5 -t 393

Python::

    from kinisot import compute_kie
    result = compute_kie(rct="reactant.out", ts="ts.out", iso="5", temperature=393.0, scale=0.961)
    result.kie_tunnel, result.zpe, result.exc, result.trpf, result.to_dict()
"""

__version__ = "2.5.0.dev0"

from .api import IsotopeEffect, IsotopologueResult, SideResult, SpeciesResult, compute_kie
from .backends import load_hessian
from .backends.ase import build_calculator, hessian_for_geometry, hessian_from_calculator, save_hessian_json
from .backends.gaussian import parse_gaussian
from .exceptions import KinisotError, KinisotInputError, KinisotParseError, KinisotWarning
from .hessian import HessianInput, mass_weight
from .isotopes import Substitution, isotope_mass, light_mass, substitute
from .Kinisot import compute_isotope_effect  # deprecated, removed in 3.0
from .projection import project_external_modes
from .scaling import ScalingChoice, choose_scaling_factor, find_scaling_factor
from .thermo import harmonic_frequencies

__all__ = [
    "__version__",
    "compute_kie",
    "IsotopeEffect",
    "SideResult",
    "IsotopologueResult",
    "SpeciesResult",
    "HessianInput",
    "load_hessian",
    "parse_gaussian",
    "hessian_from_calculator",
    "hessian_for_geometry",
    "save_hessian_json",
    "build_calculator",
    "substitute",
    "Substitution",
    "isotope_mass",
    "light_mass",
    "project_external_modes",
    "mass_weight",
    "harmonic_frequencies",
    "find_scaling_factor",
    "choose_scaling_factor",
    "ScalingChoice",
    "KinisotError",
    "KinisotInputError",
    "KinisotParseError",
    "KinisotWarning",
    "compute_isotope_effect",
]
