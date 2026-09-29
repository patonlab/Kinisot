"""Kinisot: kinetic and equilibrium isotope effects from computed Hessians.

Command line::

    kinisot --rct reactant.out --ts ts.out --iso 5 -t 393

Python::

    from kinisot import compute_kie
    result = compute_kie(rct="reactant.out", ts="ts.out", iso="5", temperature=393.0, scale=0.961)
    result.kie_tunnel, result.zpe, result.exc, result.trpf, result.to_dict()

    # a conformer ensemble: several files for one species
    result = compute_kie(rct="reactant.out", ts=[["ts_a.out", "ts_b.out"]], iso="5", temperature=393.0)
    result.kie_tunnel, result.conformers, result.kie_lowest, result.n_effective

    # transition structures in series, and parallel channels
    result = compute_kie(rct="gs.out", ts=Series(["ts_1.out", "ts_2.out"]), iso=["5", "5", "7"], temperature=393.0)
    result = channels([dict(rct="a.out", ts="ts_a.out", iso="5"), dict(rct="b.out", ts="ts_b.out", iso="3")],
                      shares=[1, 3.3])
"""

__version__ = "2.6.0"

from .api import IsotopeEffect, IsotopologueResult, SideResult, SpeciesResult, compute_kie
from .backends import load_hessian
from .backends.ase import build_calculator, hessian_for_geometry, hessian_from_calculator, save_hessian_json
from .backends.gaussian import parse_gaussian
from .ensemble import ConformerResult, Conformers, EnsembleIsotopeEffect, equivalent_positions
from .exceptions import KinisotError, KinisotInputError, KinisotParseError, KinisotWarning
from .hessian import HessianInput, mass_weight
from .isotopes import Substitution, isotope_mass, light_mass, substitute
from .Kinisot import compute_isotope_effect  # deprecated, removed in 3.0
from .pathways import ChannelIsotopeEffect, Series, SeriesIsotopeEffect, channel_kie, channels, series_kie
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
    "Conformers",
    "ConformerResult",
    "EnsembleIsotopeEffect",
    "equivalent_positions",
    "Series",
    "SeriesIsotopeEffect",
    "channels",
    "ChannelIsotopeEffect",
    "series_kie",
    "channel_kie",
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
