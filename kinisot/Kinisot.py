"""Backwards-compatibility shim for the historical ``kinisot.Kinisot`` module.

The implementation now lives in :mod:`kinisot.thermo` (physics),
:mod:`kinisot.scaling` (scaling factors) and :mod:`kinisot.cli` (command
line). This module preserves the pre-2.1 function-and-tuple API and will be
removed in a future release.
"""

import warnings

from . import thermo
from .cli import Logger, main  # noqa: F401 (re-exported)
from .Hess_to_Freq import is_linear, level_of_theory, read_hess  # noqa: F401
from .scaling import find_scaling_factor  # noqa: F401
from .thermo import (  # noqa: F401
    ATOMIC_MASS_UNIT,
    BOHR_RADIUS,
    BOLTZMANN_CONSTANT,
    ENERGY_AU,
    PLANCK_CONSTANT,
    SPEED_OF_LIGHT,
)


def _deprecated(old, new):
    warnings.warn(
        "kinisot.Kinisot.{} is deprecated; use {} instead".format(old, new),
        DeprecationWarning,
        stacklevel=3,
    )


def calc_product_factor(frequency_wn):
    return thermo.log_product_factor(frequency_wn)


def calc_zpe_factor(frequency_wn, temperature):
    return thermo.log_zpe_factor(frequency_wn, temperature)


def calc_excitation_factor(frequency_wn, temperature):
    return thermo.log_excitation_factor(frequency_wn, temperature)


class calc_rpfr:
    """Deprecated class-style RPFR; use :func:`kinisot.thermo.rpfr`."""

    def __init__(
        self, files, isomer, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0
    ):
        _deprecated("calc_rpfr", "kinisot.thermo.rpfr")
        result = thermo.rpfr(files, isomer, temperature, freq_scale_factor, freq_cutoff)
        self.PF, self.ZPE, self.EXC = result.PF, result.ZPE, result.EXC
        self.frequency_wn = result.frequency_wn
        # the historical class only had the attribute when a mode was found
        if result.im_frequency_wn is not None:
            self.im_frequency_wn = result.im_frequency_wn


def compute_isotope_effect(
    rct,
    ts=None,
    prd=None,
    label=None,
    temperature=298.15,
    freq_scale_factor=1.0,
    freq_cutoff=50.0,
):
    """Deprecated tuple-returning API; use :func:`kinisot.compute_isotope_effect`,
    which returns a structured :class:`kinisot.thermo.IsotopeEffect`."""
    _deprecated("compute_isotope_effect", "kinisot.compute_isotope_effect")
    r = thermo.compute_isotope_effect(
        rct, ts, prd, label, temperature, freq_scale_factor, freq_cutoff
    )
    return (
        list(r.rpfrs),
        r.zpe,
        r.exc,
        r.trpf,
        r.kie,
        r.kie_bell,
        r.bell_correction,
        r.v_ratio,
    )


def get_frequency_scaling(files, log):
    """Deprecated; use :func:`kinisot.scaling.get_frequency_scaling`."""
    _deprecated("get_frequency_scaling", "kinisot.scaling.get_frequency_scaling")
    from .scaling import get_frequency_scaling as _get

    factor, level, ref, warning = _get(files)
    if warning is not None:
        log.Writeonlyfile("\nWARNING: " + warning)
    if ref is not None:
        log.Write(
            "\n  Found vibrational scaling factor {} for {} level of theory".format(
                factor, level
            )
        )
        log.Write("\n  REF: " + ref)
    else:
        log.Write(
            "\n  Unable to find vibrational scaling factor for {}; using value of 1.0".format(
                level
            )
        )
    return factor


if __name__ == "__main__":
    main()
