"""Bigeleisen-Mayer thermodynamics: RPFRs and kinetic/equilibrium isotope effects.

All partition-function ingredients are accumulated in the log domain for
numerical stability; the exponentials are only taken when the final ratios
are formed.
"""

import math
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from .Hess_to_Freq import is_linear, read_hess

# PHYSICAL CONSTANTS (SI apart from speed of light)
PLANCK_CONSTANT = 6.62606957e-34  # m2 kg / s
BOLTZMANN_CONSTANT = 1.3806488e-23  # m2 kg s-2 K-1
SPEED_OF_LIGHT = 2.99792458e10  # cm s-1
ENERGY_AU = 4.35974434e-18  # J
BOHR_RADIUS = 5.2917721092e-11  # m
ATOMIC_MASS_UNIT = 1.660538921e-27  # kg

# Hartree/(Bohr^2 amu) -> (angular frequency in cm-1)^2
FORCE_TO_WAVENUMBER2 = (
    ENERGY_AU
    / (BOHR_RADIUS**2 * ATOMIC_MASS_UNIT)
    / ((SPEED_OF_LIGHT * 2 * np.pi) ** 2)
)


def log_product_factor(frequency_wn):
    """Log of the product of (scaled) vibrational frequencies, for the
    Teller-Redlich product factor. No temperature dependence."""
    f = np.asarray(frequency_wn) * SPEED_OF_LIGHT
    if f.size == 0:
        return 0.0
    return float(np.log(PLANCK_CONSTANT * f / BOLTZMANN_CONSTANT).sum())


def log_zpe_factor(frequency_wn, temperature):
    """Log of the vibrational ZPE term. The term enters the BM equation as
    exp(0.5 hv/kT); logarithmically that is just 0.5 hv/kT (which also
    avoids exp overflow at very low temperature)."""
    f = np.asarray(frequency_wn) * SPEED_OF_LIGHT
    if f.size == 0:
        return 0.0
    return float((0.5 * PLANCK_CONSTANT * f / (BOLTZMANN_CONSTANT * temperature)).sum())


def log_excitation_factor(frequency_wn, temperature):
    """Log of the excitation-factor term of the RPFR. Temperature dependent."""
    f = np.asarray(frequency_wn) * SPEED_OF_LIGHT
    if f.size == 0:
        return 0.0
    u = PLANCK_CONSTANT * f / (BOLTZMANN_CONSTANT * temperature)
    return float(np.log(-np.expm1(-u)).sum())


def skodje_truhlar_kappa(im_frequency_wn, temperature, barrier_j):
    """Skodje-Truhlar tunneling transmission coefficient for a truncated
    parabolic barrier (R. T. Skodje and D. G. Truhlar, J. Phys. Chem. 1981,
    85, 624-628), assuming an exothermic reaction.

    im_frequency_wn is the magnitude of the transition mode in cm-1 (already
    scaled), barrier_j the barrier height in Joules per molecule. Unlike the
    Bell infinite parabola, this remains defined below the crossover
    temperature: the finite barrier height caps the tunneling depth.
    """
    alpha = 2.0 * math.pi / (PLANCK_CONSTANT * SPEED_OF_LIGHT * im_frequency_wn)
    beta = 1.0 / (BOLTZMANN_CONSTANT * temperature)
    if abs(1.0 - beta / alpha) < 1e-10:
        # continuous limit at the crossover temperature
        return beta * barrier_j
    if beta < alpha:
        # above the crossover temperature; -> Bell as barrier -> infinity
        u_half = beta * math.pi / alpha
        return u_half / math.sin(u_half) - (beta / (alpha - beta)) * math.exp(
            (beta - alpha) * barrier_j
        )
    # below the crossover temperature: ground-state-dominated tunneling
    return (beta / (beta - alpha)) * math.expm1((beta - alpha) * barrier_j)


@dataclass
class RPFR:
    """Log-domain RPFR ingredients for one species (accumulated over its files)."""

    PF: float = 0.0
    ZPE: float = 0.0
    EXC: float = 0.0
    im_frequency_wn: Optional[float] = None
    # vibrational frequencies of the most recently processed file
    frequency_wn: List[float] = field(default_factory=list)


def rpfr(files, isomer, temperature=298.15, freq_scale_factor=1.0, freq_cutoff=50.0):
    """Compute the RPFR ingredients for a species (one or more files) with the
    given isotopic substitutions (one label string per file)."""
    result = RPFR()

    for i, file in enumerate(files):
        # Mass-weighted Hessian in Hartree/(amu Bohr^2), then harmonic
        # frequencies in cm-1 (imaginary modes as negative numbers)
        mw_hessmat = read_hess(file, isomer[i])
        eigs = np.linalg.eigvalsh(mw_hessmat * FORCE_TO_WAVENUMBER2)
        freqs = np.copysign(np.sqrt(np.abs(eigs)), eigs) * freq_scale_factor

        # 5 or 6 small normal modes will be removed (depending on whether the
        # molecule is linear or non-linear)
        trans_rot_modes = 5 if is_linear(file) == "linear" else 6

        # More than one large imaginary frequency indicates a higher-order
        # saddle point: the single-imaginary-mode treatment below is invalid
        n_imag = int((freqs < -freq_cutoff).sum())
        if n_imag > 1:
            warnings.warn(
                "{} imaginary frequencies larger than the cutoff found in {}; "
                "Kinisot treats only the first as the reaction coordinate".format(
                    n_imag, file
                )
            )

        # Keep a single imaginary frequency, if larger than the cutoff
        if abs(freqs[0]) > freq_cutoff:
            result.im_frequency_wn = -1.0 * float(freqs[0])
            trans_rot_modes += 1

        result.frequency_wn = [float(f) for f in freqs[trans_rot_modes:]]
        result.PF += log_product_factor(result.frequency_wn)
        result.ZPE += log_zpe_factor(result.frequency_wn, temperature)
        result.EXC += log_excitation_factor(result.frequency_wn, temperature)

    return result


@dataclass(frozen=True)
class IsotopeEffect:
    """Structured result of a KIE (or EQE) calculation for one isotopologue."""

    zpe: float  # ZPE factor ratio
    exc: float  # excitation factor ratio
    trpf: float  # Teller-Redlich product factor ratio
    v_ratio: float  # ratio of imaginary frequencies (1.0 for an EQE)
    kie: float  # uncorrected isotope effect
    bell_correction: float  # Bell infinite-parabola tunneling correction
    kie_bell: float  # kie * bell_correction
    wigner_correction: float  # Wigner tunneling correction
    kie_wigner: float  # kie * wigner_correction
    rpfrs: Tuple[
        RPFR, RPFR, RPFR, RPFR
    ]  # rct-light, rct-heavy, ts/prd-light, ts/prd-heavy
    temperature: float
    freq_scale_factor: float
    is_eqe: bool
    # Skodje-Truhlar correction (needs a barrier height; None when unavailable)
    st_correction: Optional[float] = None
    kie_st: Optional[float] = None
    barrier_light_j: Optional[float] = None

    def to_dict(self):
        return {
            "kie": self.kie,
            "kie_wigner": self.kie_wigner,
            "kie_bell": self.kie_bell,
            "zpe": self.zpe,
            "exc": self.exc,
            "trpf": self.trpf,
            "v_ratio": self.v_ratio,
            "wigner_correction": self.wigner_correction,
            "bell_correction": self.bell_correction,
            "kie_skodje_truhlar": self.kie_st,
            "skodje_truhlar_correction": self.st_correction,
            "temperature": self.temperature,
            "freq_scale_factor": self.freq_scale_factor,
            "is_eqe": self.is_eqe,
        }


def _assemble(rpfrs, is_kie, temperature, freq_scale_factor, electronic_barrier=None):
    """Build an IsotopeEffect from the four RPFRs
    (rct-light, rct-heavy, ts/prd-light, ts/prd-heavy).

    electronic_barrier (Hartree, optional) enables the Skodje-Truhlar
    correction; each isotopologue uses its own ZPE-corrected barrier,
    assembled from the electronic barrier and the RPFR ZPE sums."""
    if is_kie:
        if rpfrs[2].im_frequency_wn is None or rpfrs[3].im_frequency_wn is None:
            raise ValueError(
                "Kinisot requires a transition structure with an imaginary frequency!"
            )
        v_ratio = rpfrs[2].im_frequency_wn / rpfrs[3].im_frequency_wn
    else:
        v_ratio = 1.0

    # Application of the Bigeleisen-Mayer equation
    zpe = math.exp(rpfrs[0].ZPE - rpfrs[1].ZPE - rpfrs[2].ZPE + rpfrs[3].ZPE)
    exc = math.exp(rpfrs[0].EXC - rpfrs[1].EXC - rpfrs[2].EXC + rpfrs[3].EXC)
    trpf = math.exp(rpfrs[2].PF - rpfrs[3].PF - rpfrs[0].PF + rpfrs[1].PF)

    kie = v_ratio * zpe * exc * trpf

    # Tunneling corrections from the light/heavy TS imaginary frequencies:
    # Bell infinite parabola, kappa = (u/2)/sin(u/2), and Wigner,
    # kappa = 1 + u^2/24, with u = h*nu/kT
    if is_kie:
        tofreq = SPEED_OF_LIGHT * PLANCK_CONSTANT / BOLTZMANN_CONSTANT / temperature
        u_light = tofreq * rpfrs[2].im_frequency_wn
        u_heavy = tofreq * rpfrs[3].im_frequency_wn
        # The parabolic-barrier formula diverges at u/2 = pi and is undefined
        # beyond it, i.e. below the crossover temperature T_c = h*nu*c/(2*pi*k).
        # Reporting sin-ratio values there would be silently meaningless.
        if 0.5 * u_light >= math.pi:
            crossover = (
                SPEED_OF_LIGHT
                * PLANCK_CONSTANT
                * rpfrs[2].im_frequency_wn
                / (2.0 * math.pi * BOLTZMANN_CONSTANT)
            )
            warnings.warn(
                "temperature {:.2f} K is at or below the parabolic-barrier "
                "crossover temperature ({:.0f} K) for the {:.0f}i cm-1 "
                "transition mode: the Bell infinite-parabola correction is "
                "undefined (reported as nan); the Wigner and Skodje-Truhlar "
                "corrections stay finite but are unreliable in this "
                "deep-tunneling regime".format(
                    temperature, crossover, rpfrs[2].im_frequency_wn
                )
            )
            bell_correction = float("nan")
        else:
            bell_correction = (
                v_ratio * math.sin(0.5 * u_heavy) / math.sin(0.5 * u_light)
            )
        wigner_correction = (1.0 + u_light**2 / 24.0) / (1.0 + u_heavy**2 / 24.0)
    else:
        bell_correction = 1.0
        wigner_correction = 1.0

    # Skodje-Truhlar: truncated-parabola tunneling using per-isotopologue
    # ZPE-corrected barriers (RPFR.ZPE is the dimensionless 0.5*sum(hv/kT))
    st_correction = None
    kie_st = None
    barrier_light_j = None
    if is_kie and electronic_barrier is not None:
        kt = BOLTZMANN_CONSTANT * temperature
        elec_j = electronic_barrier * ENERGY_AU
        barrier_light_j = elec_j + (rpfrs[2].ZPE - rpfrs[0].ZPE) * kt
        barrier_heavy_j = elec_j + (rpfrs[3].ZPE - rpfrs[1].ZPE) * kt
        if barrier_light_j <= 0 or barrier_heavy_j <= 0:
            warnings.warn(
                "non-positive ZPE-corrected barrier ({:.1f} kJ/mol); "
                "skipping the Skodje-Truhlar correction".format(
                    barrier_light_j * 6.02214076e20
                )
            )
        else:
            kappa_light = skodje_truhlar_kappa(
                rpfrs[2].im_frequency_wn, temperature, barrier_light_j
            )
            kappa_heavy = skodje_truhlar_kappa(
                rpfrs[3].im_frequency_wn, temperature, barrier_heavy_j
            )
            st_correction = kappa_light / kappa_heavy
            kie_st = kie * st_correction

    return IsotopeEffect(
        zpe=zpe,
        exc=exc,
        trpf=trpf,
        v_ratio=v_ratio,
        kie=kie,
        bell_correction=bell_correction,
        kie_bell=kie * bell_correction,
        wigner_correction=wigner_correction,
        kie_wigner=kie * wigner_correction,
        rpfrs=tuple(rpfrs),
        temperature=temperature,
        freq_scale_factor=freq_scale_factor,
        is_eqe=not is_kie,
        st_correction=st_correction,
        kie_st=kie_st,
        barrier_light_j=barrier_light_j,
    )


def compute_isotope_effects(
    rct,
    ts=None,
    prd=None,
    labels=None,
    temperature=298.15,
    freq_scale_factor=1.0,
    freq_cutoff=50.0,
    electronic_barrier=None,
):
    """Compute isotope effects for several isotopologues in one pass.

    rct (and ts, or prd for an EQE) are lists of file paths. labels is a list
    of isotopologues, each holding one isotope-specification string per file
    (reactant files first). The unsubstituted (light) RPFRs are computed once
    and shared, so N isotopologues cost N+1 rather than 2N diagonalizations
    per species. Returns a list of IsotopeEffect.
    """
    if ts is not None:
        other, is_kie = ts, True
    elif prd is not None:
        other, is_kie = prd, False
    else:
        raise ValueError(
            "compute_isotope_effect needs a TS (KIE) or a product (EQE) in addition to the reactant"
        )

    light_rct = rpfr(rct, ["0"] * len(rct), temperature, freq_scale_factor, freq_cutoff)
    light_other = rpfr(
        other, ["0"] * len(other), temperature, freq_scale_factor, freq_cutoff
    )

    effects = []
    for label in labels:
        heavy_rct = rpfr(
            rct, label[0 : len(rct)], temperature, freq_scale_factor, freq_cutoff
        )
        heavy_other = rpfr(
            other, label[len(rct) :], temperature, freq_scale_factor, freq_cutoff
        )
        effects.append(
            _assemble(
                (light_rct, heavy_rct, light_other, heavy_other),
                is_kie,
                temperature,
                freq_scale_factor,
                electronic_barrier,
            )
        )
    return effects


def compute_isotope_effect(
    rct,
    ts=None,
    prd=None,
    label=None,
    temperature=298.15,
    freq_scale_factor=1.0,
    freq_cutoff=50.0,
    electronic_barrier=None,
):
    """Compute the isotope effect for one isotopologue.

    rct (and ts, or prd for an EQE) are lists of file paths; label holds one
    isotope-specification string per file, reactant files first. Returns an
    IsotopeEffect.
    """
    return compute_isotope_effects(
        rct,
        ts,
        prd,
        [label],
        temperature,
        freq_scale_factor,
        freq_cutoff,
        electronic_barrier,
    )[0]
