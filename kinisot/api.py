"""The Python API: compute_kie() and the result dataclasses it returns.

    from kinisot import compute_kie
    r = compute_kie(rct="claisen_gs.out", ts="claisen_ts.out", iso="5", temperature=393, scale=0.961)
    r.kie, r.kie_tunnel, r.zpe, r.exc, r.trpf, r.imag_ratio, r.to_dict()

Inputs may be file paths or :class:`~kinisot.hessian.HessianInput` objects,
so the same function serves every backend. The equations are in
docs/theory.md.
"""

import json
import os
import warnings
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from . import __version__
from .backends import load_hessian
from .exceptions import KinisotInputError, KinisotWarning
from .hessian import mass_weight
from .isotopes import substitute
from .projection import project_external_modes
from .scaling import ScalingChoice, choose_scaling_factor
from .thermo import (
    HARTREE_TO_KCAL_PER_MOL,
    TUNNELING_MODELS,
    harmonic_frequencies,
    log_excitation_factor,
    log_product_factor,
    log_zpe_factor,
    tunneling_correction,
)

__all__ = ["compute_kie", "IsotopeEffect", "SideResult", "IsotopologueResult", "SpeciesResult"]


def short_name(source):
    """File name without directory and extension, for tables."""
    return os.path.splitext(os.path.basename(str(source)))[0]


@dataclass(frozen=True)
class SpeciesResult:
    """One file (molecule) of one isotopologue, after diagonalization."""

    source: str
    label: str
    substitutions: Tuple  # Substitution records
    masses: Tuple[float, ...]
    frequencies: Tuple[float, ...]  # kept modes, scaled, cm-1, ascending
    discarded: Tuple[float, ...]  # external modes removed, scaled, cm-1
    imaginary: Optional[float]  # magnitude of the reaction-coordinate mode (scaled), or None
    linear: bool
    log_pf: float
    log_zpe: float
    log_exc: float
    projected: bool = False  # external modes projected out before diagonalization

    @property
    def name(self):
        return short_name(self.source)

    def to_dict(self):
        return {
            "source": self.source,
            "label": self.label,
            "substitutions": [s.to_dict() for s in self.substitutions],
            "masses": list(self.masses),
            "frequencies": list(self.frequencies),
            "discarded_modes": list(self.discarded),
            "imaginary_frequency": self.imaginary,
            "linear": self.linear,
            "projected": self.projected,
            "log_product_factor": self.log_pf,
            "log_zpe_factor": self.log_zpe,
            "log_excitation_factor": self.log_exc,
        }


@dataclass(frozen=True)
class IsotopologueResult:
    """One isotopologue of one side of the reaction (one or more files)."""

    species: Tuple[SpeciesResult, ...]

    @property
    def log_pf(self):
        return sum(s.log_pf for s in self.species)

    @property
    def log_zpe(self):
        return sum(s.log_zpe for s in self.species)

    @property
    def log_exc(self):
        return sum(s.log_exc for s in self.species)

    @property
    def imaginary(self):
        """The reaction-coordinate frequency when exactly one file has one, else None."""
        values = [s.imaginary for s in self.species if s.imaginary is not None]
        return values[0] if len(values) == 1 else None

    @property
    def substitutions(self):
        return tuple(sub for s in self.species for sub in s.substitutions)

    @property
    def labels(self):
        return tuple(s.label for s in self.species)

    @property
    def name(self):
        return " + ".join(s.name for s in self.species)

    @property
    def frequencies(self):
        return tuple(f for s in self.species for f in s.frequencies)

    def to_dict(self):
        return {
            "species": [s.to_dict() for s in self.species],
            "log_product_factor": self.log_pf,
            "log_zpe_factor": self.log_zpe,
            "log_excitation_factor": self.log_exc,
            "imaginary_frequency": self.imaginary,
        }


@dataclass(frozen=True)
class SideResult:
    """Light and heavy isotopologues of one side, with their Bigeleisen-Mayer factors."""

    light: IsotopologueResult
    heavy: IsotopologueResult

    @property
    def zpe_factor(self):
        """prod exp[(u_L - u_H) / 2]."""
        return float(np.exp(self.light.log_zpe - self.heavy.log_zpe))

    @property
    def exc_factor(self):
        """prod (1 - exp(-u_L)) / (1 - exp(-u_H))."""
        return float(np.exp(self.light.log_exc - self.heavy.log_exc))

    @property
    def trpf_factor(self):
        """prod nu_H / nu_L (Teller-Redlich product)."""
        return float(np.exp(self.heavy.log_pf - self.light.log_pf))

    @property
    def rpfr(self):
        """The reduced isotopic partition function ratio (s/s')f of this side."""
        return self.zpe_factor * self.exc_factor * self.trpf_factor

    @property
    def name(self):
        return self.light.name

    def to_dict(self):
        return {
            "name": self.name,
            "light": self.light.to_dict(),
            "heavy": self.heavy.to_dict(),
            "zpe_factor": self.zpe_factor,
            "excitation_factor": self.exc_factor,
            "product_factor": self.trpf_factor,
            "rpfr": self.rpfr,
        }


@dataclass(frozen=True)
class IsotopeEffect:
    """Result of compute_kie(): the isotope effect and everything that went into it."""

    kind: str  # "KIE" or "EQE"
    temperature: float
    scaling: ScalingChoice
    imag_cutoff: float
    tunneling: str
    reactant: SideResult
    other: SideResult  # the transition structure (KIE) or the product (EQE)
    imag_ratio: float  # nu_L / nu_H of the reaction coordinate (1 for an EQE)
    zpe: float
    exc: float
    trpf: float
    kie: float  # semiclassical, without tunnelling
    tunnel_corr: float
    kie_tunnel: float
    warnings: Tuple[str, ...] = field(default_factory=tuple)
    version: str = __version__
    project: bool = False
    barrier: Optional[float] = None  # kcal/mol, used by the Skodje-Truhlar correction
    reference: Optional["IsotopeEffect"] = None  # a second isotopologue the KIE is divided by

    @property
    def kie_relative(self):
        """KIE divided by the reference isotopologue's KIE (None without a reference)."""
        return self.kie / self.reference.kie if self.reference is not None else None

    @property
    def kie_tunnel_relative(self):
        return self.kie_tunnel / self.reference.kie_tunnel if self.reference is not None else None

    @property
    def scale_factor(self):
        return self.scaling.factor

    @property
    def transition_structure(self):
        return self.other if self.kind == "KIE" else None

    @property
    def product(self):
        return self.other if self.kind == "EQE" else None

    @property
    def species(self):
        """The four isotopologues in the order reactant light, heavy, other light, heavy."""
        return (self.reactant.light, self.reactant.heavy, self.other.light, self.other.heavy)

    def to_dict(self):
        """A JSON-serializable dictionary with every number Kinisot computed."""
        return {
            "kinisot_version": self.version,
            "kind": self.kind,
            "temperature": self.temperature,
            "scale_factor": self.scaling.factor,
            "scale_source": self.scaling.source,
            "level_of_theory": self.scaling.level,
            "imag_cutoff": self.imag_cutoff,
            "tunneling": self.tunneling,
            "barrier_kcal": self.barrier,
            "project": self.project,
            "imag_ratio": self.imag_ratio,
            "zpe": self.zpe,
            "exc": self.exc,
            "trpf": self.trpf,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "reactant": self.reactant.to_dict(),
            "transition_structure" if self.kind == "KIE" else "product": self.other.to_dict(),
            "warnings": list(self.warnings),
            "reference": self.reference.to_dict() if self.reference is not None else None,
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }

    def to_json(self, indent=2):
        return json.dumps(self.to_dict(), indent=indent)

    def summary_row(self):
        """One flat row (for CSV): the numbers of the results line plus provenance."""
        return {
            "kind": self.kind,
            "reactant": ";".join(s.source for s in self.reactant.light.species),
            "transition_structure_or_product": ";".join(s.source for s in self.other.light.species),
            "labels": ";".join(self.reactant.heavy.labels + self.other.heavy.labels),
            "temperature": self.temperature,
            "scale_factor": self.scaling.factor,
            "project": self.project,
            "tunneling": self.tunneling,
            "barrier_kcal": self.barrier,
            "imag_light": self.other.light.imaginary,
            "imag_heavy": self.other.heavy.imaginary,
            "imag_ratio": self.imag_ratio,
            "zpe": self.zpe,
            "exc": self.exc,
            "trpf": self.trpf,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "reference_labels": ";".join(self.reference.reactant.heavy.labels + self.reference.other.heavy.labels)
            if self.reference is not None
            else "",
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }


def evaluate_species(data, label, temperature, scale, imag_cutoff, warnings_out, project=False):
    """Diagonalize one HessianInput for one isotope label."""
    notes = []
    masses, applied = substitute(data, label, notes)
    for message in notes:
        if message not in warnings_out:
            warnings_out.append(message)
            warnings.warn(message, KinisotWarning, stacklevel=4)

    mw_hessian = mass_weight(data.hessian, masses)
    if project:
        if data.positions is None:
            raise KinisotInputError(
                "%s: projecting the external modes needs the geometry, which this input does not carry" % data.source
            )
        mw_hessian, n_external = project_external_modes(mw_hessian, data.positions, masses)
        freqs = harmonic_frequencies(mw_hessian) * scale
        # the external modes are now ~0: remove them by magnitude, keep the rest in order
        order = np.argsort(np.abs(freqs))
        discarded = np.sort(freqs[order[:n_external]])
        remaining = np.sort(freqs[order[n_external:]])
    else:
        freqs = harmonic_frequencies(mw_hessian) * scale
        n_external = 5 if data.linear else 6
        discarded = None
        remaining = None

    if project:
        imaginary_modes = remaining[remaining < -imag_cutoff]
        imaginary = -float(remaining[0]) if len(imaginary_modes) else None
        kept = remaining[1:] if imaginary is not None else remaining
        extra_discarded = np.array([])
    else:
        # 5 or 6 external modes are removed (linear / non-linear molecule), plus one
        # reaction-coordinate mode when it is imaginary beyond the cutoff
        imaginary_modes = freqs[freqs < -imag_cutoff]
        imaginary = -float(freqs[0]) if len(imaginary_modes) else None
        n_drop = n_external + (1 if imaginary is not None else 0)
        discarded = freqs[(1 if imaginary is not None else 0) : n_drop]
        kept = freqs[n_drop:]
        extra_discarded = np.array([])
    if len(imaginary_modes) > 1:
        message = (
            "%s has %d imaginary frequencies beyond the %.1f cm-1 cutoff (%s); only the largest is "
            "treated as the reaction coordinate and the others %s. Check the structure."
            % (
                data.source,
                len(imaginary_modes),
                imag_cutoff,
                ", ".join("%.1fi" % -f for f in imaginary_modes),
                "stay in the vibrational product as negative frequencies, which is an error below"
                if project
                else "are discarded with the external modes, which leaves one low-frequency external mode in the "
                "vibrational product",
            )
        )
        warnings_out.append(message)
        warnings.warn(message, KinisotWarning, stacklevel=4)
    if np.any(kept <= 0):
        raise KinisotInputError(
            "%s: %d non-positive frequencies remain after removing %d external modes (%s); the structure "
            "is not a stationary point Kinisot can use"
            % (data.source, int(np.sum(kept <= 0)), n_external, ", ".join("%.1f" % f for f in kept[kept <= 0]))
        )
    del extra_discarded

    if not applied and data.program_frequencies is not None:
        # Self-check for the unsubstituted species: Kinisot must reproduce the program's frequencies
        # (which are projected, hence the 1 cm-1 tolerance). Catches unit and mass-convention mistakes.
        theirs = np.sort(np.asarray(data.program_frequencies, dtype=float))
        if project:
            mine = np.sort(remaining) / scale
        else:
            # Drop the 5/6 modes closest to zero regardless of the cutoff, so the check is cutoff independent.
            mine = np.sort(freqs[np.argsort(np.abs(freqs))[n_external:]]) / scale
        largest = float(np.abs(mine - theirs).max()) if len(mine) == len(theirs) else float("nan")
        if len(mine) != len(theirs) or largest > 1.0:
            message = (
                "%s: the frequencies Kinisot computes from the Hessian differ from the %d the program printed "
                "(largest difference %.2f cm-1); check that the file is a completed frequency job and that the "
                "Hessian and the masses belong to the same geometry" % (data.source, len(theirs), largest)
            )
            warnings_out.append(message)
            warnings.warn(message, KinisotWarning, stacklevel=4)
    return SpeciesResult(
        source=str(data.source),
        label=str(label),
        substitutions=tuple(applied),
        masses=tuple(masses),
        frequencies=tuple(float(f) for f in kept),
        discarded=tuple(float(f) for f in discarded),
        imaginary=imaginary,
        linear=data.linear,
        log_pf=log_product_factor(kept),
        log_zpe=log_zpe_factor(kept, temperature),
        log_exc=log_excitation_factor(kept, temperature),
        projected=project,
    )


def evaluate_isotopologue(inputs, labels, temperature, scale, imag_cutoff, warnings_out, project=False):
    """Evaluate one isotopologue of one side (one label per input)."""
    species = tuple(
        evaluate_species(data, label, temperature, scale, imag_cutoff, warnings_out, project)
        for data, label in zip(inputs, labels)
    )
    with_imaginary = [s for s in species if s.imaginary is not None]
    if len(with_imaginary) > 1:
        raise KinisotInputError(
            "more than one file on the same side of the reaction has an imaginary frequency: %s"
            % ", ".join("%s (%.1fi cm-1)" % (s.source, s.imaginary) for s in with_imaginary)
        )
    return IsotopologueResult(species)


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (str, bytes, os.PathLike)) or not hasattr(value, "__iter__"):
        return [value]
    return list(value)


def normalize_labels(iso, n_reactants, n_other):
    """One label per file. A single label is used for both sides when there are two files."""
    if iso is None:
        raise KinisotInputError("iso is required: the atom number(s) to substitute (see --help)")
    labels = [str(x) for x in ([iso] if isinstance(iso, (str, int)) else list(iso))]
    n_files = n_reactants + n_other
    if len(labels) == 1 and n_files == 2:
        labels = labels * 2
    if len(labels) != n_files:
        raise KinisotInputError(
            "%d files were given but %d isotope labels; give one label per file (0 for no substitution)"
            % (n_files, len(labels))
        )
    return labels


def _describe(substitutions):
    return ", ".join(str(s) for s in substitutions) if substitutions else "none"


def _check_substitution_balance(reactant, other, side_name):
    """Both sides must carry the same isotopic substitutions (same elements)."""
    left = sorted(s.symbol for s in reactant.heavy.substitutions)
    right = sorted(s.symbol for s in other.heavy.substitutions)
    if not left and not right:
        raise KinisotInputError("no isotopic substitution requested: every --iso label is '0'")
    if left != right:
        raise KinisotInputError(
            "the isotopic substitutions differ between the reactant side [%s] and the %s side [%s]. "
            "An isotope effect compares the same isotopologue on both sides; check the atom numbering "
            "in each file." % (_describe(reactant.heavy.substitutions), side_name, _describe(other.heavy.substitutions))
        )


def compute_kie(
    rct,
    ts=None,
    prd=None,
    iso=None,
    temperature=298.15,
    scale=1.0,
    imag_cutoff=50.0,
    tunneling="bell",
    scale_type="zpe",
    project=False,
    barrier=None,
    reference=None,
):
    """Compute a kinetic (``ts``) or equilibrium (``prd``) isotope effect.

    Parameters
    ----------
    rct, ts, prd : path, HessianInput, or a list of them
        Reactant(s) and either the transition structure(s) or the product(s).
    iso : str or list of str
        Atom numbers to substitute, one label per file in the order of the
        reactants followed by the transition structure/product ('0' for a
        file without substitution, '7,8' for two atoms). A single label is
        applied to both files when there are exactly two.
    temperature : float, K.
    scale : float or None
        Vibrational scaling factor; None looks it up in the Truhlar database
        for the detected level of theory (1.0 if not found).
    imag_cutoff : float, cm-1
        A mode below -imag_cutoff is the reaction coordinate.
    tunneling : 'bell' (default), 'wigner', 'skodje' (Skodje-Truhlar) or 'none'.
    scale_type : which Truhlar factor to use when ``scale`` is None: 'zpe'
        (default), 'harm' or 'fund'.
    project : project translations and rotations out of the Hessian before
        diagonalizing (needs geometries) instead of discarding the 5/6
        lowest modes.
    barrier : barrier height in kcal/mol for the Skodje-Truhlar correction;
        by default the electronic energies read from the files
        (E(TS) - sum E(reactants)).
    reference : isotope label(s) of a second isotopologue; the result's
        ``reference`` holds its IsotopeEffect and ``kie_relative`` /
        ``kie_tunnel_relative`` the ratios.

    Returns
    -------
    IsotopeEffect

    Raises
    ------
    KinisotInputError, KinisotParseError
    Non-fatal problems are emitted as KinisotWarning when found and are also listed in ``result.warnings``.
    """
    if (ts is None) == (prd is None):
        raise KinisotInputError(
            "give either transition structure files (KIE) or product files (EQE), not both and not neither"
        )
    if temperature <= 0:
        raise KinisotInputError("temperature must be positive (got %s K)" % temperature)
    if scale is not None and scale <= 0:
        raise KinisotInputError("the vibrational scaling factor must be positive (got %s)" % scale)
    if tunneling not in TUNNELING_MODELS:
        raise KinisotInputError(
            "unknown tunnelling model %r (choose from %s)" % (tunneling, ", ".join(TUNNELING_MODELS))
        )

    kind = "KIE" if ts is not None else "EQE"
    side_name = "transition structure" if kind == "KIE" else "product"
    reactants = [load_hessian(x) for x in _as_list(rct)]
    others = [load_hessian(x) for x in _as_list(ts if kind == "KIE" else prd)]
    if not reactants or not others:
        raise KinisotInputError("at least one reactant and one %s file are required" % side_name)
    labels = normalize_labels(iso, len(reactants), len(others))
    scaling = choose_scaling_factor(reactants + others, scale, scale_type)

    collected = []
    n = len(reactants)
    reactant = SideResult(
        evaluate_isotopologue(reactants, ["0"] * n, temperature, scaling.factor, imag_cutoff, collected, project),
        evaluate_isotopologue(reactants, labels[:n], temperature, scaling.factor, imag_cutoff, collected, project),
    )
    other = SideResult(
        evaluate_isotopologue(
            others, ["0"] * len(others), temperature, scaling.factor, imag_cutoff, collected, project
        ),
        evaluate_isotopologue(others, labels[n:], temperature, scaling.factor, imag_cutoff, collected, project),
    )

    # Reactants (and products) must be minima
    minima = [("reactant", reactant.light)] + ([("product", other.light)] if kind == "EQE" else [])
    for role, isotopologue in minima:
        for s in isotopologue.species:
            if s.imaginary is not None:
                raise KinisotInputError(
                    "%s has an imaginary frequency (%.1fi cm-1) but was given as a %s. Reactants and products "
                    "must be minima; re-optimize the structure or, if this is a spurious low mode, raise "
                    "--imag-cutoff above %.1f" % (s.source, s.imaginary, role, s.imaginary)
                )

    if kind == "KIE":
        if other.light.imaginary is None or other.heavy.imaginary is None:
            raise KinisotInputError(
                "Kinisot requires a transition structure with an imaginary frequency beyond the %.1f cm-1 "
                "cutoff (--imag-cutoff)! Files given as --ts: %s" % (imag_cutoff, ", ".join(d.source for d in others))
            )
        imag_ratio = other.light.imaginary / other.heavy.imaginary
    else:
        imag_ratio = 1.0

    _check_substitution_balance(reactant, other, side_name)

    # Application of the Bigeleisen-Mayer equation
    zpe = reactant.zpe_factor / other.zpe_factor
    exc = reactant.exc_factor / other.exc_factor
    trpf = reactant.trpf_factor / other.trpf_factor
    kie = imag_ratio * zpe * exc * trpf
    if kind == "KIE":
        if tunneling == "skodje" and barrier is None:
            energies = [d.energy for d in reactants + others]
            if all(e is not None for e in energies):
                barrier = (sum(d.energy for d in others) - sum(d.energy for d in reactants)) * HARTREE_TO_KCAL_PER_MOL
        tunnel_corr = tunneling_correction(
            tunneling, other.light.imaginary, other.heavy.imaginary, temperature, barrier
        )
    else:
        tunnel_corr = 1.0

    reference_result = None
    if reference is not None:
        reference_result = compute_kie(
            reactants, others if kind == "KIE" else None, others if kind == "EQE" else None, iso=reference,
            temperature=temperature, scale=scaling.factor, imag_cutoff=imag_cutoff, tunneling=tunneling,
            scale_type=scale_type, project=project, barrier=barrier,
        )  # fmt: skip

    return IsotopeEffect(
        kind=kind,
        temperature=float(temperature),
        scaling=scaling,
        imag_cutoff=float(imag_cutoff),
        tunneling=tunneling if kind == "KIE" else "none",
        reactant=reactant,
        other=other,
        imag_ratio=float(imag_ratio),
        zpe=float(zpe),
        exc=float(exc),
        trpf=float(trpf),
        kie=float(kie),
        tunnel_corr=float(tunnel_corr),
        kie_tunnel=float(kie * tunnel_corr),
        warnings=tuple(collected),
        project=project,
        barrier=barrier if tunneling == "skodje" else None,
        reference=reference_result,
    )
