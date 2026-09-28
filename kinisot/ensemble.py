"""Conformer ensembles: isotope effects from several conformers of each species.

    from kinisot import Conformers, compute_kie
    r = compute_kie(rct="gs.out", ts=Conformers(["ts1.out", "ts2.out"]), iso="5", temperature=298.15)
    r.kie_tunnel, r.conformers, r.kie_lowest, r.n_effective

With transition-state theory and conformers in fast equilibrium (Curtin-Hammett),
isotopologue X reacts with

    k^X = (k_B T / h) sum_j kappa_j^X Q'_j^X exp(-E'_j / kT) / sum_i Q_i^X exp(-E_i / kT)

and, since every conformer of a species has the same atoms, the Teller-Redlich
mass factor cancels conformer by conformer:

    KIE = sum_i x_i rho_i / sum_j y_j rho'_j

- rho_i is the reduced isotopic partition-function ratio (s/s')f of reactant conformer i.
- rho'_j is that of transition-structure conformer j divided by (nu'_L/nu'_H)(kappa_L/kappa_H),
  so that rho_i / rho'_j is the ordinary pairwise KIE.
- x_i is the Boltzmann population of reactant conformer i for the light isotopologue.
- y_j is conformer j's share of the light isotopologue's rate, proportional to kappa_L,j exp(-G_j / RT).

Weights from the light isotopologue alone are exact, and only the free energies
within each ensemble matter. Several species on one side multiply. For an EQE the
product conformers take the place of the transition structures, with population
weights. docs/theory.md, section 7, has the derivation.
"""

import math
import warnings
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from . import __version__
from .api import IsotopeEffect, _as_list, compute_kie, evaluate_species, normalize_labels, short_name
from .backends import load_hessian
from .exceptions import KinisotInputError, KinisotWarning
from .isotopes import element_symbol
from .scaling import ScalingChoice, choose_scaling_factor
from .thermo import (
    ATOMIC_MASS_UNIT,
    BOHR_RADIUS,
    BOHR_TO_ANGSTROM,
    BOLTZMANN_CONSTANT,
    HARTREE_TO_KCAL_PER_MOL,
    PLANCK_CONSTANT,
    TUNNELING_MODELS,
    tunneling_kappa,
)

__all__ = ["Conformers", "ConformerResult", "EnsembleIsotopeEffect", "compute_ensemble", "equivalent_positions"]

WEIGHTS = ("qrrho", "rrho", "user", "lowest", "equal")
GAS_CONSTANT_KCAL = 1.987204259e-3  # kcal / (mol K)
QRRHO_CUTOFF = 100.0  # cm-1, Grimme's entropy interpolation, as in GoodVibes
ENERGY_UNITS = {"kcal/mol": 1.0, "kj/mol": 1.0 / 4.184, "hartree": HARTREE_TO_KCAL_PER_MOL}
# single-bond covalent radii (Angstrom; Cordero et al., Dalton Trans. 2008, 2832), for the bond-graph check
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Br": 1.20, "I": 1.39, "Li": 1.28, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "K": 2.03, "Ca": 1.76,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Sn": 1.39,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36,
}  # fmt: skip


@dataclass(frozen=True)
class Conformers:
    """Several conformers of one species, all with the same atoms in the same order.

    ``files`` are paths or HessianInput objects. ``free_energies`` (optional, one
    per conformer, in ``energy_unit``: 'kcal/mol', 'kJ/mol' or 'hartree'; only
    differences matter) replace the computed weights. ``degeneracy`` (optional)
    multiplies each conformer's weight, for example 2 for a conformer whose
    mirror image is not in the list.
    """

    files: tuple
    free_energies: Optional[tuple] = None
    degeneracy: Optional[tuple] = None
    energy_unit: str = "kcal/mol"

    def __post_init__(self):
        files = tuple(_as_list(self.files))
        if not files:
            raise KinisotInputError("a conformer ensemble needs at least one file")
        object.__setattr__(self, "files", files)
        unit = str(self.energy_unit).lower()
        if unit not in ENERGY_UNITS:
            raise KinisotInputError("unknown energy unit %r (choose kcal/mol, kJ/mol or hartree)" % self.energy_unit)
        object.__setattr__(self, "energy_unit", unit)
        for name in ("free_energies", "degeneracy"):
            values = getattr(self, name)
            if values is None:
                continue
            values = tuple(float(v) for v in _as_list(values))
            if len(values) != len(files):
                raise KinisotInputError("%d conformers but %d %s" % (len(files), len(values), name.replace("_", " ")))
            if name == "degeneracy" and min(values) <= 0:
                raise KinisotInputError("degeneracies must be positive")
            object.__setattr__(self, name, values)

    def __len__(self):
        return len(self.files)


@dataclass(frozen=True)
class ConformerResult:
    """One conformer of one species: its weight, its isotope ratio and its pairwise isotope effect."""

    role: str  # 'reactant', 'transition structure' or 'product'
    species: int  # index of the species on its side of the reaction
    source: str
    label: str
    free_energy: Optional[float]  # kcal/mol above the lowest conformer of the species (None for 'equal')
    weight_source: str  # 'qrrho', 'rrho', 'user' or 'equal'
    degeneracy: float
    population: float  # share of the light isotopologue (for a TS: of its rate, with tunnelling)
    population_semiclassical: float  # the same without the tunnelling factor
    rho: float  # (s/s')f of the conformer
    imaginary_light: Optional[float]
    imaginary_heavy: Optional[float]
    kappa_light: float
    kappa_heavy: float
    kie: float  # this conformer against the ensemble on the other side, semiclassical
    kie_tunnel: float
    light: object = field(repr=False, compare=False, default=None)  # SpeciesResult
    heavy: object = field(repr=False, compare=False, default=None)

    @property
    def name(self):
        return short_name(self.source)

    @property
    def rho_semiclassical(self):
        """rho' without tunnelling: rho / (nu_L / nu_H) for a transition structure, rho otherwise."""
        return self.rho / (self.imaginary_light / self.imaginary_heavy) if self.imaginary_light else self.rho

    @property
    def rho_tunnel(self):
        return self.rho_semiclassical / (self.kappa_light / self.kappa_heavy)

    def to_dict(self):
        return {
            "role": self.role,
            "species": self.species,
            "source": self.source,
            "label": self.label,
            "free_energy_kcal": self.free_energy,
            "weight_source": self.weight_source,
            "degeneracy": self.degeneracy,
            "population": self.population,
            "population_semiclassical": self.population_semiclassical,
            "rho": self.rho,
            "imaginary_light": self.imaginary_light,
            "imaginary_heavy": self.imaginary_heavy,
            "kappa_light": self.kappa_light,
            "kappa_heavy": self.kappa_heavy,
            "kie": self.kie,
            "kie_tunnel": self.kie_tunnel,
        }


@dataclass(frozen=True)
class EnsembleIsotopeEffect:
    """Result of compute_kie() when a species is given as several conformers.

    It carries the attributes the command line, the CSV writer and the benchmark
    runner read from an IsotopeEffect (kie, kie_tunnel, kie_relative, ...), plus
    the conformer table and diagnostics. There is no ensemble breakdown into
    ZPE, EXC and TRPF factors: those exist per conformer only.
    """

    kind: str
    temperature: float
    scaling: ScalingChoice
    imag_cutoff: float
    tunneling: str
    weights: str
    conformers: Tuple[ConformerResult, ...]
    kie: float
    kie_tunnel: float
    kie_lowest: Optional[float]  # from the lowest conformer of every species (None if they cannot be ranked)
    kie_tunnel_lowest: Optional[float]
    rho_reactant: float  # prod over reactant species of the ensemble-averaged rho
    rho_other: float  # the same for the other side, rho' without tunnelling
    rho_other_tunnel: float  # and with it
    n_effective: float  # 1 / sum y_j^2 of the ensemble on the transition-structure (or product) side
    kie_tunnel_range: Tuple[float, float]  # when each free energy moves by +/- weight_uncertainty
    weight_uncertainty: float
    warnings: Tuple[str, ...] = field(default_factory=tuple)
    version: str = __version__
    project: bool = False
    barrier: Optional[float] = None
    reference: Optional["EnsembleIsotopeEffect"] = None

    @property
    def tunnel_corr(self):
        return self.kie_tunnel / self.kie

    @property
    def kie_relative(self):
        return self.kie / self.reference.kie if self.reference is not None else None

    @property
    def kie_tunnel_relative(self):
        return self.kie_tunnel / self.reference.kie_tunnel if self.reference is not None else None

    @property
    def scale_factor(self):
        return self.scaling.factor

    def rows(self, role):
        return tuple(c for c in self.conformers if c.role == role)

    @property
    def labels(self):
        """The isotope label of each species, reactants first."""
        seen = {}
        for c in self.conformers:
            seen.setdefault((c.role != "reactant", c.species), c.label)
        return tuple(seen[k] for k in sorted(seen))

    def to_dict(self):
        return {
            "kinisot_version": self.version,
            "kind": self.kind,
            "ensemble": True,
            "temperature": self.temperature,
            "scale_factor": self.scaling.factor,
            "scale_source": self.scaling.source,
            "level_of_theory": self.scaling.level,
            "imag_cutoff": self.imag_cutoff,
            "tunneling": self.tunneling,
            "barrier_kcal": self.barrier,
            "project": self.project,
            "weights": self.weights,
            "weight_uncertainty_kcal": self.weight_uncertainty,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "kie_lowest": self.kie_lowest,
            "kie_tunnel_lowest": self.kie_tunnel_lowest,
            "rho_reactant": self.rho_reactant,
            "rho_other": self.rho_other,
            "rho_other_tunnel": self.rho_other_tunnel,
            "n_effective": self.n_effective,
            "kie_tunnel_range": list(self.kie_tunnel_range),
            "conformers": [c.to_dict() for c in self.conformers],
            "warnings": list(self.warnings),
            "reference": self.reference.to_dict() if self.reference is not None else None,
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }

    def to_json(self, indent=2):
        import json

        return json.dumps(self.to_dict(), indent=indent)

    def summary_row(self):
        """One flat row (for CSV)."""
        other_role = "transition structure" if self.kind == "KIE" else "product"
        return {
            "kind": self.kind,
            "reactant": ";".join(c.source for c in self.rows("reactant")),
            "transition_structure_or_product": ";".join(c.source for c in self.rows(other_role)),
            "labels": ";".join(self.labels),
            "temperature": self.temperature,
            "scale_factor": self.scaling.factor,
            "project": self.project,
            "tunneling": self.tunneling,
            "weights": self.weights,
            "conformers": len(self.conformers),
            "n_effective": self.n_effective,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "kie_lowest": self.kie_lowest,
            "kie_tunnel_lowest": self.kie_tunnel_lowest,
            "kie_tunnel_min": self.kie_tunnel_range[0],
            "kie_tunnel_max": self.kie_tunnel_range[1],
            "reference_labels": ";".join(self.reference.labels) if self.reference is not None else "",
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }


def is_ensemble(value):
    """True when ``value`` (an rct/ts/prd argument) gives some species as several conformers."""
    return any(isinstance(item, (Conformers, list, tuple)) for item in _as_list(value))


def as_species(value):
    """Each species of an rct/ts/prd argument as Conformers (a nested list is one species' conformers)."""
    species = []
    for item in _as_list(value):
        if isinstance(item, Conformers):
            species.append(item)
        elif isinstance(item, (list, tuple)):
            species.append(Conformers(tuple(item)))
        else:
            species.append(Conformers((item,)))
    return species


def _sides(result):
    """(rho of the reactant side, rho' of the other side without and with tunnelling) of a result."""
    if isinstance(result, EnsembleIsotopeEffect):
        return result.rho_reactant, result.rho_other, result.rho_other_tunnel
    other = result.other.rpfr / result.imag_ratio
    return result.reactant.rpfr, other, other / result.tunnel_corr


def equivalent_positions(results):
    """Isotope effect of positions made equivalent by fast motion, exactly.

    ``results`` holds one result of compute_kie per placement of the label over
    the equivalent positions (the three hydrogens of a rotating methyl group,
    the two oxygens of a nitro group), all from the same files. The placements
    are conformers with equal weights, so the isotope effect is the mean of rho
    on one side over the mean on the other: for a difference on one side only,
    the harmonic mean of the separate isotope effects. Ensemble results
    average the same way. Returns (kie, kie_tunnel).
    """
    sides = np.array([_sides(r) for r in results], dtype=float)
    if not len(sides):
        raise KinisotInputError("equivalent_positions needs at least one result")
    reactant, semiclassical, tunnel = sides.mean(axis=0)
    return float(reactant / semiclassical), float(reactant / tunnel)


# ------------------------------------------------------------------ validation


def _distances(positions):
    xyz = np.asarray(positions, dtype=float) * BOHR_TO_ANGSTROM
    return np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)


def _check_conformers(inputs, role, warnings_out):
    """Same atoms in the same order (error); same bonds and no duplicates (warnings)."""
    first = inputs[0]
    for data in inputs[1:]:
        if len(data.atomic_numbers) != len(first.atomic_numbers):
            raise KinisotInputError(
                "%s conformers must have the same atoms in the same order: %s has %d atoms but %s has %d"
                % (role, first.source, len(first.atomic_numbers), data.source, len(data.atomic_numbers))
            )
        for i, (a, b) in enumerate(zip(first.atomic_numbers, data.atomic_numbers)):
            if a != b:
                raise KinisotInputError(
                    "%s conformers must have the same atoms in the same order, so that a label means the same atom "
                    "in each: atom %d is %s in %s but %s in %s"
                    % (role, i + 1, element_symbol(a), first.source, element_symbol(b), data.source)
                )
    geometries = [d for d in inputs if d.positions is not None]
    if len(geometries) < 2:
        return
    radii = np.array([COVALENT_RADII.get(element_symbol(z), np.nan) for z in first.atomic_numbers])
    reach = radii[:, None] + radii[None, :]
    distances = [_distances(d.positions) for d in geometries]
    bonded = [d < 1.1 * reach for d in distances]
    apart = [d > 1.4 * reach for d in distances]
    for data, bond_d, apart_d in zip(geometries[1:], bonded[1:], apart[1:]):
        # a bond in one conformer that is plainly broken in the other: a renumbered atom, not a conformer
        moved = np.argwhere(np.triu((bonded[0] & apart_d) | (bond_d & apart[0]), 1))
        if len(moved):
            i, j = moved[0]
            message = (
                "%s conformers %s and %s differ in bonding (atoms %d %s and %d %s, and %d more pair(s)): check that "
                "both use the same atom numbering"
                % (role, geometries[0].source, data.source, i + 1, element_symbol(first.atomic_numbers[i]), j + 1,
                   element_symbol(first.atomic_numbers[j]), len(moved) - 1)
            )  # fmt: skip
            warnings_out.append(message)
            warnings.warn(message, KinisotWarning, stacklevel=5)
    for a in range(len(geometries)):
        for b in range(a + 1, len(geometries)):
            same_shape = np.abs(distances[a] - distances[b]).max() < 0.01
            ea, eb = geometries[a].energy, geometries[b].energy
            same_energy = ea is None or eb is None or abs(ea - eb) < 1e-5
            if same_shape and same_energy:
                message = (
                    "%s conformers %s and %s are the same structure (or mirror images), which doubles its weight; "
                    "remove one or give the other a degeneracy of 2"
                    % (role, geometries[a].source, geometries[b].source)
                )
                warnings_out.append(message)
                warnings.warn(message, KinisotWarning, stacklevel=5)
    levels = {d.level_of_theory for d in inputs if d.level_of_theory}
    if len(levels) > 1:
        message = "%s conformers come from different levels of theory (%s)" % (role, ", ".join(sorted(levels)))
        warnings_out.append(message)
        warnings.warn(message, KinisotWarning, stacklevel=5)


# ------------------------------------------------------------------ free energies


def _rotational_temperatures(positions, masses, linear):
    """Principal rotational temperatures (K) from a geometry in Bohr and masses in amu."""
    xyz = np.asarray(positions, dtype=float) * BOHR_RADIUS
    m = np.asarray(masses, dtype=float) * ATOMIC_MASS_UNIT
    xyz = xyz - (m[:, None] * xyz).sum(axis=0) / m.sum()
    tensor = np.einsum("i,ij,ik->jk", m, xyz, xyz)
    moments = np.linalg.eigvalsh(np.trace(tensor) * np.eye(3) - tensor)
    moments = moments[moments > 1e-50]
    temperatures = PLANCK_CONSTANT**2 / (8 * math.pi**2 * moments * BOLTZMANN_CONSTANT)
    return [float(max(temperatures))] if linear else [float(t) for t in temperatures]


def free_energy(data, light, temperature, scheme):
    """Relative free energy (kcal/mol) of one conformer's light isotopologue for 'qrrho' or 'rrho' weights.

    E + ZPE + thermal vibrational energy - T (S_vib + S_rot), from the scaled
    frequencies Kinisot kept (the reaction-coordinate mode excluded). 'qrrho'
    interpolates the vibrational entropy towards a free rotor below 100 cm-1
    (Grimme, as in GoodVibes); 'rrho' keeps it harmonic. Translation, the
    electronic entropy and the rotational energy are the same for every
    conformer of a species and are left out; so is the rotational symmetry
    number.
    """
    from goodvibes import thermo as gv

    if data.energy is None:
        raise KinisotInputError(
            "%s: no electronic energy was read, which '%s' weights need; give the free energies instead "
            "(Conformers(..., free_energies=...) or --energies) or use --weights equal" % (data.source, scheme)
        )
    frequencies = [f for f in light.frequencies if f > 0]
    energy = gv.calc_vibrational_energy(temperature, frequencies)  # J/mol, zero-point and thermal
    rrho = np.array(gv.calc_rrho_entropy(temperature, frequencies))
    if scheme == "qrrho":
        damp = np.array(gv.calc_damp(frequencies, QRRHO_CUTOFF))
        rotor = np.array(gv.calc_freerot_entropy(temperature, frequencies))
        vibrational_entropy = float(np.sum(damp * rrho + (1 - damp) * rotor))
    else:
        vibrational_entropy = float(np.sum(rrho))
    rotational_entropy = 0.0
    if data.positions is not None and len(light.masses) > 1:
        rotational = _rotational_temperatures(data.positions, light.masses, data.linear)
        rotational_entropy = gv.calc_rotational_entropy(temperature, rotational, linear=data.linear)
    joules = energy - temperature * (vibrational_entropy + rotational_entropy)
    return data.energy * HARTREE_TO_KCAL_PER_MOL + joules / 4184.0


# ------------------------------------------------------------------ evaluation


def _rho(light, heavy):
    """(s/s')f of one conformer: ZPE x EXC x TRPF factors, light over heavy."""
    return float(np.exp(light.log_zpe - heavy.log_zpe + light.log_exc - heavy.log_exc + heavy.log_pf - light.log_pf))


def _populations(energies, degeneracy, kappas, temperature, scheme):
    """Normalized weights g kappa exp(-G/RT) (or all on the lowest, or g alone for 'equal')."""
    g = np.asarray(degeneracy, dtype=float)
    if scheme == "equal":
        weights = g.copy()
    elif scheme == "lowest":
        weights = np.zeros(len(g))
        weights[int(np.argmin(energies))] = 1.0
    else:
        e = np.asarray(energies, dtype=float)
        weights = g * np.asarray(kappas, dtype=float) * np.exp(-(e - e.min()) / (GAS_CONSTANT_KCAL * temperature))
    return weights / weights.sum()


def _side_values(sides, energies, temperature, scheme, tunnel):
    """(prod over reactant species of <rho>, prod over the other side's species of <rho'>), with given energies."""
    values = []
    for side in ("reactant", "other"):
        value = 1.0
        for sp in sides[side]:
            kappa = [c["kappa_light"] if tunnel else 1.0 for c in sp]
            x = _populations([energies[id(c)] for c in sp], [c["degeneracy"] for c in sp], kappa, temperature, scheme)
            value *= float(np.dot(x, [c["rho_tunnel"] if tunnel else c["rho_semiclassical"] for c in sp]))
        values.append(value)
    return tuple(values)


def _ensemble_value(sides, energies, temperature, scheme, tunnel):
    """The isotope effect: prod over reactant species <rho> / prod over the other side's species <rho'>."""
    reactant, other = _side_values(sides, energies, temperature, scheme, tunnel)
    return reactant / other


def compute_ensemble(
    rct, ts=None, prd=None, iso=None, temperature=298.15, scale=1.0, imag_cutoff=50.0, tunneling="bell",
    scale_type="zpe", project=None, barrier=None, reference=None, calculator=None, delta=0.01, weights=None,
    weight_uncertainty=0.5,
):  # fmt: skip
    """compute_kie() for inputs with several conformers of a species; see the module docstring.

    ``weights``: 'qrrho' (default), 'rrho', 'user', 'lowest' or 'equal'. Free
    energies given with Conformers(...) are used whenever present, except with
    'equal'; 'user' requires them for every species with more than one conformer.
    """
    if (ts is None) == (prd is None):
        raise KinisotInputError(
            "give either transition structure files (KIE) or product files (EQE), not both and not neither"
        )
    if tunneling not in TUNNELING_MODELS:
        raise KinisotInputError(
            "unknown tunnelling model %r (choose from %s)" % (tunneling, ", ".join(TUNNELING_MODELS))
        )
    scheme = "qrrho" if weights is None else str(weights).lower()
    if scheme not in WEIGHTS:
        raise KinisotInputError("unknown weights %r (choose from %s)" % (weights, ", ".join(WEIGHTS)))
    if temperature <= 0:
        raise KinisotInputError("temperature must be positive (got %s K)" % temperature)
    kind = "KIE" if ts is not None else "EQE"
    other_role = "transition structure" if kind == "KIE" else "product"
    species = {"reactant": as_species(rct), "other": as_species(ts if kind == "KIE" else prd)}
    if not species["reactant"] or not species["other"]:
        raise KinisotInputError("at least one reactant and one %s are required" % other_role)

    # one conformer of every species and no free energies: the ordinary calculation, unchanged
    if all(len(c) == 1 for side in species.values() for c in side) and scheme not in ("user",):
        files = {side: [c.files[0] for c in species[side]] for side in species}
        return compute_kie(
            files["reactant"], files["other"] if kind == "KIE" else None, files["other"] if kind == "EQE" else None,
            iso=iso, temperature=temperature, scale=scale, imag_cutoff=imag_cutoff, tunneling=tunneling,
            scale_type=scale_type, project=project, barrier=barrier, reference=reference, calculator=calculator,
            delta=delta,
        )  # fmt: skip

    labels = normalize_labels(iso, len(species["reactant"]), len(species["other"]))
    collected = []
    inputs = {side: [[load_hessian(f, calculator, delta) for f in c.files] for c in species[side]] for side in species}
    everything = [d for side in inputs.values() for sp in side for d in sp]
    if project is None:
        project = any(str(d.program).lower().startswith("ase") for d in everything)
    project = bool(project)
    scaling = choose_scaling_factor(everything, scale, scale_type)
    for side in inputs:
        for k, sp in enumerate(inputs[side]):
            role = "reactant" if side == "reactant" else other_role
            _check_conformers(sp, "%s %d" % (role, k + 1) if len(inputs[side]) > 1 else role, collected)

    # evaluate every conformer: light once, heavy for the label
    n_rct = len(species["reactant"])
    evaluated = {"reactant": [], "other": []}
    for side in ("reactant", "other"):
        role = "reactant" if side == "reactant" else other_role
        for k, (conf, members) in enumerate(zip(species[side], inputs[side])):
            label = labels[k] if side == "reactant" else labels[n_rct + k]
            rows = []
            for i, data in enumerate(members):
                light = evaluate_species(data, "0", temperature, scaling.factor, imag_cutoff, collected, project)
                heavy = evaluate_species(data, label, temperature, scaling.factor, imag_cutoff, collected, project)
                rows.append(
                    {
                        "role": role, "species": k, "data": data, "label": label, "light": light, "heavy": heavy,
                        "rho": _rho(light, heavy),
                        "degeneracy": conf.degeneracy[i] if conf.degeneracy else 1.0,
                        "user_energy": conf.free_energies[i] * ENERGY_UNITS[conf.energy_unit]
                        if conf.free_energies else None,
                    }
                )  # fmt: skip
            evaluated[side].append(rows)

    # minima and transition structures
    for side, rows_by_species in evaluated.items():
        for rows in rows_by_species:
            with_mode = [r for r in rows if r["light"].imaginary is not None]
            if side == "reactant" or kind == "EQE":
                if with_mode:
                    r = with_mode[0]
                    raise KinisotInputError(
                        "%s has an imaginary frequency (%.1fi cm-1) but was given as a %s conformer; reactants and "
                        "products must be minima" % (r["data"].source, r["light"].imaginary, r["role"])
                    )
            elif with_mode and len(with_mode) != len(rows):
                r = next(r for r in rows if r["light"].imaginary is None)
                raise KinisotInputError(
                    "%s has no imaginary frequency beyond the %.1f cm-1 cutoff, but the other conformers of its "
                    "species do: every transition-structure conformer needs exactly one"
                    % (r["data"].source, imag_cutoff)
                )
    if kind == "KIE":
        with_ts = [rows for rows in evaluated["other"] if rows[0]["light"].imaginary is not None]
        if len(with_ts) != 1:
            raise KinisotInputError(
                "Kinisot requires exactly one transition-structure species with an imaginary frequency beyond the "
                "%.1f cm-1 cutoff; found %d" % (imag_cutoff, len(with_ts))
            )
    left = sorted(s.symbol for rows in evaluated["reactant"] for s in rows[0]["heavy"].substitutions)
    right = sorted(s.symbol for rows in evaluated["other"] for s in rows[0]["heavy"].substitutions)
    if not left and not right:
        raise KinisotInputError("no isotopic substitution requested: every --iso label is '0'")
    if left != right:
        raise KinisotInputError(
            "the isotopic substitutions differ between the reactant side and the %s side; an isotope effect compares "
            "the same isotopologue on both sides" % other_role
        )

    # tunnelling factors of each transition-structure conformer
    # (Skodje-Truhlar: each conformer's own barrier above the lowest reactant conformers, as compute_kie does)
    baseline = None
    if tunneling == "skodje" and barrier is None and kind == "KIE":
        beside = [rows for rows in evaluated["other"] if rows[0]["light"].imaginary is None]
        if all(r["data"].energy is not None for rows in evaluated["reactant"] + beside for r in rows):
            baseline = sum(min(r["data"].energy for r in rows) for rows in evaluated["reactant"]) - sum(
                min(r["data"].energy for r in rows) for rows in beside
            )
    for rows in evaluated["reactant"] + evaluated["other"]:
        for r in rows:
            r["kappa_light"] = r["kappa_heavy"] = 1.0
            if r["light"].imaginary is not None:
                nu_l, nu_h = r["light"].imaginary, r["heavy"].imaginary
                v = barrier
                if v is None and baseline is not None and r["data"].energy is not None:
                    v = (r["data"].energy - baseline) * HARTREE_TO_KCAL_PER_MOL
                r["kappa_light"] = tunneling_kappa(tunneling, nu_l, temperature, v)
                r["kappa_heavy"] = tunneling_kappa(tunneling, nu_h, temperature, v)
                r["rho_semiclassical"] = r["rho"] / (nu_l / nu_h)
            else:
                r["rho_semiclassical"] = r["rho"]
            r["rho_tunnel"] = r["rho_semiclassical"] / (r["kappa_light"] / r["kappa_heavy"])

    # free energies of the light isotopologues (relative within each species)
    energies, sources = {}, {}
    for rows in evaluated["reactant"] + evaluated["other"]:
        several = len(rows) > 1
        if scheme == "user" and several and rows[0]["user_energy"] is None:
            raise KinisotInputError(
                "--weights user needs a free energy for every conformer (Conformers(..., free_energies=...) or "
                "--energies); %s has none" % rows[0]["data"].source
            )
        for r in rows:
            if scheme == "equal" or not several:
                value, source = 0.0, "equal" if several else "single"
            elif r["user_energy"] is not None:
                value, source = r["user_energy"], "user"
            else:
                method = "rrho" if scheme == "rrho" else "qrrho"
                value, source = free_energy(r["data"], r["light"], temperature, method), method
            energies[id(r)], sources[id(r)] = value, source
        low = min(energies[id(r)] for r in rows)
        for r in rows:
            energies[id(r)] -= low

    sides = {side: evaluated[side] for side in ("reactant", "other")}
    rho_reactant, rho_other = _side_values(sides, energies, temperature, scheme, tunnel=False)
    rho_other_tunnel = _side_values(sides, energies, temperature, scheme, tunnel=True)[1]
    kie, kie_tunnel = rho_reactant / rho_other, rho_reactant / rho_other_tunnel
    # the lowest conformers: 'equal' weights set every free energy to zero, so rank by the given ones or by
    # computed qRRHO ones; without either (no electronic energies), there is no lowest conformer to report
    ranking = energies
    if scheme == "equal":
        ranking = {}
        for rows in evaluated["reactant"] + evaluated["other"]:
            if len(rows) < 2:
                ranking.update({id(r): 0.0 for r in rows})
            elif all(r["user_energy"] is not None for r in rows):
                ranking.update({id(r): r["user_energy"] for r in rows})
            elif all(r["data"].energy is not None for r in rows):
                ranking.update({id(r): free_energy(r["data"], r["light"], temperature, "qrrho") for r in rows})
            else:
                ranking = None
                break
    kie_lowest = kie_tunnel_lowest = None
    if ranking is not None:
        kie_lowest = _ensemble_value(sides, ranking, temperature, "lowest", tunnel=False)
        kie_tunnel_lowest = _ensemble_value(sides, ranking, temperature, "lowest", tunnel=True)

    # sensitivity: move each free energy by +/- weight_uncertainty in turn
    values = [kie_tunnel]
    if scheme not in ("equal", "lowest") and weight_uncertainty:
        for rows in evaluated["reactant"] + evaluated["other"]:
            if len(rows) < 2:
                continue
            for r in rows:
                for shift in (-weight_uncertainty, weight_uncertainty):
                    moved = dict(energies)
                    moved[id(r)] += shift
                    values.append(_ensemble_value(sides, moved, temperature, scheme, tunnel=True))

    # the conformer table: populations and each conformer's KIE against the other side's ensemble
    def side_average(side, tunnel, skip=None):
        total = 1.0
        for rows in evaluated[side]:
            if rows is skip:
                continue
            kappa = [r["kappa_light"] if tunnel else 1.0 for r in rows]
            x = _populations([energies[id(r)] for r in rows], [r["degeneracy"] for r in rows], kappa, temperature,
                             scheme)  # fmt: skip
            total *= float(np.dot(x, [r["rho_tunnel"] if tunnel else r["rho_semiclassical"] for r in rows]))
        return total

    conformers = []
    n_effective = 1.0
    for side in ("reactant", "other"):
        for rows in evaluated[side]:
            share = {
                t: _populations([energies[id(r)] for r in rows], [r["degeneracy"] for r in rows],
                                [r["kappa_light"] if t else 1.0 for r in rows], temperature, scheme)
                for t in (False, True)
            }  # fmt: skip
            if side == "other" and (kind == "EQE" or rows[0]["light"].imaginary is not None):
                n_effective = float(1.0 / np.sum(share[True] ** 2))
            for i, r in enumerate(rows):
                pair = {}
                for t in (False, True):
                    rho = r["rho_tunnel"] if t else r["rho_semiclassical"]
                    if side == "reactant":
                        pair[t] = rho * side_average("reactant", t, rows) / side_average("other", t)
                    else:
                        pair[t] = side_average("reactant", t) / (rho * side_average("other", t, rows))
                conformers.append(
                    ConformerResult(
                        role=r["role"], species=r["species"], source=str(r["data"].source), label=r["label"],
                        free_energy=None if sources[id(r)] in ("equal", "single") else energies[id(r)],
                        weight_source=sources[id(r)], degeneracy=r["degeneracy"],
                        population=float(share[True][i]), population_semiclassical=float(share[False][i]),
                        rho=r["rho"], imaginary_light=r["light"].imaginary, imaginary_heavy=r["heavy"].imaginary,
                        kappa_light=r["kappa_light"], kappa_heavy=r["kappa_heavy"], kie=pair[False],
                        kie_tunnel=pair[True], light=r["light"], heavy=r["heavy"],
                    )
                )  # fmt: skip

    reference_result = None
    if reference is not None:
        reference_result = compute_ensemble(
            species["reactant"], species["other"] if kind == "KIE" else None,
            species["other"] if kind == "EQE" else None, iso=reference, temperature=temperature,
            scale=scaling.factor, imag_cutoff=imag_cutoff, tunneling=tunneling, scale_type=scale_type,
            project=project, barrier=barrier, calculator=calculator, delta=delta, weights=scheme,
            weight_uncertainty=weight_uncertainty,
        )  # fmt: skip
        if isinstance(reference_result, IsotopeEffect):  # pragma: no cover - cannot happen with several conformers
            raise KinisotInputError("the reference must use the same conformers")

    return EnsembleIsotopeEffect(
        kind=kind,
        temperature=float(temperature),
        scaling=scaling,
        imag_cutoff=float(imag_cutoff),
        tunneling=tunneling if kind == "KIE" else "none",
        weights=scheme,
        conformers=tuple(conformers),
        kie=float(kie),
        kie_tunnel=float(kie_tunnel if kind == "KIE" else kie),
        kie_lowest=None if kie_lowest is None else float(kie_lowest),
        kie_tunnel_lowest=None if kie_lowest is None else float(kie_tunnel_lowest if kind == "KIE" else kie_lowest),
        rho_reactant=float(rho_reactant),
        rho_other=float(rho_other),
        rho_other_tunnel=float(rho_other_tunnel if kind == "KIE" else rho_other),
        n_effective=n_effective,
        kie_tunnel_range=(float(min(values)), float(max(values))),
        weight_uncertainty=float(weight_uncertainty or 0.0),
        warnings=tuple(collected),
        project=project,
        barrier=barrier if tunneling == "skodje" else None,
        reference=reference_result,
    )
