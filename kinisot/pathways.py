"""Transition structures in series and parallel channels (IMPLEMENTATION_PLAN.md, Phase 10).

    from kinisot import Series, channels, compute_kie
    r = compute_kie(rct=["aldehyde.out", "ylide.out"], ts=Series(["ts_4.out", "ts_6.out"]), iso=["5", "0", "12", "15"])
    r = channels([dict(rct="3.out", ts="ts_S.out", iso=["1", "40"]), dict(rct="3.out", ts="ts_R.out", iso=["3", "38"])],
                 shares=[3.3, 1])

**Series.** Steps R <=> I1 <=> I2 ... -> P, none of which alone commits the
substrate. At steady state the inverse rate constants add, 1/k = sum_n 1/k_n,
with k_n the rate constant from the starting material over transition
structure n. Hence

    KIE = sum_n w_n KIE_n,   w_n proportional to 1/k_n^L = exp(+G_n / RT) / kappa_L,n

- KIE_n is the ordinary KIE from the starting material to transition structure n.
- The weights are those of the light isotopologue, which is exact.
- The highest transition structure counts most, the opposite of conformers.
- For two steps this is KIE = (KIE_2 + C_f KIE_1) / (1 + C_f), with the
  commitment factor C_f = k_2 / k_-1 = exp[(G_1 - G_2) / RT].

**Channels.** Parallel routes, each with its own reactant, transition structure
or label (two enantiomers reacting through diastereomeric transition
structures to one product). At low conversion

    1 / KIE = sum_c y_c / KIE_c

where y_c is channel c's share of the light isotopologue's rate. For two
channels with selectivity s = y_2 / y_1 this is KIE = (1 + s) KIE_1 KIE_2 / (KIE_2 + s KIE_1).

A step or a channel may itself be a conformer ensemble. docs/theory.md,
section 7, has the derivations.
"""

import json
import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from . import __version__
from .api import IsotopeEffect, _as_list, compute_kie
from .backends import load_hessian
from .ensemble import (
    ENERGY_UNITS,
    GAS_CONSTANT_KCAL,
    WEIGHTS,
    Conformers,
    EnsembleIsotopeEffect,
    as_species,
    free_energy,
)
from .exceptions import KinisotInputError
from .thermo import tunneling_kappa

__all__ = [
    "Series",
    "SeriesIsotopeEffect",
    "ChannelIsotopeEffect",
    "channels",
    "compute_series",
    "series_kie",
    "channel_kie",
]


# ------------------------------------------------------------------ the formulas


def _series_weights(free_energies, temperature):
    """w_n proportional to exp(+G_n / RT), normalized."""
    x = np.asarray(free_energies, dtype=float) / (GAS_CONSTANT_KCAL * temperature)
    w = np.exp(x - x.max())
    return w / w.sum()


def _commitment_energies(commitment, temperature):
    """Step free energies (kcal/mol) equivalent to a commitment factor C_f = exp[(G_1 - G_2) / RT]."""
    return np.array([GAS_CONSTANT_KCAL * temperature * math.log(commitment), 0.0])


def series_kie(kies, free_energies=None, commitment=None, temperature=298.15, energy_unit="kcal/mol"):
    """Isotope effect of transition structures in series from the KIE of each step.

    ``kies`` are the KIEs from the starting material to each transition
    structure, in order. Give either their free energies (any common zero, in
    ``energy_unit``; KIE = sum_n w_n KIE_n with w_n proportional to
    exp(+G_n / RT)) or, for two steps, the commitment factor C_f = k_2 / k_-1
    (KIE = (KIE_2 + C_f KIE_1) / (1 + C_f)).
    """
    kies = np.asarray(_as_list(kies), dtype=float)
    return float(np.dot(_given_series_weights(len(kies), free_energies, commitment, temperature, energy_unit), kies))


def _given_series_weights(n, free_energies, commitment, temperature, energy_unit):
    if (free_energies is None) == (commitment is None):
        raise KinisotInputError("give the steps' free energies or, for two steps, a commitment factor (one of them)")
    if commitment is not None:
        if n != 2 or commitment <= 0:
            raise KinisotInputError("a commitment factor needs exactly two steps and must be positive")
        return _series_weights(_commitment_energies(commitment, temperature), temperature)
    energies = np.asarray(_as_list(free_energies), dtype=float) * _unit(energy_unit)
    if len(energies) != n:
        raise KinisotInputError("%d steps but %d free energies" % (n, len(energies)))
    return _series_weights(energies, temperature)


def channel_kie(kies, shares):
    """Isotope effect of parallel channels: 1 / KIE = sum_c y_c / KIE_c, with y_c the normalized ``shares``."""
    kies = np.asarray(_as_list(kies), dtype=float)
    y = _normalized(shares, len(kies), "shares")
    return float(1.0 / np.dot(y, 1.0 / kies))


def _normalized(values, n, name):
    values = np.asarray(_as_list(values), dtype=float)
    if len(values) != n:
        raise KinisotInputError("%d channels but %d %s" % (n, len(values), name))
    if (values < 0).any() or values.sum() <= 0:
        raise KinisotInputError("%s must be non-negative and not all zero" % name)
    return values / values.sum()


def _unit(energy_unit):
    unit = str(energy_unit).lower()
    if unit not in ENERGY_UNITS:
        raise KinisotInputError("unknown energy unit %r (choose kcal/mol, kJ/mol or hartree)" % energy_unit)
    return ENERGY_UNITS[unit]


# ------------------------------------------------------------------ free energies of species in results


def _formula(data):
    return tuple(sorted(data.atomic_numbers))


def _members(result, role, index, inputs):
    """(data, light SpeciesResult, degeneracy, dG within the species or None, kappa_L) per conformer."""
    if isinstance(result, EnsembleIsotopeEffect):
        rows = [c for c in result.conformers if c.role == role and c.species == index]
        return [(d, c.light, c.degeneracy, c.free_energy, c.kappa_light) for d, c in zip(inputs, rows)]
    side = result.reactant if role == "reactant" else result.other
    light = side.light.species[index]
    kappa = 1.0
    if light.imaginary is not None and result.kind == "KIE":
        kappa = tunneling_kappa(result.tunneling, light.imaginary, result.temperature, result.barrier)
    return [(inputs[0], light, 1.0, None, kappa)]


def _anchor(members):
    """The conformer the species' free energy refers to: its lowest one, or the first without energies."""
    energies = [m[3] for m in members]
    known = [i for i, e in enumerate(energies) if e is not None]
    return min(known, key=lambda i: energies[i]) if known else 0


def _log_sum(members, temperature, scheme, tunnel):
    """ln sum_j g_j kappa_j exp(-dG_j / RT) over the conformers of a species, as its weights count them."""
    anchor = _anchor(members)
    if scheme == "lowest":
        return math.log(members[anchor][4]) if tunnel else 0.0
    total = 0.0
    for _data, _light, g, dg, kappa in members:
        term = g * (kappa if tunnel else 1.0)
        if scheme != "equal" and dg is not None:
            term *= math.exp(-dg / (GAS_CONSTANT_KCAL * temperature))
        total += term
    return math.log(total)


def _computed_anchor(members, temperature, scheme, context):
    data, light = members[_anchor(members)][:2]
    try:
        return free_energy(data, light, temperature, "rrho" if scheme == "rrho" else "qrrho")
    except KinisotInputError as err:
        raise KinisotInputError("%s; or give %s" % (err, context)) from None


def _effective(anchor_energy, members, temperature, scheme):
    """(G_eff semiclassical, G_eff with tunnelling) = G_anchor - RT ln sum_j g_j kappa_j exp(-dG_j / RT)."""
    rt = GAS_CONSTANT_KCAL * temperature
    return tuple(anchor_energy - rt * _log_sum(members, temperature, scheme, t) for t in (False, True))


# ------------------------------------------------------------------ inputs


def _loaded(species, calculator, delta):
    """Conformers whose files are loaded once (HessianInput), so every step and channel reuses them."""
    return Conformers(
        tuple(load_hessian(f, calculator, delta) for f in species.files),
        free_energies=species.free_energies,
        degeneracy=species.degeneracy,
        energy_unit=species.energy_unit,
    )


def _argument(species):
    """What compute_kie takes for loaded species: plain inputs when every species is one conformer."""
    if all(len(c) == 1 and c.free_energies is None for c in species):
        return [c.files[0] for c in species]
    return list(species)


def _labels(iso, n):
    labels = [str(x) for x in ([iso] if isinstance(iso, (str, int)) else _as_list(iso))]
    if len(labels) == 1:
        labels = labels * n
    if len(labels) != n:
        raise KinisotInputError(
            "%d species were given but %d isotope labels; give one per reactant and one per step (or a single label "
            "when the numbering is the same in every file)" % (n, len(labels))
        )
    return labels


def _result_labels(result):
    if isinstance(result, EnsembleIsotopeEffect):
        return list(result.labels)
    if isinstance(result, IsotopeEffect):
        return list(result.reactant.heavy.labels + result.other.heavy.labels)
    return list(result.labels)


# ------------------------------------------------------------------ series


@dataclass(frozen=True)
class Series:
    """Transition structures in series, passed as ``ts``: ``compute_kie(rct, ts=Series([...]), iso=...)``.

    ``steps`` are the transition structures in the order the reaction meets
    them; each is a file, a HessianInput, a list of conformers or Conformers.
    Their weights come from (in order of precedence) ``commitment`` (two steps:
    C_f = k_2 / k_-1, for example from trajectories), ``free_energies`` (one per
    step, any common zero, in ``energy_unit``: the free energy of the step's
    lowest conformer), or free energies Kinisot computes (qRRHO, or RRHO with
    weights='rrho'), which needs the same atoms in every step.
    """

    steps: tuple
    free_energies: Optional[tuple] = None
    commitment: Optional[float] = None
    energy_unit: str = "kcal/mol"

    def __post_init__(self):
        steps = tuple(_as_list(self.steps))
        if not steps:
            raise KinisotInputError("a series needs at least one step")
        object.__setattr__(self, "steps", steps)
        _unit(self.energy_unit)
        if self.free_energies is not None and self.commitment is not None:
            raise KinisotInputError("give the steps' free energies or a commitment factor, not both")
        if self.free_energies is not None:
            energies = tuple(float(v) for v in _as_list(self.free_energies))
            if len(energies) != len(steps):
                raise KinisotInputError("%d steps but %d free energies" % (len(steps), len(energies)))
            object.__setattr__(self, "free_energies", energies)
        if self.commitment is not None:
            if len(steps) != 2 or not float(self.commitment) > 0:
                raise KinisotInputError("a commitment factor needs exactly two steps and must be positive")
            object.__setattr__(self, "commitment", float(self.commitment))

    def __len__(self):
        return len(self.steps)


class _Combined:
    """What SeriesIsotopeEffect and ChannelIsotopeEffect share."""

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
    def scaling(self):
        return self.parts[0].scaling

    @property
    def scale_factor(self):
        return self.scaling.factor

    @property
    def tunneling(self):
        return self.parts[0].tunneling

    @property
    def project(self):
        return self.parts[0].project

    @property
    def warnings(self):
        seen = []
        for part in self.parts:
            seen.extend(w for w in part.warnings if w not in seen)
        return tuple(seen)

    def to_json(self, indent=2):
        return json.dumps(self.to_dict(), indent=indent)

    def _common(self):
        return {
            "kinisot_version": self.version,
            "kind": self.kind,
            "temperature": self.temperature,
            "scale_factor": self.scaling.factor,
            "scale_source": self.scaling.source,
            "level_of_theory": self.scaling.level,
            "tunneling": self.tunneling,
            "project": self.project,
            "weights": self.weights,
            "weight_uncertainty_kcal": self.weight_uncertainty,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "kie_tunnel_range": list(self.kie_tunnel_range) if self.kie_tunnel_range is not None else None,
        }

    def _tail(self):
        return {
            "warnings": list(self.warnings),
            "reference": self.reference.to_dict() if self.reference is not None else None,
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }


@dataclass(frozen=True)
class SeriesIsotopeEffect(_Combined):
    """Result of compute_kie() with ``ts=Series(...)``: the combined KIE and the KIE of each step."""

    temperature: float
    steps: Tuple  # IsotopeEffect or EnsembleIsotopeEffect from the starting material to each transition structure
    free_energies: Tuple[float, ...]  # effective free energy of each step (kcal/mol, the lowest set to 0)
    weight_source: str  # 'commitment', 'user' or the computed scheme ('qrrho' or 'rrho')
    shares: Tuple[float, ...]  # w_n, with tunnelling
    shares_semiclassical: Tuple[float, ...]
    kie: float
    kie_tunnel: float
    kie_tunnel_range: Optional[Tuple[float, float]]  # each step's free energy moved by +/- weight_uncertainty
    weights: str
    weight_uncertainty: float
    reference: Optional["SeriesIsotopeEffect"] = None
    version: str = __version__
    kind: str = "KIE"

    @property
    def parts(self):
        return self.steps

    @property
    def commitment(self):
        """C_f = k_2 / k_-1 of a two-step series (with tunnelling), else None."""
        return self.shares[0] / self.shares[1] if len(self.steps) == 2 else None

    @property
    def labels(self):
        return [_result_labels(s) for s in self.steps]

    def kie_at(self, commitment, tunnel=True, relative=False):
        """The KIE of a two-step series for another commitment factor C_f."""
        k1, k2 = self._pair(tunnel, relative)
        return float((k2[0] + commitment * k1[0]) / (k2[1] + commitment * k1[1]))

    def commitment_for(self, value, tunnel=True, relative=False):
        """The commitment factor C_f of a two-step series that gives ``value`` (for example a measured KIE)."""
        (k1, r1), (k2, r2) = self._pair(tunnel, relative)
        commitment = (value * r2 - k2) / (k1 - value * r1)
        if not commitment >= 0:
            low, high = sorted((k1 / r1, k2 / r2))
            raise KinisotInputError("%.4f is outside what the two steps can give (%.4f to %.4f)" % (value, low, high))
        return float(commitment)

    def _pair(self, tunnel, relative):
        if len(self.steps) != 2:
            raise KinisotInputError("defined for two steps only")
        if relative and self.reference is None:
            raise KinisotInputError("relative values need a reference isotopologue")
        pairs = []
        for i, step in enumerate(self.steps):
            k = step.kie_tunnel if tunnel else step.kie
            r = 1.0
            if relative:
                ref = self.reference.steps[i]
                r = ref.kie_tunnel if tunnel else ref.kie
            pairs.append((k, r))
        return pairs

    def to_dict(self):
        data = self._common()
        data.update(
            {
                "series": True,
                "weight_source": self.weight_source,
                "free_energies_kcal": list(self.free_energies),
                "shares": list(self.shares),
                "shares_semiclassical": list(self.shares_semiclassical),
                "commitment": self.commitment,
                "steps": [s.to_dict() for s in self.steps],
            }
        )
        data.update(self._tail())
        return data

    def summary_row(self):
        return {
            "kind": "KIE (series)",
            "steps": ";".join("+".join(_sources(s, "other")) for s in self.steps),
            "reactant": ";".join(_sources(self.steps[0], "reactant")),
            "labels": " | ".join(";".join(labels) for labels in self.labels),
            "temperature": self.temperature,
            "scale_factor": self.scale_factor,
            "project": self.project,
            "tunneling": self.tunneling,
            "weight_source": self.weight_source,
            "shares": ";".join("%.6f" % y for y in self.shares),
            "commitment": self.commitment,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "kie_tunnel_min": self.kie_tunnel_range[0] if self.kie_tunnel_range else None,
            "kie_tunnel_max": self.kie_tunnel_range[1] if self.kie_tunnel_range else None,
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }


def _sources(result, side):
    if isinstance(result, EnsembleIsotopeEffect):
        role = "reactant" if side == "reactant" else ("transition structure" if result.kind == "KIE" else "product")
        return [c.source for c in result.rows(role)]
    if isinstance(result, IsotopeEffect):
        return [s.source for s in (result.reactant if side == "reactant" else result.other).light.species]
    return [p for part in result.parts for p in _sources(part, side)]


def compute_series(
    rct, series, iso=None, temperature=298.15, scale=1.0, imag_cutoff=50.0, tunneling="bell", scale_type="zpe",
    project=None, barrier=None, reference=None, calculator=None, delta=0.01, weights=None, weight_uncertainty=0.5,
):  # fmt: skip
    """compute_kie() with ``ts=Series(...)``; see the module docstring.

    ``iso`` gives one label per reactant species and then one per step (a
    single label when the numbering is the same everywhere); so does
    ``reference``.
    """
    scheme = "qrrho" if weights is None else str(weights).lower()
    if scheme not in WEIGHTS:
        raise KinisotInputError("unknown weights %r (choose from %s)" % (weights, ", ".join(WEIGHTS)))
    if temperature <= 0:
        raise KinisotInputError("temperature must be positive (got %s K)" % temperature)
    reactants = [_loaded(s, calculator, delta) for s in as_species(rct)]
    steps = []
    for step in series.steps:
        species = as_species([step] if not isinstance(step, (list, tuple)) else [list(step)])
        steps.append(_loaded(species[0], calculator, delta))
    n_rct = len(reactants)
    if not n_rct:
        raise KinisotInputError("at least one reactant is required")
    labels = _labels(iso, n_rct + len(steps))
    references = _labels(reference, n_rct + len(steps)) if reference is not None else None
    options = dict(
        temperature=temperature, scale=scale, imag_cutoff=imag_cutoff, tunneling=tunneling, scale_type=scale_type,
        project=project, barrier=barrier, calculator=calculator, delta=delta, weights=weights,
        weight_uncertainty=weight_uncertainty,
    )  # fmt: skip
    results = []
    for n, step in enumerate(steps):
        results.append(
            compute_kie(
                rct=_argument(reactants),
                ts=_argument([step]),
                iso=labels[:n_rct] + [labels[n_rct + n]],
                reference=references[:n_rct] + [references[n_rct + n]] if references else None,
                **options,
            )  # fmt: skip
        )
    if any(r.kind != "KIE" for r in results):  # pragma: no cover - compute_kie with ts always gives a KIE
        raise KinisotInputError("a series is made of transition structures")

    # effective free energy of each step: its lowest conformer's G, less RT ln of its conformer sum
    members = [_members(r, "transition structure", 0, s.files) for r, s in zip(results, steps)]
    if series.commitment is not None:
        source = "commitment"
        semiclassical = tunnel = _commitment_energies(series.commitment, temperature)
    else:
        if series.free_energies is not None:
            source = "user"
            anchors = np.array(series.free_energies) * _unit(series.energy_unit)
        elif scheme == "user":
            raise KinisotInputError(
                "--weights user needs the steps' free energies (Series(..., free_energies=...)) or a commitment factor"
            )
        else:
            source = "rrho" if scheme == "rrho" else "qrrho"
            formulas = {_formula(s.files[0]) for s in steps}
            if len(formulas) > 1:
                raise KinisotInputError(
                    "the steps of a series have different atoms, so their computed free energies cannot be compared; "
                    "give Series(..., free_energies=...) or a commitment factor"
                )
            anchors = np.array(
                [_computed_anchor(m, temperature, scheme, "Series(..., free_energies=...)") for m in members]
            )
        effective = np.array([_effective(a, m, temperature, scheme) for a, m in zip(anchors, members)])
        semiclassical, tunnel = effective[:, 0], effective[:, 1]
    share = {False: _series_weights(semiclassical, temperature), True: _series_weights(tunnel, temperature)}
    kie = float(np.dot(share[False], [r.kie for r in results]))
    kie_tunnel = float(np.dot(share[True], [r.kie_tunnel for r in results]))

    values = [kie_tunnel]
    if weight_uncertainty:
        for n in range(len(results) if len(results) > 1 else 0):
            for shift in (-weight_uncertainty, weight_uncertainty):
                moved = np.array(tunnel, dtype=float)
                moved[n] += shift
                values.append(float(np.dot(_series_weights(moved, temperature), [r.kie_tunnel for r in results])))

    reference_result = None
    if references is not None:
        reference_result = SeriesIsotopeEffect(
            temperature=float(temperature), steps=tuple(r.reference for r in results),
            free_energies=tuple(float(g) for g in np.asarray(tunnel) - np.min(tunnel)), weight_source=source,
            shares=tuple(float(y) for y in share[True]), shares_semiclassical=tuple(float(y) for y in share[False]),
            kie=float(np.dot(share[False], [r.reference.kie for r in results])),
            kie_tunnel=float(np.dot(share[True], [r.reference.kie_tunnel for r in results])),
            kie_tunnel_range=None, weights=scheme, weight_uncertainty=float(weight_uncertainty or 0.0),
        )  # fmt: skip
    return SeriesIsotopeEffect(
        temperature=float(temperature),
        steps=tuple(results),
        free_energies=tuple(float(g) for g in np.asarray(tunnel) - np.min(tunnel)),
        weight_source=source,
        shares=tuple(float(y) for y in share[True]),
        shares_semiclassical=tuple(float(y) for y in share[False]),
        kie=kie,
        kie_tunnel=kie_tunnel,
        kie_tunnel_range=(float(min(values)), float(max(values))),
        weights=scheme,
        weight_uncertainty=float(weight_uncertainty or 0.0),
        reference=reference_result,
    )


# ------------------------------------------------------------------ channels


@dataclass(frozen=True)
class ChannelIsotopeEffect(_Combined):
    """Result of channels(): the combined KIE and the KIE of each channel."""

    temperature: float
    names: Tuple[str, ...]
    channels: Tuple  # the result of each channel (IsotopeEffect, EnsembleIsotopeEffect or SeriesIsotopeEffect)
    share_source: str  # 'given', 'barriers' or the computed scheme ('qrrho' or 'rrho')
    shares: Tuple[float, ...]  # y_c, share of the light isotopologue's rate, with tunnelling
    shares_semiclassical: Tuple[float, ...]
    kie: float
    kie_tunnel: float
    kie_tunnel_range: Optional[Tuple[float, float]]  # each channel's barrier moved by +/- weight_uncertainty
    weights: str
    weight_uncertainty: float
    reference: Optional["ChannelIsotopeEffect"] = None
    version: str = __version__
    kind: str = "KIE"

    @property
    def parts(self):
        return self.channels

    @property
    def selectivity(self):
        """s = y_2 / y_1 of two channels (with tunnelling), else None."""
        return self.shares[1] / self.shares[0] if len(self.channels) == 2 else None

    @property
    def labels(self):
        return [_result_labels(c) for c in self.channels]

    def kie_at(self, shares, tunnel=True, relative=False):
        """The KIE for other channel shares (normalized here)."""
        kies = [c.kie_tunnel if tunnel else c.kie for c in self.channels]
        value = channel_kie(kies, shares)
        if relative:
            if self.reference is None:
                raise KinisotInputError("relative values need a reference isotopologue")
            value /= self.reference.kie_at(shares, tunnel)
        return value

    def selectivity_for(self, value, tunnel=True, relative=False):
        """The selectivity s = y_2 / y_1 of two channels that gives ``value`` (for example a measured KIE)."""
        if len(self.channels) != 2:
            raise KinisotInputError("defined for two channels only")
        k1, k2 = (c.kie_tunnel if tunnel else c.kie for c in self.channels)
        r1 = r2 = 1.0
        if relative:
            if self.reference is None:
                raise KinisotInputError("relative values need a reference isotopologue")
            r1, r2 = (c.kie_tunnel if tunnel else c.kie for c in self.reference.channels)
        s = (value / k1 - 1.0 / r1) / (1.0 / r2 - value / k2)
        if not s >= 0:
            low, high = sorted((k1 / r1, k2 / r2))
            raise KinisotInputError(
                "%.4f is outside what the two channels can give (%.4f to %.4f)" % (value, low, high)
            )
        return float(s)

    def to_dict(self):
        data = self._common()
        data.update(
            {
                "channels": True,
                "names": list(self.names),
                "share_source": self.share_source,
                "shares": list(self.shares),
                "shares_semiclassical": list(self.shares_semiclassical),
                "selectivity": self.selectivity,
                "results": [c.to_dict() for c in self.channels],
            }
        )
        data.update(self._tail())
        return data

    def summary_row(self):
        return {
            "kind": "KIE (channels)",
            "channels": ";".join(self.names),
            "labels": " | ".join(";".join(labels) for labels in self.labels),
            "temperature": self.temperature,
            "scale_factor": self.scale_factor,
            "project": self.project,
            "tunneling": self.tunneling,
            "share_source": self.share_source,
            "shares": ";".join("%.6f" % y for y in self.shares),
            "selectivity": self.selectivity,
            "kie": self.kie,
            "tunnel_corr": self.tunnel_corr,
            "kie_tunnel": self.kie_tunnel,
            "kie_tunnel_min": self.kie_tunnel_range[0] if self.kie_tunnel_range else None,
            "kie_tunnel_max": self.kie_tunnel_range[1] if self.kie_tunnel_range else None,
            "kie_relative": self.kie_relative,
            "kie_tunnel_relative": self.kie_tunnel_relative,
        }


CHANNEL_KEYS = ("rct", "ts", "iso", "reference", "name", "amount")


def channels(
    jobs, shares=None, barriers=None, amounts=None, energy_unit="kcal/mol", temperature=298.15, scale=1.0,
    imag_cutoff=50.0, tunneling="bell", scale_type="zpe", project=None, barrier=None, calculator=None, delta=0.01,
    weights=None, weight_uncertainty=0.5,
):  # fmt: skip
    """Isotope effect of parallel channels: 1 / KIE = sum_c y_c / KIE_c.

    ``jobs`` holds one dictionary per channel with ``rct``, ``ts`` (a file,
    conformers or a Series), ``iso`` and optionally ``reference``, ``name`` and
    ``amount``. The remaining keywords are those of compute_kie and apply to
    every channel.

    The shares y_c (of the light isotopologue's rate) come from, in order:

    - ``shares``: given directly, for example [1, 3.3] for a measured
      selectivity of 3.3 in favour of the second channel;
    - ``barriers``: each channel's free-energy barrier (kcal/mol, or
      ``energy_unit``) from its lowest reactant conformers to its lowest
      transition-structure conformer; y_c is proportional to
      amount_c kappa_c exp(-barrier_c / RT), with the conformer sums added;
    - otherwise free energies Kinisot computes (qRRHO, or RRHO with
      weights='rrho'), which needs the same molecules in every channel.

    ``amounts`` (or ``amount`` in a job) are the relative amounts of the
    channels' reactants, which do not interconvert (1 each by default; 1 and
    1 for the two enantiomers of a racemate).
    """
    jobs = list(jobs)
    if not jobs:
        raise KinisotInputError("channels needs at least one channel")
    if shares is not None and barriers is not None:
        raise KinisotInputError("give the channels' shares or their barriers, not both")
    scheme = "qrrho" if weights is None else str(weights).lower()
    if scheme not in WEIGHTS:
        raise KinisotInputError("unknown weights %r (choose from %s)" % (weights, ", ".join(WEIGHTS)))
    for job in jobs:
        unknown = set(job) - set(CHANNEL_KEYS)
        if unknown or "rct" not in job or "ts" not in job or "iso" not in job:
            raise KinisotInputError(
                "each channel needs rct, ts and iso (and may give %s); got %s"
                % (", ".join(k for k in CHANNEL_KEYS[3:]), ", ".join(sorted(job)))
            )
    if amounts is not None:
        amounts = _normalized(amounts, len(jobs), "amounts") * len(jobs)
    else:
        amounts = np.array([float(job.get("amount", 1.0)) for job in jobs])
        if (amounts <= 0).any():
            raise KinisotInputError("amounts must be positive")
    if barriers is not None:
        barriers = np.asarray(_as_list(barriers), dtype=float) * _unit(energy_unit)
        if len(barriers) != len(jobs):
            raise KinisotInputError("%d channels but %d barriers" % (len(jobs), len(barriers)))
    if (shares is None and barriers is None) and scheme == "user":
        raise KinisotInputError("--weights user needs the channels' shares or barriers")
    options = dict(
        temperature=temperature, scale=scale, imag_cutoff=imag_cutoff, tunneling=tunneling, scale_type=scale_type,
        project=project, barrier=barrier, calculator=calculator, delta=delta, weights=weights,
        weight_uncertainty=weight_uncertainty,
    )  # fmt: skip

    if any((job.get("reference") is None) != (jobs[0].get("reference") is None) for job in jobs):
        raise KinisotInputError("give a reference for every channel or for none")

    # every channel, with its inputs loaded once
    loaded, results = [], []
    for job in jobs:
        reactants = [_loaded(s, calculator, delta) for s in as_species(job["rct"])]
        ts = job["ts"]
        if isinstance(ts, (list, tuple)) and len(ts) == 1 and isinstance(ts[0], Series):
            ts = ts[0]
        other = None if isinstance(ts, Series) else [_loaded(s, calculator, delta) for s in as_species(ts)]
        loaded.append((reactants, other))
        results.append(
            compute_kie(
                rct=_argument(reactants),
                ts=ts if other is None else _argument(other),
                iso=job["iso"],
                reference=job.get("reference"),
                **options,
            )  # fmt: skip
        )
    names = tuple(str(job.get("name", "channel %d" % (c + 1))) for c, job in enumerate(jobs))

    # effective barriers: ln y_c = ln amount_c - barrier_c / RT, the barrier from the reactants' effective free
    # energy to the transition structures', each G_anchor - RT ln sum_j g_j kappa_j exp(-dG_j / RT)
    rt = GAS_CONSTANT_KCAL * temperature
    effective = {False: np.zeros(len(jobs)), True: np.zeros(len(jobs))}
    if shares is not None:
        source = "given"
    else:
        source = "barriers" if barriers is not None else ("rrho" if scheme == "rrho" else "qrrho")
        if barriers is None:
            composition = {
                (tuple(sorted(_formula(s.files[0]) for s in r)), tuple(sorted(_formula(s.files[0]) for s in o or [])))
                for r, o in loaded
            }
            if any(o is None for r, o in loaded) or len(composition) > 1:
                raise KinisotInputError(
                    "the channels do not have the same molecules (or one is a series), so their computed free "
                    "energies cannot be compared; give their shares or barriers"
                )
        for c, (result, (reactants, other)) in enumerate(zip(results, loaded)):
            if other is None:  # a series channel: its barrier is used as given
                effective[False][c] = effective[True][c] = barriers[c]
                continue
            for role, side, sign in (("reactant", reactants, -1.0), ("transition structure", other, 1.0)):
                for k, species in enumerate(side):
                    members = _members(result, role, k, species.files)
                    anchor = 0.0
                    if barriers is None:
                        anchor = _computed_anchor(members, temperature, scheme, "the channels' shares or barriers")
                    semiclassical, tunnel = _effective(anchor, members, temperature, scheme)
                    effective[False][c] += sign * semiclassical
                    effective[True][c] += sign * tunnel
            if barriers is not None:
                effective[False][c] += barriers[c]
                effective[True][c] += barriers[c]

    def shares_for(barrier):
        x = np.log(amounts) - np.asarray(barrier, dtype=float) / rt
        w = np.exp(x - x.max())
        return w / w.sum()

    if shares is not None:
        y = {t: _normalized(shares, len(jobs), "shares") for t in (False, True)}
    else:
        y = {t: shares_for(effective[t]) for t in (False, True)}
    kie = channel_kie([r.kie for r in results], y[False])
    kie_tunnel = channel_kie([r.kie_tunnel for r in results], y[True])
    values = [kie_tunnel]
    if weight_uncertainty and source != "given" and len(jobs) > 1:
        for c in range(len(jobs)):
            for shift in (-weight_uncertainty, weight_uncertainty):
                moved = np.array(effective[True], dtype=float)
                moved[c] += shift
                values.append(channel_kie([r.kie_tunnel for r in results], shares_for(moved)))

    reference_result = None
    if jobs[0].get("reference") is not None:
        reference_result = ChannelIsotopeEffect(
            temperature=float(temperature), names=names, channels=tuple(r.reference for r in results),
            share_source=source, shares=tuple(float(v) for v in y[True]),
            shares_semiclassical=tuple(float(v) for v in y[False]),
            kie=channel_kie([r.reference.kie for r in results], y[False]),
            kie_tunnel=channel_kie([r.reference.kie_tunnel for r in results], y[True]),
            kie_tunnel_range=None, weights=scheme, weight_uncertainty=float(weight_uncertainty or 0.0),
        )  # fmt: skip
    return ChannelIsotopeEffect(
        temperature=float(temperature),
        names=names,
        channels=tuple(results),
        share_source=source,
        shares=tuple(float(v) for v in y[True]),
        shares_semiclassical=tuple(float(v) for v in y[False]),
        kie=kie,
        kie_tunnel=kie_tunnel,
        kie_tunnel_range=(float(min(values)), float(max(values))),
        weights=scheme,
        weight_uncertainty=float(weight_uncertainty or 0.0),
        reference=reference_result,
    )
