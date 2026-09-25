"""Choosing the vibrational scaling factor for a calculation.

The factors come from the Truhlar group database (version 5) as shipped by
GoodVibes (``goodvibes.vib_scale_factors``), whose ``canonicalize_level``
maps program-specific spellings (Gaussian ``PBE1PBE``, ORCA ``M06-2X``,
``6-31G*`` versus ``6-31G(d)``) onto the database keys.
"""

from dataclasses import dataclass
from typing import Optional, Tuple

from goodvibes.vib_scale_factors import FUNCTIONAL_ALIASES, canonicalize_level, scaling_data_dict, scaling_refs

from .exceptions import KinisotInputError

__all__ = ["find_scaling_factor", "choose_scaling_factor", "ScalingChoice", "SCALE_TYPES"]

# Spellings the GoodVibes 4.4 alias table does not cover yet (submitted upstream);
# harmless once GoodVibes carries them.
for _alias, _canonical in {
    "MN15-L": "MN15L",
    "MN12-L": "MN12L",
    "MN12-SX": "MN12SX",
    "M06-L(DKH2)": "M06L(DKH2)",
}.items():
    FUNCTIONAL_ALIASES.setdefault(_alias, _canonical)

# --scale-type name -> (factor field, reference field) of a GoodVibes ScalingData entry
SCALE_TYPES = {
    "zpe": ("zpe_fac", "zpe_ref"),
    "harm": ("harm_fac", "harm_ref"),
    "fund": ("fund_fac", "fund_ref"),
}


def _candidates(level):
    """Canonical keys to try: as given, and with Gaussian's R/U/RO prefix removed."""
    keys = [canonicalize_level(level)]
    for prefix in ("RO", "R", "U"):
        if level.upper().startswith(prefix):
            keys.append(canonicalize_level(level[len(prefix) :]))
    return keys


def find_scaling_factor(level, scale_type="zpe"):
    """Look up the vibrational scaling factor for a level of theory.

    The level string (as written by the program, e.g. ``RM062X/MG3S`` or
    ``M06-2X/def2-TZVP``) is matched against the Truhlar database after
    canonicalization and after stripping Gaussian's R/U/RO spin prefix.
    ``scale_type`` selects the ZPE (default), harmonic or fundamental factor.
    Returns (factor, reference) or (None, None) if the level is not listed.
    """
    try:
        factor_field, reference_field = SCALE_TYPES[scale_type]
    except KeyError:
        raise KinisotInputError(
            "unknown scale type %r (choose from %s)" % (scale_type, ", ".join(SCALE_TYPES))
        ) from None
    for key in _candidates(level):
        entry = scaling_data_dict.get(key)
        if entry is not None:
            return getattr(entry, factor_field), scaling_refs[getattr(entry, reference_field)]
    return None, None


@dataclass(frozen=True)
class ScalingChoice:
    """The scaling factor used in a calculation and where it came from.

    ``source`` is 'user' (given on the command line or in the API), 'truhlar'
    (looked up for the detected level of theory) or 'default' (1.0 because
    the level is unknown, not in the database, or differs between files).
    ``scale_type`` is the kind of Truhlar factor ('zpe', 'harm', 'fund').
    """

    factor: float
    source: str
    level: Optional[str] = None
    reference: Optional[str] = None
    messages: Tuple[str, ...] = ()
    scale_type: str = "zpe"


def choose_scaling_factor(inputs, user_factor=None, scale_type="zpe"):
    """Decide the scaling factor for a set of HessianInputs.

    A user-supplied factor wins. Otherwise all inputs must share a level of
    theory that is in the Truhlar database; if not, the factor is 1.0 and the
    messages say why.
    """
    if scale_type not in SCALE_TYPES:
        raise KinisotInputError("unknown scale type %r (choose from %s)" % (scale_type, ", ".join(SCALE_TYPES)))
    if user_factor is not None:
        return ScalingChoice(float(user_factor), "user", scale_type=scale_type)

    levels = {}
    for item in inputs:
        levels.setdefault(item.level_of_theory or "unknown", []).append(item.source)

    if len(levels) > 1:
        messages = ["WARNING: the files were not computed at the same level of theory:"]
        messages += ["   %s: %s" % (level, ", ".join(names)) for level, names in levels.items()]
        messages.append("Unable to assign a vibrational scaling factor; using 1.0 (override with -s)")
        return ScalingChoice(1.0, "default", None, None, tuple(messages), scale_type)

    level = next(iter(levels))
    factor, reference = find_scaling_factor(level, scale_type) if level != "unknown" else (None, None)
    if factor is None:
        return ScalingChoice(
            1.0,
            "default",
            level,
            None,
            ("Unable to find vibrational scaling factor for %s; using value of 1.0 (override with -s)" % level,),
            scale_type,
        )
    return ScalingChoice(
        factor,
        "truhlar",
        level,
        reference,
        (
            "Found vibrational scaling factor %s for %s level of theory (%s factor, Truhlar database v5 via GoodVibes)"
            % (factor, level, {"zpe": "ZPE", "harm": "harmonic", "fund": "fundamental"}[scale_type]),
            "REF: " + reference,
        ),
        scale_type,
    )
