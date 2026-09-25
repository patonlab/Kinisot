"""Choosing the vibrational scaling factor for a calculation."""

from dataclasses import dataclass
from typing import Optional, Tuple

from .vib_scale_factors import REFERENCES, SCALING_FACTORS

__all__ = ["find_scaling_factor", "choose_scaling_factor", "ScalingChoice"]


def _normalize(level):
    return level.upper().replace("-", "")


_INDEX = {_normalize(level): level for level in SCALING_FACTORS}


def find_scaling_factor(level):
    """Look up the ZPE vibrational scaling factor for a level of theory.

    The level string (as written in the Gaussian archive, e.g. RM062X/MG3S)
    is matched exactly against the Truhlar database entries after normalizing
    case and hyphens and stripping Gaussian's R/U/RO spin prefix. Returns
    (factor, reference) or (None, None) if the level is not in the database.
    """
    candidates = [_normalize(level)]
    for prefix in ("RO", "R", "U"):
        if level.upper().startswith(prefix):
            candidates.append(_normalize(level[len(prefix) :]))
    for candidate in candidates:
        key = _INDEX.get(candidate)
        if key is not None:
            entry = SCALING_FACTORS[key]
            return entry.zpe, REFERENCES[entry.zpe_ref]
    return None, None


@dataclass(frozen=True)
class ScalingChoice:
    """The scaling factor used in a calculation and where it came from.

    ``source`` is 'user' (given on the command line or in the API), 'truhlar'
    (looked up for the detected level of theory) or 'default' (1.0 because
    the level is unknown, not in the database, or differs between files).
    """

    factor: float
    source: str
    level: Optional[str] = None
    reference: Optional[str] = None
    messages: Tuple[str, ...] = ()


def choose_scaling_factor(inputs, user_factor=None):
    """Decide the scaling factor for a set of HessianInputs.

    A user-supplied factor wins. Otherwise all inputs must share a level of
    theory that is in the Truhlar database; if not, the factor is 1.0 and the
    messages say why.
    """
    if user_factor is not None:
        return ScalingChoice(float(user_factor), "user")

    levels = {}
    for item in inputs:
        levels.setdefault(item.level_of_theory or "unknown", []).append(item.source)

    if len(levels) > 1:
        messages = ["WARNING: the files were not computed at the same level of theory:"]
        messages += ["   %s: %s" % (level, ", ".join(names)) for level, names in levels.items()]
        messages.append("Unable to assign a vibrational scaling factor; using 1.0 (override with -s)")
        return ScalingChoice(1.0, "default", None, None, tuple(messages))

    level = next(iter(levels))
    factor, reference = find_scaling_factor(level) if level != "unknown" else (None, None)
    if factor is None:
        return ScalingChoice(
            1.0,
            "default",
            level,
            None,
            ("Unable to find vibrational scaling factor for %s; using value of 1.0 (override with -s)" % level,),
        )
    return ScalingChoice(
        factor,
        "truhlar",
        level,
        reference,
        (
            "Found vibrational scaling factor %s for %s level of theory (ZPE factor, Truhlar database)"
            % (factor, level),
            "REF: " + reference,
        ),
    )
