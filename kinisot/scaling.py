"""Vibrational frequency scaling-factor lookup.

Scaling factors come from the Truhlar group database (v5) bundled with
GoodVibes; level-of-theory detection comes from GoodVibes' parsers.
"""

from goodvibes.io import level_of_theory
from goodvibes.vib_scale_factors import (
    canonicalize_level,
    scaling_data_dict,
    scaling_refs,
)


def find_scaling_factor(level):
    """
    Look up the ZPE vibrational scaling factor for a level of theory in the
    Truhlar group database (v5), via GoodVibes. Level strings are matched
    exactly after normalization by goodvibes' canonicalize_level (case and
    cross-program functional aliases); a leading R/U/RO spin prefix, as
    written in Gaussian archives, is also tolerated. Returns
    (factor, reference) or (None, None) if the level is not in the database.
    """
    candidates = [level]
    for prefix in ("RO", "R", "U"):
        if level.upper().startswith(prefix):
            candidates.append(level[len(prefix) :])
    for cand in candidates:
        entry = scaling_data_dict.get(canonicalize_level(cand))
        if entry is not None:
            return entry.zpe_fac, scaling_refs[entry.zpe_ref]
    return None, None


def get_frequency_scaling(files):
    """Detect the level of theory across the input files and look up the
    matching ZPE scaling factor.

    Returns (factor, level, reference, warning): factor is 1.0 with a level
    of "unknown" when the files disagree or the level is not in the
    database; warning carries a message when the files disagree.
    """
    l_o_t = [level_of_theory(file) for file in files]
    if len(set(l_o_t)) > 1:
        warning = "found different levels of theory across input files: " + " / ".join(
            l_o_t
        )
        return 1.0, "unknown", None, warning

    level = l_o_t[0]
    factor, ref = find_scaling_factor(level)
    if factor is None:
        return 1.0, level, None, None
    return factor, level, ref, None
