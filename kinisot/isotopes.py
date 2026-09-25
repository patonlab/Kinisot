"""Isotope masses and the substitution of atoms by their heavy isotopes.

Kinisot 2.x substitutes one heavy isotope per element, selected by atom
number on the command line. Masses are those of the pure isotopes as used
by Gaussian for the light species (e.g. 12C = 12.00000), so that a
substitution changes only the atoms requested. Phase 7 of the
implementation plan replaces this table with a full isotope list and an
explicit ``--iso 5:13C`` syntax.
"""

import re
from dataclasses import dataclass

import numpy as np

from .exceptions import KinisotInputError

ELEMENT_SYMBOLS = (
    "X",
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)  # fmt: skip

# element symbol -> (light isotope, light mass / amu, heavy isotope, heavy mass / amu)
# The light masses must match what the QC program used for the unsubstituted
# atom; the heavy masses are the values Kinisot has always used.
SUBSTITUTIONS = {
    "H": ("1H", 1.00783, "2H", 2.0141),
    "C": ("12C", 12.00000, "13C", 13.00335),
    "O": ("16O", 15.99491, "17O", 16.9991),
}


def element_symbol(atomic_number):
    """Return the element symbol for an atomic number (``'X'`` if unknown)."""
    try:
        return ELEMENT_SYMBOLS[int(atomic_number)]
    except (IndexError, ValueError):
        return "X"


def supported_substitutions():
    """Human-readable list of the substitutions Kinisot can make."""
    return ", ".join("%s -> %s" % (light, heavy) for light, _, heavy, _ in SUBSTITUTIONS.values())


@dataclass(frozen=True)
class Substitution:
    """One isotopic substitution applied to an atom of a species."""

    file: str
    atom: int  # 1-based atom number, as on the command line
    symbol: str
    light_mass: float
    heavy_mass: float

    def __str__(self):
        return "%s (%s atom %d)" % (self.symbol, self.file, self.atom)

    def to_dict(self):
        return {
            "atom": self.atom,
            "element": self.symbol,
            "light_mass": self.light_mass,
            "heavy_mass": self.heavy_mass,
        }


_LABEL_SPLIT = re.compile(r"[,\s]+")


def parse_label(label, natoms, file):
    """Turn an ``--iso`` label ('0', '5', '7,8') into 0-based atom indices."""
    text = str(label).strip()
    if text in ("", "0"):
        return []
    indices = []
    for token in _LABEL_SPLIT.split(text):
        if not token:
            continue
        try:
            atom = int(token)
        except ValueError:
            raise KinisotInputError(
                "--iso label '%s' for %s: '%s' is not an atom number" % (label, file, token)
            ) from None
        if atom == 0:
            raise KinisotInputError(
                "--iso label '%s' for %s: 0 means 'no substitution' and cannot be combined with atom numbers"
                % (label, file)
            )
        if atom < 1 or atom > natoms:
            raise KinisotInputError(
                "--iso label '%s' for %s: atom %d is out of range (the file has %d atoms)" % (label, file, atom, natoms)
            )
        if atom - 1 in indices:
            raise KinisotInputError("--iso label '%s' for %s: atom %d is listed twice" % (label, file, atom))
        indices.append(atom - 1)
    return indices


def substitute(data, label):
    """Apply the heavy-isotope substitutions requested by ``label`` to ``data``.

    ``data`` is a :class:`~kinisot.hessian.HessianInput` (anything with
    ``masses``, ``atomic_numbers`` and ``source``). Returns
    ``(masses, substitutions)``: the per-atom mass list with the substituted
    atoms replaced by their heavy isotope, and a list of Substitution records.
    Every requested atom must be an element Kinisot can substitute and must
    still carry the light isotope mass, otherwise KinisotInputError is raised.
    """
    masses = list(data.masses)
    applied = []
    for idx in parse_label(label, len(masses), data.source):
        symbol = element_symbol(data.atomic_numbers[idx])
        entry = SUBSTITUTIONS.get(symbol)
        if entry is None:
            raise KinisotInputError(
                "%s: atom %d is %s, for which Kinisot has no isotopic substitution (supported: %s)"
                % (data.source, idx + 1, symbol, supported_substitutions())
            )
        light_name, light_mass, heavy_name, heavy_mass = entry
        if not np.isclose(masses[idx], light_mass):
            raise KinisotInputError(
                "%s: atom %d (%s) has mass %.5f, not the %s mass %.5f expected before substitution with %s; "
                "was an isotope already set in the Gaussian input?"
                % (data.source, idx + 1, symbol, masses[idx], light_name, light_mass, heavy_name)
            )
        masses[idx] = heavy_mass
        applied.append(Substitution(data.source, idx + 1, symbol, light_mass, heavy_mass))
    return masses, applied
