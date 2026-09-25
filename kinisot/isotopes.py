"""Isotope masses and the substitution of atoms by heavier isotopes.

Both isotopologues of every species are built from the isotope table in
:mod:`kinisot.isotope_data` (AME 2020 masses): the light one uses the most
abundant isotope of every element, which is also the convention of Gaussian
(¹²C 12.000, ¹H 1.00783), and the heavy one replaces the requested atoms.
The masses a program reports are used only to check that the atoms were not
already substituted in the program's own input.

Label syntax (``--iso``)::

    5           atom 5 with its default heavy isotope (H->2H, C->13C, N->15N, O->18O, ...)
    5:13C       atom 5 as carbon-13 (explicit isotope; must match the atom's element)
    7:D, 7:T    deuterium and tritium shorthands
    5:13.5      an explicit mass in amu
    7,8 / 7 8   several atoms, comma or space separated
    0           no substitution in this file
"""

import re
from dataclasses import dataclass

from .exceptions import KinisotInputError
from .isotope_data import ISOTOPE_MASSES, MOST_ABUNDANT

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

# Heavy isotope meant by a bare atom number, for the elements commonly labelled.
# Other elements need an explicit isotope (e.g. 5:29Si).
DEFAULT_HEAVY = {"H": 2, "C": 13, "N": 15, "O": 18, "S": 34, "Cl": 37, "Br": 81, "Si": 29}

# Shorthands accepted after the colon
ISOTOPE_ALIASES = {"D": ("H", 2), "T": ("H", 3)}

_ELEMENT_INDEX = {symbol.upper(): symbol for symbol in ELEMENT_SYMBOLS[1:]}
_LABEL_SPLIT = re.compile(r"[,\s]+")
_ISOTOPE = re.compile(r"^(\d{1,3})([A-Za-z]{1,2})$")


def element_symbol(atomic_number):
    """Return the element symbol for an atomic number (``'X'`` if unknown)."""
    try:
        return ELEMENT_SYMBOLS[int(atomic_number)]
    except (IndexError, ValueError):
        return "X"


def light_mass(symbol):
    """Mass of the most abundant isotope of ``symbol`` (the light isotopologue)."""
    try:
        return ISOTOPE_MASSES[symbol][MOST_ABUNDANT[symbol]]
    except KeyError:
        raise KinisotInputError("no isotope data for element %r" % symbol) from None


def isotope_mass(symbol, mass_number):
    """Mass of the isotope ``mass_number`` of ``symbol``."""
    try:
        return ISOTOPE_MASSES[symbol][int(mass_number)]
    except KeyError:
        known = ISOTOPE_MASSES.get(symbol, {})
        raise KinisotInputError(
            "no mass for %d%s in the isotope table (known: %s)"
            % (int(mass_number), symbol, ", ".join("%d%s" % (a, symbol) for a in sorted(known)) or "none")
        ) from None


def default_heavy_isotope(symbol):
    """Mass number substituted for a bare atom number, or None if the element has no default."""
    return DEFAULT_HEAVY.get(symbol)


def supported_substitutions():
    """Human-readable list of the default substitutions."""
    return ", ".join("%s -> %d%s" % (symbol, a, symbol) for symbol, a in DEFAULT_HEAVY.items())


@dataclass(frozen=True)
class Substitution:
    """One isotopic substitution applied to an atom of a species."""

    file: str
    atom: int  # 1-based atom number, as on the command line
    symbol: str
    light_mass: float
    heavy_mass: float
    isotope: str = ""  # e.g. '13C', '2H', or 'm=13.5' for an explicit mass

    def __str__(self):
        return "%s (%s atom %d%s)" % (self.symbol, self.file, self.atom, " as " + self.isotope if self.isotope else "")

    def to_dict(self):
        return {
            "atom": self.atom,
            "element": self.symbol,
            "isotope": self.isotope,
            "light_mass": self.light_mass,
            "heavy_mass": self.heavy_mass,
        }


def _parse_isotope_spec(spec, label, file):
    """'13C' -> ('C', 13, None); 'D' -> ('H', 2, None); '13.5' -> (None, None, 13.5)."""
    text = spec.strip()
    if text.upper() in ISOTOPE_ALIASES:
        symbol, mass_number = ISOTOPE_ALIASES[text.upper()]
        return symbol, mass_number, None
    match = _ISOTOPE.match(text)
    if match:
        symbol = _ELEMENT_INDEX.get(match.group(2).upper())
        if symbol is None:
            raise KinisotInputError("--iso label '%s' for %s: unknown element in '%s'" % (label, file, text))
        return symbol, int(match.group(1)), None
    if "." in text:
        try:
            return None, None, float(text)
        except ValueError:
            pass
    raise KinisotInputError(
        "--iso label '%s' for %s: cannot read the isotope '%s' (use e.g. 13C, 18O, D, T, or a mass like 13.5)"
        % (label, file, text)
    )


def parse_label(label, natoms, file):
    """Turn an ``--iso`` label into a list of ``(index, symbol, mass_number, mass)``.

    ``index`` is 0-based; ``symbol``/``mass_number`` are None for a bare atom
    number (default heavy isotope) and ``mass`` is set only for an explicit
    mass. Returns [] for '0'.
    """
    text = str(label).strip()
    if text in ("", "0"):
        return []
    entries = []
    seen = set()
    for token in _LABEL_SPLIT.split(text):
        if not token:
            continue
        atom_text, _, spec = token.partition(":")
        try:
            atom = int(atom_text)
        except ValueError:
            raise KinisotInputError(
                "--iso label '%s' for %s: '%s' is not an atom number" % (label, file, atom_text)
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
        if atom in seen:
            raise KinisotInputError("--iso label '%s' for %s: atom %d is listed twice" % (label, file, atom))
        seen.add(atom)
        symbol, mass_number, mass = _parse_isotope_spec(spec, label, file) if spec else (None, None, None)
        entries.append((atom - 1, symbol, mass_number, mass))
    return entries


def light_masses(data):
    """Light-isotopologue masses for a HessianInput (most abundant isotope of every atom).

    Checks that the masses the program reported are consistent with those
    (an atom already substituted in the program's input is refused).
    """
    masses = []
    for idx, (z, reported) in enumerate(zip(data.atomic_numbers, data.masses)):
        symbol = element_symbol(z)
        mass = light_mass(symbol)
        if abs(reported - mass) > max(0.5, 0.015 * mass):
            raise KinisotInputError(
                "%s: atom %d (%s) has mass %.5f in the program's output but the light isotope %d%s weighs %.5f; "
                "Kinisot builds both isotopologues itself, so do not substitute isotopes in the program's input"
                % (data.source, idx + 1, symbol, reported, MOST_ABUNDANT[symbol], symbol, mass)
            )
        masses.append(mass)
    return masses


def substitute(data, label, warnings_out=None):
    """Apply the substitutions requested by ``label`` to the light masses of ``data``.

    ``data`` is a :class:`~kinisot.hessian.HessianInput`. Returns
    ``(masses, substitutions)``: the per-atom mass list of the isotopologue
    and a list of Substitution records. Every requested atom must be an
    element with a default heavy isotope (or carry an explicit isotope that
    matches its element), otherwise KinisotInputError is raised.
    """
    masses = light_masses(data)
    applied = []
    for idx, symbol, mass_number, mass in parse_label(label, len(masses), data.source):
        element = element_symbol(data.atomic_numbers[idx])
        if mass is not None:
            heavy, isotope = float(mass), "m=%g" % mass
        else:
            if symbol is None:
                mass_number = default_heavy_isotope(element)
                if mass_number is None:
                    raise KinisotInputError(
                        "%s: atom %d is %s, which has no default heavy isotope; give one explicitly, e.g. %d:%d%s "
                        "(defaults: %s)"
                        % (
                            data.source,
                            idx + 1,
                            element,
                            idx + 1,
                            max(ISOTOPE_MASSES.get(element, {MOST_ABUNDANT.get(element, 0): 0})),
                            element,
                            supported_substitutions(),
                        )
                    )
                if element == "O" and warnings_out is not None:
                    warnings_out.append(
                        "%s: atom %d is oxygen and a bare atom number now means 18O (Kinisot 2.3 and earlier used "
                        "17O); write %d:17O for the old behaviour. This note disappears in the next release."
                        % (data.source, idx + 1, idx + 1)
                    )
            elif symbol != element:
                raise KinisotInputError(
                    "%s: atom %d is %s but the label asks for %d%s"
                    % (data.source, idx + 1, element, mass_number, symbol)
                )
            heavy, isotope = isotope_mass(element, mass_number), "%d%s" % (mass_number, element)
        if abs(heavy - masses[idx]) < 1e-9:
            raise KinisotInputError(
                "%s: atom %d (%s) is already %s; a substitution must change the mass"
                % (data.source, idx + 1, element, isotope)
            )
        applied.append(Substitution(data.source, idx + 1, element, masses[idx], heavy, isotope))
        masses[idx] = heavy
    return masses, applied
