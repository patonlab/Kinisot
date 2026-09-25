"""Isotope masses used for substitutions.

Kinisot 2.1 substitutes one heavy isotope per element, selected by atom
number on the command line. Masses are those of the pure isotopes as used
by Gaussian for the light species (e.g. 12C = 12.00000), so that a
substitution changes only the atoms requested. Phase 7 of the
implementation plan replaces this table with a full isotope list and an
explicit ``--iso 5:13C`` syntax.
"""

ELEMENT_SYMBOLS = (
    'X',
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
)

# element symbol -> (light isotope, light mass / amu, heavy isotope, heavy mass / amu)
# The light masses must match what the QC program used for the unsubstituted
# atom; the heavy masses are the values Kinisot has always used.
SUBSTITUTIONS = {
    'H': ('1H', 1.00783, '2H', 2.0141),
    'C': ('12C', 12.00000, '13C', 13.00335),
    'O': ('16O', 15.99491, '17O', 16.9991),
}


def element_symbol(atomic_number):
    """Return the element symbol for an atomic number (``'X'`` if unknown)."""
    try:
        return ELEMENT_SYMBOLS[int(atomic_number)]
    except (IndexError, ValueError):
        return 'X'


def supported_substitutions():
    """Human-readable list of the substitutions Kinisot can make."""
    return ', '.join('%s -> %s' % (light, heavy)
                     for light, _, heavy, _ in SUBSTITUTIONS.values())
