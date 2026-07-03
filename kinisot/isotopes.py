"""Isotope masses for substitutions.

The light isotopologue is mass-weighted with principal isotope masses
(matching what Gaussian prints); heavy substitutions use the exact isotopic
masses below. Elements are identified from the mass the QC program reported,
via a window wide enough to catch both isotopic (12.00000) and
abundance-averaged (12.011) conventions.
"""

# element symbol -> (window low, window high, principal isotope mass)
ELEMENTS = {
    "H": (0.9, 1.5, 1.00783),
    "C": (11.5, 12.5, 12.00000),
    "N": (13.5, 14.5, 14.00307),
    "O": (15.5, 16.5, 15.99491),
}

# isotope label -> (element it replaces, mass)
HEAVY_ISOTOPES = {
    "2D": ("H", 2.0141),
    "3T": ("H", 3.01605),
    "13C": ("C", 13.00335),
    "14C": ("C", 14.00324),
    "15N": ("N", 15.00011),
    "17O": ("O", 16.9991),
    "18O": ("O", 17.99916),
}

# plain atom numbers (no :isotope suffix) get the traditional default
DEFAULT_HEAVY = {"H": "2D", "C": "13C", "N": "15N", "O": "17O"}


def element_from_mass(mass):
    """Identify the element from the mass the QC program reported, or None."""
    for symbol, (low, high, _principal) in ELEMENTS.items():
        if low < mass < high:
            return symbol
    return None
