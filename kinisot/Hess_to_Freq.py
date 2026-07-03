#!/usr/bin/python

# Comments and/or additions are welcome (send e-mail to:
# robert.paton@colostate.edu

import numpy as np

# File parsing is delegated to GoodVibes (same research group), which reads
# the Cartesian Hessian and isotope-aware per-atom masses from Gaussian
# archive blocks and ORCA .hess files, and provides program-agnostic
# level-of-theory / linearity detection.
from goodvibes.io import parse_hessian, parse_qcdata
from goodvibes.io import level_of_theory  # noqa: F401 (re-exported)

from .isotopes import DEFAULT_HEAVY, ELEMENTS, HEAVY_ISOTOPES, element_from_mass


def normalize_principal_masses(mass_list):
    # The Bigeleisen-Mayer equation compares isotopically pure species, so
    # the light isotopologue must be mass-weighted with principal isotope
    # masses (1H = 1.00783, 12C = 12.00000, ...). Gaussian prints exactly
    # these; ORCA .hess files store abundance-averaged atomic masses
    # (C = 12.011) instead, which would bias KIEs at the 1e-4 level. Map the
    # substitutable elements onto the principal masses; other elements keep
    # the program's value, which cancels almost exactly in the RPFR ratios.
    normalized = []
    for mass in mass_list:
        element = element_from_mass(mass)
        if element is not None:
            mass = ELEMENTS[element][2]
        normalized.append(mass)
    return normalized


def substitute_isotopes(mass_list, iso):
    """Apply the isotopic substitutions in ``iso`` to a list of masses.

    ``iso`` is a comma-separated list of 1-based atom numbers, each with an
    optional ``:isotope`` suffix, e.g. ``"5"`` (default heavy isotope: 2D,
    13C, 15N or 17O depending on the element) or ``"5:18O,7:2D"``. ``"0"``
    means no substitution in this file.
    """
    mass_list = list(mass_list)
    for item in iso.split(","):
        atom, _, isotope = item.partition(":")
        try:
            i = int(atom) - 1
        except ValueError:
            raise ValueError(
                "Invalid atom number {!r} in isotope specification {!r}".format(
                    atom, iso
                )
            )
        if i == -1:
            continue
        if i < -1 or i >= len(mass_list):
            raise ValueError(
                "Atom number {} in isotope specification {!r} is out of range: "
                "the molecule has {} atoms".format(atom, iso, len(mass_list))
            )

        element = element_from_mass(mass_list[i])
        if element is None:
            raise ValueError(
                "No isotope substitution is defined for atom {} (mass {}); "
                "substitutable elements: {}".format(
                    atom, mass_list[i], ", ".join(sorted(ELEMENTS))
                )
            )

        if not isotope:
            isotope = DEFAULT_HEAVY[element]
        if isotope not in HEAVY_ISOTOPES:
            raise ValueError(
                "Unknown isotope {!r} for atom {}; available: {}".format(
                    isotope, atom, ", ".join(sorted(HEAVY_ISOTOPES))
                )
            )
        iso_element, iso_mass = HEAVY_ISOTOPES[isotope]
        if iso_element != element:
            raise ValueError(
                "Isotope {} cannot replace atom {}, which is {} (mass {})".format(
                    isotope, atom, element, mass_list[i]
                )
            )
        mass_list[i] = iso_mass
    return mass_list


def read_hess(file, iso):
    # The Cartesian force constant matrix (Hartree/Bohr^2) and per-atom masses
    # are obtained via GoodVibes; the matrix values are then mass-weighted
    # according to the isotopic masses.
    # Vibrational scaling factors are not applied to matrix elements at this
    # stage: the resulting frequencies can be scaled after diagonalization
    hess_data = parse_hessian(file)
    if not hess_data.masses:
        raise ValueError("No atomic masses found for " + file)

    try:
        mass_list = substitute_isotopes(
            normalize_principal_masses(hess_data.masses), iso
        )
    except ValueError as e:
        raise ValueError("{} (file {})".format(e, file))

    masses_per_coord = np.repeat(mass_list, 3)
    return hess_data.hessian / np.sqrt(np.outer(masses_per_coord, masses_per_coord))


def is_linear(file):
    # Check if a molecule is linear or not, as detected by the QC program
    # This affects the number of rotational d.o.f.
    if parse_qcdata(file).linear_mol:
        return "linear"
    return "none"


def describe_substitutions(mass_list, iso):
    """Human-readable substitution list for an iso spec.

    E.g. "26" on a hydrogen -> ["H26 → 2D"]; "5:18O,7" -> ["O5 → 18O",
    "H7 → 2D"]. "0" (no substitution) gives []. Assumes the spec has
    already been validated by substitute_isotopes.
    """
    parts = []
    for item in iso.split(","):
        atom, _, isotope = item.partition(":")
        index = int(atom) - 1
        if index == -1:
            continue
        element = element_from_mass(mass_list[index])
        if not isotope:
            isotope = DEFAULT_HEAVY[element]
        parts.append("{}{} → {}".format(element, atom, isotope))
    return parts
