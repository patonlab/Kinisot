#!/usr/bin/python

# Comments and/or additions are welcome (send e-mail to:
# robert.paton@colostate.edu

import numpy as np

# File parsing is delegated to GoodVibes (same research group), which reads
# the Cartesian Hessian and isotope-aware per-atom masses from Gaussian
# archive blocks and ORCA .hess files, and provides program-agnostic
# level-of-theory / linearity detection.
from goodvibes.io import parse_hessian, parse_qcdata
from goodvibes.io import level_of_theory  # noqa: F401 (re-exported for Kinisot.py)


def substitute_isotopes(mass_list, iso):
    # Swap the requested atoms (1-based numbers in the comma-separated iso
    # string; '0' means no substitution) for their heavier isotopes.
    # Substitution considers 1H/2H, 12C/13C and 16O/17O. Elements are
    # identified by a mass window so that both Gaussian's isotopic masses
    # (12.00000) and ORCA's average atomic masses (12.011) are recognized.
    # More can be added, but it hasn't been necessary so far...
    mass_list = list(mass_list)
    for atom in iso.split(","):
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
        if 0.9 < mass_list[i] < 1.5:
            mass_list[i] = 2.0141  # 1H - 2D
        if 11.5 < mass_list[i] < 12.5:
            mass_list[i] = 13.00335  # 12C - 13C
        if 15.5 < mass_list[i] < 16.5:
            mass_list[i] = 16.9991  # 16O - 17O
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
        mass_list = substitute_isotopes(hess_data.masses, iso)
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
