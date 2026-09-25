"""Read ORCA frequency jobs into a HessianInput.

ORCA writes the Hessian to ``name.hess`` next to ``name.out``; either path is
accepted and the other is located by name. The ``$hessian`` block is read
with GoodVibes' parser, the ``$atoms`` block (symbols, masses, Bohr
coordinates in the Hessian's frame) and ``$vibrational_frequencies`` here,
and the level of theory from the output with GoodVibes' ORCA-aware
``level_of_theory``.

Masses: ORCA reports standard atomic weights (C 12.011) rather than pure
isotopes. For the elements Kinisot can substitute the light isotopologue is
built from the pure-isotope masses of :mod:`kinisot.isotopes` so that the
same Hessian gives the same numbers whichever program produced it; other
elements keep ORCA's masses until the full isotope table of Phase 7.
"""

import os

from goodvibes.io import level_of_theory as _gv_level_of_theory
from goodvibes.io import parse_hessian as _gv_parse_hessian

from ..exceptions import KinisotParseError
from ..hessian import HessianInput, linear_from_geometry
from ..isotopes import ELEMENT_SYMBOLS, SUBSTITUTIONS

__all__ = ["parse_orca", "orca_paths"]

_ATOMIC_NUMBERS = {symbol.upper(): z for z, symbol in enumerate(ELEMENT_SYMBOLS) if z}


def orca_paths(file):
    """(output path or None, .hess path) for an ORCA job given either file."""
    stub, ext = os.path.splitext(file)
    if ext == ".hess":
        hess = file
        out = next((stub + e for e in (".out", ".log") if os.path.exists(stub + e)), None)
    else:
        out = file
        hess = stub + ".hess"
        if not os.path.exists(hess):
            raise KinisotParseError(
                "%s: ORCA stores the Hessian in a separate file; expected %s next to the output" % (file, hess)
            )
    return out, hess


def _read_atoms_and_frequencies(hess_path):
    try:
        with open(hess_path, encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()
    except OSError as err:
        raise KinisotParseError("cannot read %s: %s" % (hess_path, err.strerror or err)) from None
    symbols, masses, positions, frequencies = [], [], [], []
    i = 0
    while i < len(lines):
        key = lines[i].strip()
        if key == "$atoms":
            n = int(lines[i + 1].split()[0])
            for line in lines[i + 2 : i + 2 + n]:
                tokens = line.split()
                symbols.append(tokens[0])
                masses.append(float(tokens[1]))
                positions.append([float(x) for x in tokens[2:5]])
            i += 2 + n
        elif key == "$vibrational_frequencies":
            n = int(lines[i + 1].split()[0])
            frequencies = [float(line.split()[1]) for line in lines[i + 2 : i + 2 + n]]
            i += 2 + n
        else:
            i += 1
    if not symbols:
        raise KinisotParseError("%s: no $atoms section found" % hess_path)
    return symbols, masses, positions, frequencies


def parse_orca(file):
    """Parse an ORCA frequency job (``.out`` or ``.hess`` path) into a HessianInput."""
    out, hess_path = orca_paths(file)
    try:
        data = _gv_parse_hessian(hess_path)
    except (ValueError, OSError) as err:
        raise KinisotParseError("%s: %s" % (hess_path, err)) from None
    symbols, orca_masses, positions, frequencies = _read_atoms_and_frequencies(hess_path)
    if data.hessian.shape[0] != 3 * len(symbols):
        raise KinisotParseError(
            "%s: $hessian is %dx%d but $atoms lists %d atoms" % (hess_path, *data.hessian.shape, len(symbols))
        )

    atomic_numbers = []
    masses = []
    for symbol, mass in zip(symbols, orca_masses):
        z = _ATOMIC_NUMBERS.get(symbol.upper())
        if z is None:
            raise KinisotParseError("%s: unknown element symbol %r in $atoms" % (hess_path, symbol))
        atomic_numbers.append(z)
        entry = SUBSTITUTIONS.get(ELEMENT_SYMBOLS[z])
        masses.append(entry[1] if entry is not None else mass)

    level = None
    if out is not None:
        try:
            level = _gv_level_of_theory(out)
        except (OSError, ValueError, IndexError):
            level = None
        if level in ("none/none", "", None):
            level = None
        elif level.endswith("/none"):
            level = level[: -len("/none")]

    vibrational = [f for f in frequencies if abs(f) > 1e-6] if frequencies else []
    return HessianInput(
        hessian=data.hessian,
        masses=tuple(masses),
        atomic_numbers=tuple(atomic_numbers),
        source=out or hess_path,
        program="Orca",
        level_of_theory=level,
        linear=linear_from_geometry(positions, masses),
        positions=positions,
        program_frequencies=tuple(vibrational) if vibrational else None,
    )
