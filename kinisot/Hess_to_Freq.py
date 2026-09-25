#!/usr/bin/python
"""Read Gaussian frequency output and build mass-weighted Hessians.

Comments and/or additions are welcome (send e-mail to
robert.paton@colostate.edu).
"""

import re
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .exceptions import KinisotInputError, KinisotParseError
from .isotopes import SUBSTITUTIONS, element_symbol, supported_substitutions


@dataclass(frozen=True)
class FrequencyData:
    """Everything Kinisot needs from one frequency calculation.

    ``hessian`` holds the Cartesian force constants (3N x 3N, Hartree/Bohr^2)
    and ``masses`` the per-atom masses the program used (amu), which respect
    isotope keywords given to the program itself.
    """

    file: str
    atomic_numbers: Tuple[int, ...]
    masses: Tuple[float, ...]
    hessian: np.ndarray
    level_of_theory: Optional[str]
    linear: bool

    @property
    def natoms(self):
        return len(self.masses)


@dataclass(frozen=True)
class Substitution:
    """One isotopic substitution applied to an atom of a file."""

    file: str
    atom: int  # 1-based atom number, as on the command line
    symbol: str
    light_mass: float
    heavy_mass: float

    def __str__(self):
        return "%s (%s atom %d)" % (self.symbol, self.file, self.atom)


def _read_lines(file):
    try:
        with open(file, encoding="utf-8", errors="replace") as handle:
            return handle.readlines()
    except OSError as err:
        raise KinisotParseError("cannot read %s: %s" % (file, err.strerror or err)) from None


def _linear_from_line(line):
    """Linear molecules have a zero (or absent) first rotational constant,
    e.g. ``Rotational constants (GHZ): 0.00000 11.69 11.69``."""
    constants = []
    for token in line.split(":", 1)[1].split():
        try:
            constants.append(float(token))
        except ValueError:
            pass  # Gaussian prints ******** for an infinite constant
    return len(constants) < 3 or min(abs(c) for c in constants) < 1e-4


def _find_archive(lines):
    """Return the text of the last archive entry (joined, unwrapped) or None."""
    start = end = None
    for i, raw in enumerate(lines):
        line = raw.strip()
        if line.startswith("1\\1\\") or line.startswith("1|1|"):
            start, end = i, None
        elif start is not None and end is None and line.endswith("@"):
            end = i
    if start is None or end is None:
        return None
    archive = "".join(l.strip() for l in lines[start : end + 1])
    if archive.startswith("1|1|"):
        # Windows builds of Gaussian use | instead of \ as the separator
        archive = archive.replace("|", "\\")
    return archive


def _level_from_archive(archive):
    fields = archive.split("\\")
    if len(fields) > 5 and fields[4] and fields[5]:
        return fields[4] + "/" + fields[5]
    return None


def parse_gaussian(file):
    """Parse a normally terminated Gaussian frequency job into FrequencyData.

    The force constants are read from the archive entry at the end of the
    output (the ``NImag=`` section), which Gaussian writes only on normal
    termination of a ``freq`` job. Per-atom masses come from the
    ``Atom N has atomic number Z and mass M`` lines of the same job.

    Raises KinisotParseError when anything needed is missing.
    """
    lines = _read_lines(file)
    natoms = None
    atomic_numbers, masses = [], []
    linear = False

    for raw in lines:
        line = raw.strip()
        if line.startswith("NAtoms="):
            try:
                natoms = int(line.split()[1])
            except (IndexError, ValueError):
                raise KinisotParseError("%s: cannot read the atom count from '%s'" % (file, line)) from None
        elif line.startswith("Atom") and "has atomic number" in line and "and mass" in line:
            tokens = line.split()
            try:
                atomic_numbers.append(int(tokens[5]))
                masses.append(float(tokens[8]))
            except (IndexError, ValueError):
                raise KinisotParseError("%s: cannot read atomic number and mass from '%s'" % (file, line)) from None
        elif "Rotational constants (GHZ):" in line:
            linear = _linear_from_line(line)

    if natoms is None:
        raise KinisotParseError("%s: no 'NAtoms=' line found; is this a Gaussian output file?" % file)

    archive = _find_archive(lines)
    if archive is None or "NImag" not in archive:
        raise KinisotParseError(
            "%s: no archive entry with force constants found. Kinisot needs a normally "
            "terminated Gaussian frequency job (the Hessian is read from the archive block "
            "at the end of the output)" % file
        )
    try:
        triangle = [float(x) for x in archive.split("NImag")[1].split("\\")[2].split(",")]
    except (IndexError, ValueError):
        raise KinisotParseError("%s: cannot parse the force constants in the archive entry" % file) from None

    dof = 3 * natoms
    if len(triangle) != dof * (dof + 1) // 2:
        raise KinisotParseError(
            "%s: the archive entry holds %d force constants but %d atoms need %d"
            % (file, len(triangle), natoms, dof * (dof + 1) // 2)
        )
    if len(masses) < natoms:
        raise KinisotParseError("%s: found masses for %d atoms but NAtoms=%d" % (file, len(masses), natoms))
    # Composite jobs print several mass blocks; the last one belongs to the frequency job
    atomic_numbers, masses = atomic_numbers[-natoms:], masses[-natoms:]

    hessian = np.zeros((dof, dof))
    hessian[np.tril_indices(dof)] = triangle
    hessian = hessian + hessian.T - np.diag(np.diag(hessian))

    return FrequencyData(
        file=file,
        atomic_numbers=tuple(atomic_numbers),
        masses=tuple(masses),
        hessian=hessian,
        level_of_theory=_level_from_archive(archive),
        linear=linear,
    )


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
                "--iso label '%s' for %s: 0 means 'no substitution' and cannot be combined "
                "with atom numbers" % (label, file)
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
    """Apply the heavy-isotope substitutions requested by ``label``.

    Returns ``(masses, substitutions)``: the per-atom mass list with the
    substituted atoms replaced by their heavy isotope, and a list of
    Substitution records describing what was changed. Every requested atom
    must be an element Kinisot can substitute and must still carry the light
    isotope mass, otherwise KinisotInputError is raised.
    """
    masses = list(data.masses)
    applied = []
    for idx in parse_label(label, data.natoms, data.file):
        symbol = element_symbol(data.atomic_numbers[idx])
        entry = SUBSTITUTIONS.get(symbol)
        if entry is None:
            raise KinisotInputError(
                "%s: atom %d is %s, for which Kinisot has no isotopic substitution "
                "(supported: %s)" % (data.file, idx + 1, symbol, supported_substitutions())
            )
        light_name, light_mass, heavy_name, heavy_mass = entry
        if not np.isclose(masses[idx], light_mass):
            raise KinisotInputError(
                "%s: atom %d (%s) has mass %.5f, not the %s mass %.5f expected before "
                "substitution with %s; was an isotope already set in the Gaussian input?"
                % (data.file, idx + 1, symbol, masses[idx], light_name, light_mass, heavy_name)
            )
        masses[idx] = heavy_mass
        applied.append(Substitution(data.file, idx + 1, symbol, light_mass, heavy_mass))
    return masses, applied


def mass_weight(hessian, masses):
    """Mass-weight a Cartesian Hessian: H_ij / sqrt(m_i m_j)."""
    weights = np.repeat(np.asarray(masses, dtype=float) ** -0.5, 3)
    return hessian * weights[:, None] * weights[None, :]


def read_hess(file, iso):
    """Mass-weighted Hessian of ``file`` with the substitutions in ``iso`` applied.

    Kept for backwards compatibility; new code should use parse_gaussian(),
    substitute() and mass_weight() so that the substitutions are available.
    """
    data = parse_gaussian(file)
    masses, _ = substitute(data, iso)
    return mass_weight(data.hessian, masses)


def level_of_theory(file):
    """Level of theory and basis set of the archived job, e.g. 'RB3LYP/6-31G(d)'.

    Returns None when the file has no archive entry (single point, truncated
    or crashed job).
    """
    archive = _find_archive(_read_lines(file))
    return _level_from_archive(archive) if archive is not None else None


def is_linear(file):
    """'linear' or 'none', from the rotational constants in the Gaussian output.

    This affects the number of external (rotational) degrees of freedom to
    remove: 2 rather than 3.
    """
    linear = False
    for line in _read_lines(file):
        if "Rotational constants (GHZ):" in line:
            linear = _linear_from_line(line)
    return "linear" if linear else "none"
