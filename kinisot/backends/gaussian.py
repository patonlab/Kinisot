"""Read Gaussian frequency output into a HessianInput.

Kinisot needs a normally terminated ``freq`` job: the atom count and per-atom
masses printed by the frequency job, and the archive entry at the end of the
output, which carries the level of theory, the geometry and the
lower-triangular force-constant matrix (the ``NImag=`` section).
"""

import numpy as np

from ..exceptions import KinisotParseError
from ..hessian import HessianInput
from ..thermo import BOHR_TO_ANGSTROM

__all__ = ["parse_gaussian", "level_of_theory", "is_linear"]


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
    archive = "".join(line.strip() for line in lines[start : end + 1])
    if archive.startswith("1|1|"):
        # Windows builds of Gaussian use | instead of \ as the separator
        archive = archive.replace("|", "\\")
    return archive


def _level_from_archive(archive):
    fields = archive.split("\\")
    if len(fields) > 5 and fields[4] and fields[5]:
        return fields[4] + "/" + fields[5]
    return None


def _positions_from_archive(archive, natoms):
    """Cartesian coordinates (Bohr) from the geometry section of the archive, or None.

    The section after the title holds ``charge,multiplicity`` followed by one
    ``Symbol,x,y,z`` (or ``Symbol,0,x,y,z``) entry per atom, in Angstrom.
    """
    sections = archive.split("\\\\")
    for section in sections:
        entries = section.split("\\")
        if len(entries) != natoms + 1:
            continue
        head = entries[0].split(",")
        if len(head) != 2:
            continue
        try:
            int(head[0]), int(head[1])
            coordinates = [[float(x) for x in entry.split(",")[-3:]] for entry in entries[1:]]
        except ValueError:
            continue
        return np.array(coordinates) / BOHR_TO_ANGSTROM
    return None


def parse_gaussian(file):
    """Parse a normally terminated Gaussian frequency job into a HessianInput.

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

    return HessianInput(
        hessian=hessian,
        masses=tuple(masses),
        atomic_numbers=tuple(atomic_numbers),
        source=file,
        program="Gaussian",
        level_of_theory=_level_from_archive(archive),
        linear=linear,
        positions=_positions_from_archive(archive, natoms),
    )


def level_of_theory(file):
    """Level of theory and basis set of the archived job, e.g. 'RB3LYP/6-31G(d)'.

    Returns None when the file has no archive entry (single point, truncated
    or crashed job).
    """
    archive = _find_archive(_read_lines(file))
    return _level_from_archive(archive) if archive is not None else None


def is_linear(file):
    """'linear' or 'none', from the rotational constants in the Gaussian output."""
    linear = False
    for line in _read_lines(file):
        if "Rotational constants (GHZ):" in line:
            linear = _linear_from_line(line)
    return "linear" if linear else "none"
