"""The program-independent input to an isotope-effect calculation.

Every backend (Gaussian today; ORCA and ASE in later phases) produces a
:class:`HessianInput`. The physics in :mod:`kinisot.api` only ever sees this
type, so adding a program means adding a parser, not touching the
Bigeleisen-Mayer code.
"""

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .isotopes import element_symbol


@dataclass(frozen=True)
class HessianInput:
    """Cartesian Hessian plus the atom data needed to mass-weight it.

    Attributes
    ----------
    hessian : (3N, 3N) array, Hartree/Bohr^2, symmetric.
    masses : per-atom masses in amu of the *light* isotopologue, i.e. the
        masses the program used (pure most-abundant isotopes in Gaussian).
    atomic_numbers : per-atom atomic numbers.
    source : where the data came from (file path or a description).
    program : 'Gaussian', 'Orca', 'ase', ... (informational).
    level_of_theory : e.g. 'RB3LYP/6-31G(d)', used to look up a scaling
        factor; None when unknown (machine-learned potentials).
    linear : whether the molecule is linear (5 rather than 6 external modes).
    positions : optional (N, 3) Cartesian coordinates in Bohr, in the same
        frame as the Hessian (needed for projection of external modes).
    program_frequencies : optional vibrational frequencies (cm-1, unscaled,
        negative for imaginary) as printed by the program, used only to
        check that Kinisot reproduces them for the unsubstituted species.
    energy : optional electronic energy in Hartree (for the Skodje-Truhlar
        barrier height).
    """

    hessian: np.ndarray
    masses: Tuple[float, ...]
    atomic_numbers: Tuple[int, ...]
    source: str = ""
    program: str = ""
    level_of_theory: Optional[str] = None
    linear: bool = False
    positions: Optional[np.ndarray] = None
    program_frequencies: Optional[Tuple[float, ...]] = None
    energy: Optional[float] = None

    def __post_init__(self):
        hessian = np.array(self.hessian, dtype=float)
        masses = tuple(float(m) for m in self.masses)
        atomic_numbers = tuple(int(z) for z in self.atomic_numbers)
        n = len(masses)
        if len(atomic_numbers) != n:
            raise ValueError("%d masses but %d atomic numbers" % (n, len(atomic_numbers)))
        if hessian.shape != (3 * n, 3 * n):
            raise ValueError("Hessian has shape %s but %d atoms need (%d, %d)" % (hessian.shape, n, 3 * n, 3 * n))
        if any(m <= 0 for m in masses):
            raise ValueError("all masses must be positive")
        positions = self.positions
        if positions is not None:
            positions = np.array(positions, dtype=float)
            if positions.shape != (n, 3):
                raise ValueError("positions have shape %s but %d atoms need (%d, 3)" % (positions.shape, n, n))
        object.__setattr__(self, "hessian", hessian)
        object.__setattr__(self, "masses", masses)
        object.__setattr__(self, "atomic_numbers", atomic_numbers)
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "linear", bool(self.linear))
        if self.program_frequencies is not None:
            object.__setattr__(self, "program_frequencies", tuple(float(f) for f in self.program_frequencies))

    @property
    def natoms(self):
        return len(self.masses)

    @property
    def symbols(self):
        return tuple(element_symbol(z) for z in self.atomic_numbers)

    @property
    def file(self):
        """Alias of ``source`` (kept for the 2.1 FrequencyData name)."""
        return self.source

    def mass_weighted(self, masses=None):
        """Mass-weighted Hessian H_ij / sqrt(m_i m_j) with ``masses`` (default: the light ones)."""
        return mass_weight(self.hessian, self.masses if masses is None else masses)


def mass_weight(hessian, masses):
    """Mass-weight a Cartesian Hessian: H_ij / sqrt(m_i m_j)."""
    weights = np.repeat(np.asarray(masses, dtype=float) ** -0.5, 3)
    return np.asarray(hessian) * weights[:, None] * weights[None, :]


def linear_from_geometry(positions, masses, tolerance=1e-4):
    """Whether a molecule is linear: its smallest principal moment of inertia is ~0.

    ``positions`` in any consistent length unit (N, 3), ``masses`` in amu; the
    test is on the smallest moment relative to the largest.
    """
    positions = np.asarray(positions, dtype=float)
    masses = np.asarray(masses, dtype=float)
    if len(masses) < 3:
        return True
    center = (masses[:, None] * positions).sum(axis=0) / masses.sum()
    r = positions - center
    inertia = np.zeros((3, 3))
    for m, (x, y, z) in zip(masses, r):
        inertia += m * np.array(
            [[y * y + z * z, -x * y, -x * z], [-x * y, x * x + z * z, -y * z], [-x * z, -y * z, x * x + y * y]]
        )
    moments = np.linalg.eigvalsh(inertia)
    return bool(moments[0] < tolerance * max(moments[-1], 1e-300))
