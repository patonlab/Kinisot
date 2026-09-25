"""Eckart projection of translations and rotations out of a mass-weighted Hessian.

Quantum-chemistry programs do this before printing frequencies; Kinisot's
default is to diagonalize the raw mass-weighted Hessian and discard the
5/6 lowest modes, which works on tightly converged geometries but not on
noisy finite-difference Hessians. With projection the external modes come
out at exactly zero (to numerical precision) and are removed by value.
"""

import numpy as np

__all__ = ["external_mode_vectors", "project_external_modes"]


def external_mode_vectors(positions, masses, tolerance=1e-8):
    """Orthonormal basis of the translations and rotations in mass-weighted Cartesians.

    ``positions`` (N, 3) in any length unit, ``masses`` (N,) in amu. Returns
    a (3N, n_external) array with n_external = 6, or 5 for a linear molecule
    (the rotation about the molecular axis has no amplitude and drops out).
    """
    positions = np.asarray(positions, dtype=float)
    masses = np.asarray(masses, dtype=float)
    n = len(masses)
    sqrt_m = np.sqrt(masses)
    center = (masses[:, None] * positions).sum(axis=0) / masses.sum()
    r = positions - center
    vectors = []
    for k in range(3):
        v = np.zeros((n, 3))
        v[:, k] = sqrt_m
        vectors.append(v.ravel())
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1.0
        vectors.append((np.cross(axis, r) * sqrt_m[:, None]).ravel())
    basis = np.array(vectors).T  # (3N, 6)
    # Orthonormalize; a linear molecule leaves one rotation with zero norm
    u, s, _ = np.linalg.svd(basis, full_matrices=False)
    keep = s > tolerance * s.max()
    return u[:, keep]


def project_external_modes(mw_hessian, positions, masses):
    """Return (projected mass-weighted Hessian, number of external modes removed)."""
    q = external_mode_vectors(positions, masses)
    projector = np.eye(mw_hessian.shape[0]) - q @ q.T
    projected = projector @ np.asarray(mw_hessian) @ projector
    return 0.5 * (projected + projected.T), q.shape[1]
