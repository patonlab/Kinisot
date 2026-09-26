"""Re-optimize the Claisen reactant and transition structure with GFN2-xTB and compute their Hessians.

    pip install "kinisot[ase]" tblite sella
    python scripts/make_xtb_claisen.py

Starts from the B3LYP/6-31G(d) geometries in tests/data/gaussian, minimizes
the reactant (BFGS) and refines the saddle point (Sella, first-order) with
GFN2-xTB through its ASE calculator, then writes to tests/data/xtb/:

    claisen_gs.xyz, claisen_ts.xyz                  optimized geometries
    claisen_gs.hessian.json, claisen_ts.hessian.json Hessians (VibrationsData JSON)

The calculator is exactly what ``--calc xtb`` builds (GFN2-xTB, SCF
accuracy 0.01); the Hessians use central differences with 0.005 Angstrom
and four displacements per coordinate.
"""

import os

import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from sella import Sella

from kinisot import build_calculator, hessian_from_calculator, parse_gaussian, save_hessian_json
from kinisot.thermo import BOHR_TO_ANGSTROM

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SOURCE = os.path.join(ROOT, "tests", "data", "gaussian")
TARGET = os.path.join(ROOT, "tests", "data", "xtb")
FMAX = 1e-4  # eV/Angstrom
DELTA, NFREE = 0.005, 4


def calculator():
    # exactly what `--calc xtb` builds: GFN2-xTB with a tight SCF (accuracy 0.01), because
    # finite-difference Hessians are sensitive to noise in the forces
    return build_calculator("xtb")


def dft_geometry(name):
    data = parse_gaussian(os.path.join(SOURCE, name + ".out"))
    return Atoms(numbers=list(data.atomic_numbers), positions=np.asarray(data.positions) * BOHR_TO_ANGSTROM)


def main():
    os.makedirs(TARGET, exist_ok=True)
    for name, saddle in (("claisen_gs", False), ("claisen_ts", True)):
        atoms = dft_geometry(name)
        atoms.calc = calculator()
        if saddle:
            Sella(atoms, order=1, internal=True, logfile=None).run(fmax=FMAX, steps=500)
        else:
            BFGS(atoms, logfile=None).run(fmax=FMAX, steps=500)
        fmax = np.abs(atoms.get_forces()).max()
        atoms.info["level_of_theory"] = "GFN2-xTB"
        write(os.path.join(TARGET, name + ".xyz"), atoms)
        data = hessian_from_calculator(atoms, calculator(), delta=DELTA, nfree=NFREE, source=name)
        path = save_hessian_json(data, os.path.join(TARGET, name + ".hessian.json"), atoms=atoms, calc_spec="xtb")
        freqs = np.asarray(data.program_frequencies)
        print(
            "%s: max force %.1e eV/A, energy %.6f Eh, %d imaginary, lowest %s -> %s"
            % (name, fmax, data.energy, int((freqs < -50).sum()), np.round(np.sort(freqs)[:3], 1), path)
        )


if __name__ == "__main__":
    main()
