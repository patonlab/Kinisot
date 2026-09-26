"""Re-optimize the Claisen reactant and transition structure with any ASE calculator and compute their Hessians.

    pip install "kinisot[ase]" sella            # plus the calculator's package
    python scripts/make_claisen_structures.py --calc xtb --out tests/data/xtb
    python scripts/make_claisen_structures.py --calc mace_off:medium --out tests/data/mace_off23_medium

Starts from the B3LYP/6-31G(d) geometries in tests/data/gaussian, minimizes
the reactant (BFGS) and refines the saddle point (Sella, first order, internal
coordinates) with the calculator that ``--calc`` builds in Kinisot, then writes
to ``--out``:

    claisen_gs.xyz, claisen_ts.xyz                     optimized geometries
    claisen_gs.hessian.json, claisen_ts.hessian.json   Hessians (VibrationsData JSON)

Hessians use the calculator's analytic ``get_hessian`` when it has one
(MACE) and otherwise central differences (0.005 Angstrom, four displacements
per coordinate), exactly as ``kinisot --calc`` does. ``--model-file`` loads
local weights while recording the portable ``--calc`` specification in the
output (for machines that cannot download the weights).

Every structure is checked before it is written: the reactant must be a
minimum, and the saddle point must be the Claisen transition structure (one
imaginary mode, C1-C6 and C4-O3 both partial bonds, and those two stretches
dominating the imaginary mode). A saddle search can otherwise wander onto a
different reaction (C-O dissociation, ring closure) and yield plausible
looking but meaningless isotope effects; see examples/mlip_claisen/README.md.
If either structure fails, the script exits with status 1 and writes
neither, so an old file cannot end up paired with a new one
(``--keep-invalid`` writes both anyway, for inspection).
"""

import argparse
import itertools
import os
import sys

import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS

from kinisot import build_calculator, hessian_from_calculator, parse_gaussian, save_hessian_json
from kinisot.hessian import mass_weight
from kinisot.isotopes import light_masses
from kinisot.projection import project_external_modes
from kinisot.thermo import BOHR_TO_ANGSTROM, HESSIAN_TO_WAVENUMBER_SQ

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SOURCE = os.path.join(ROOT, "tests", "data", "gaussian")
FMAX = 1e-4  # eV/Angstrom
DELTA, NFREE = 0.005, 4
# Modes below -IMAGINARY cm-1 count as imaginary; after projection the noise at converged geometries is far smaller
IMAGINARY = 20.0

HEAVY = {1: "C1", 2: "C2", 3: "O3", 4: "C4", 5: "C5", 6: "C6"}
# Partial-bond windows, generous around B3LYP (2.31 / 1.90 A) and GFN2-xTB (1.94 / 1.59 A)
FORMING_C1_C6 = (1.75, 2.9)
BREAKING_C4_O3 = (1.5, 2.5)


def normal_modes(data):
    """Projected frequencies (cm-1, ascending) and Cartesian displacement of the lowest mode."""
    masses = np.asarray(light_masses(data))
    projected, n_external = project_external_modes(mass_weight(data.hessian, masses), data.positions, masses)
    eigenvalues, vectors = np.linalg.eigh(projected * HESSIAN_TO_WAVENUMBER_SQ)
    order = np.argsort(np.abs(eigenvalues))[n_external:]
    eigenvalues, vectors = eigenvalues[order], vectors[:, order]
    ascending = np.argsort(eigenvalues)
    frequencies = np.copysign(np.sqrt(np.abs(eigenvalues[ascending])), eigenvalues[ascending])
    lowest = (vectors[:, ascending[0]] / np.repeat(np.sqrt(masses), 3)).reshape(-1, 3)
    return frequencies, lowest / np.linalg.norm(lowest)


def check_structure(data, saddle):
    """Return a list of problems (empty when the structure is what the script is meant to produce)."""
    frequencies, mode = normal_modes(data)
    imaginary = frequencies[frequencies < -IMAGINARY]
    if not saddle:
        return ["reactant has imaginary modes %s" % np.round(imaginary, 1)] if len(imaginary) else []
    problems = []
    if len(imaginary) != 1:
        problems.append("%d imaginary modes %s (one expected)" % (len(imaginary), np.round(imaginary, 1)))
    x = np.asarray(data.positions) * BOHR_TO_ANGSTROM

    def distance(i, j):
        return float(np.linalg.norm(x[i - 1] - x[j - 1]))

    for (i, j), (low, high), what in (
        ((1, 6), FORMING_C1_C6, "forming C1-C6"),
        ((3, 4), BREAKING_C4_O3, "breaking C4-O3"),
    ):
        if not low <= distance(i, j) <= high:
            problems.append(
                "%s is %.2f A, outside the partial-bond window %.2f-%.2f A" % (what, distance(i, j), low, high)
            )

    def stretch(i, j):
        unit = (x[i - 1] - x[j - 1]) / distance(i, j)
        return abs(float(unit @ (mode[i - 1] - mode[j - 1])))

    ranked = sorted(((stretch(i, j), i, j) for i, j in itertools.combinations(HEAVY, 2)), reverse=True)
    top_two = {tuple(sorted((i, j))) for _, i, j in ranked[:2]}
    if len(imaginary) and top_two != {(1, 6), (3, 4)}:
        problems.append(
            "the imaginary mode is dominated by %s, not by the C1-C6 and C4-O3 stretches"
            % ", ".join("%s-%s (%.2f)" % (HEAVY[i], HEAVY[j], value) for value, i, j in ranked[:3])
        )
    return problems


def dft_geometry(name):
    data = parse_gaussian(os.path.join(SOURCE, name + ".out"))
    return Atoms(numbers=list(data.atomic_numbers), positions=np.asarray(data.positions) * BOHR_TO_ANGSTROM)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--calc", required=True, help="Kinisot --calc specification, e.g. xtb or mace_off:medium")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--label", help="level of theory stored in the files (default: from --calc)")
    parser.add_argument("--model-file", help="local model weights to load instead of downloading (MACE)")
    parser.add_argument("--keep-invalid", action="store_true", help="write structures that fail the checks")
    options = parser.parse_args(argv)
    from sella import Sella

    name, _, model = options.calc.partition(":")
    spec = options.calc if not options.model_file else "%s:%s" % (name, options.model_file)
    label = options.label or {"xtb": "GFN2-xTB"}.get(options.calc, options.calc)

    def calculator():
        return build_calculator(spec)

    results = []
    for structure, saddle in (("claisen_gs", False), ("claisen_ts", True)):
        atoms = dft_geometry(structure)
        atoms.calc = calculator()
        if saddle:
            Sella(atoms, order=1, internal=True, logfile=None).run(fmax=FMAX, steps=1000)
        else:
            BFGS(atoms, logfile=None).run(fmax=FMAX, steps=1000)
        fmax = np.abs(atoms.get_forces()).max()
        atoms.info["level_of_theory"] = label
        data = hessian_from_calculator(atoms, calculator(), delta=DELTA, nfree=NFREE, source=structure)
        problems = (["not converged: max force %.1e eV/A" % fmax] if fmax > 10 * FMAX else []) + check_structure(
            data, saddle
        )
        freqs = normal_modes(data)[0]
        summary = "%s %s: max force %.1e eV/A, energy %.6f Eh, %d imaginary, lowest %s cm-1" % (
            label, structure, fmax, data.energy, int((freqs < -IMAGINARY).sum()), np.round(freqs[:3], 1),
        )  # fmt: skip
        print(summary + ("\n  REJECTED: " + "; ".join(problems) if problems else ""))
        results.append((structure, atoms, data, problems))

    failed = any(problems for *_, problems in results)
    if failed and not options.keep_invalid:
        print("nothing written to %s" % options.out)
        return 1
    os.makedirs(options.out, exist_ok=True)
    for structure, atoms, data, _ in results:
        write(os.path.join(options.out, structure + ".xyz"), atoms)
        path = save_hessian_json(
            data, os.path.join(options.out, structure + ".hessian.json"), atoms=atoms, calc_spec=options.calc
        )
        print("wrote " + path)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
