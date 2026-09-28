"""Conformers of allyl vinyl ether and the chair and boat Claisen transition structures at GFN2-xTB.

    pip install "kinisot[ase]" sella tblite
    python scripts/make_conformer_example.py --out examples/conformers

Makes the structures of the conformer-ensemble worked example
(examples/conformers/README.md), in about a minute:

1. **Reactant conformers.** Starting from the GFN2-xTB reactant in
   tests/data/xtb, every combination of the three rotatable dihedrals
   (C1=C2-O3-C4 at 0 and 180 degrees, C2-O3-C4-C5 at 60, 180 and 300,
   O3-C4-C5=C6 at 0, 120 and 240) is minimized (BFGS, after a small seeded
   displacement so that no start is held at a symmetric saddle point).
   Structures that are the same, or mirror images of each other, up to
   renumbering of hydrogens on one carbon, are kept once. Conformers within
   ``--window`` kcal/mol (electronic energy) of the lowest are kept.
2. **Transition structures.** The chair is the GFN2-xTB saddle point in
   tests/data/xtb. The boat starts from it with the allyl fragment
   (C4, C5, C6 and their hydrogens) reflected through the plane of C1, O3,
   C4 and C6, which moves C5 to the same side as C2, and is refined with
   Sella.
3. **Checks.** Every Hessian (central differences, 0.005 Angstrom, four
   displacements, as ``kinisot --calc`` does) must show a minimum, or
   exactly one imaginary mode. Each saddle point must connect allyl vinyl
   ether with 4-pentenal: displaced a little along its imaginary mode either
   way and minimized, it must give the reactant's bonds on one side and the
   product's (C1-C6 made, C4-O3 broken) on the other. The chair and boat
   must be what their names say (C2 and C5 on opposite or on the same side
   of the plane of the four terminal atoms). The boat is not required to
   pass the stricter test of scripts/make_claisen_structures.py, that the
   C1-C6 and C4-O3 stretches dominate its imaginary mode: at GFN2-xTB it is
   asynchronous (C4-O3 1.51 Angstrom, against 1.59 in the chair), and C1-C6
   formation leads.
4. **Degeneracies.** A conformer without a mirror plane has a mirror image
   with the same energy that the list does not repeat, so its degeneracy is
   2; one that is its own mirror image (up to renumbering hydrogens on one
   carbon) has 1. They go into ``degeneracies.txt`` (the ``--energies``
   format, free energies left to Kinisot).

Every structure keeps the atom numbering of tests/data/xtb, so an isotope
label means the same atom in all of them. Written to ``--out``:
gs_<n>.xyz and gs_<n>.hessian.json (n in order of electronic energy),
ts_chair.* and ts_boat.*, and degeneracies.txt.
"""

import argparse
import itertools
import os
import sys

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS

from kinisot import build_calculator, hessian_from_calculator, save_hessian_json
from kinisot.thermo import HARTREE_TO_KCAL_PER_MOL

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_claisen_structures import DELTA, FMAX, IMAGINARY, NFREE, check_structure, normal_modes  # noqa: E402

ROOT = os.path.join(HERE, "..")
START = os.path.join(ROOT, "tests", "data", "xtb")
LABEL = "GFN2-xTB"
EV_TO_KCAL = 23.060547830619
# 0-based atom indices: C1 C2 O3 C4 C5 C6 are 0-5
DIHEDRALS = {
    (0, 1, 2, 3): (0.0, 180.0),  # C1=C2-O3-C4: vinyl ether s-cis or s-trans
    (1, 2, 3, 4): (60.0, 180.0, 300.0),  # C2-O3-C4-C5
    (2, 3, 4, 5): (0.0, 120.0, 240.0),  # O3-C4-C5=C6
}
TERMINAL = (0, 2, 3, 5)  # C1, O3, C4, C6
RADII = {1: 0.31, 6: 0.76, 8: 0.66}


def bonds(atoms):
    x, z = atoms.positions, atoms.numbers
    return {
        (i, j)
        for i, j in itertools.combinations(range(len(atoms)), 2)
        if np.linalg.norm(x[i] - x[j]) < 1.2 * (RADII[z[i]] + RADII[z[j]])
    }


def moving_side(bond_set, n, a, b):
    """Atoms on b's side of the a-b bond (the part a dihedral rotation about a-b moves)."""
    neighbours = {i: set() for i in range(n)}
    for i, j in bond_set:
        if {i, j} != {a, b}:
            neighbours[i].add(j)
            neighbours[j].add(i)
    seen, stack = {b}, [b]
    while stack:
        for k in neighbours[stack.pop()]:
            if k not in seen:
                seen.add(k)
                stack.append(k)
    return [i in seen for i in range(n)]


def hydrogen_groups(atoms):
    """Hydrogens bonded to the same carbon: the renumberings that leave a structure the same."""
    groups = {}
    for i, j in bonds(atoms):
        h, c = (i, j) if atoms.numbers[i] == 1 else (j, i)
        if atoms.numbers[h] == 1 and atoms.numbers[c] != 1:
            groups.setdefault(c, []).append(h)
    return [sorted(g) for g in groups.values() if len(g) > 1]


def rmsd_after_alignment(x, y):
    """Kabsch RMSD of y onto x (rotation only, no reflection)."""
    x, y = x - x.mean(0), y - y.mean(0)
    u, _, vt = np.linalg.svd(y.T @ x)
    d = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(((y @ rotation - x) ** 2).sum(1).mean()))


def same_structure(a, b, groups, mirror=False, tolerance=0.05):
    """True when b (or its mirror image) superimposes on a after renumbering hydrogens on one carbon."""
    y = b.positions * (np.array([1.0, 1.0, -1.0]) if mirror else 1.0)
    for choice in itertools.product(*[list(itertools.permutations(g)) for g in groups]):
        order = list(range(len(a)))
        for group, permuted in zip(groups, choice):
            for old, new in zip(group, permuted):
                order[old] = new
        if rmsd_after_alignment(a.positions, y[order]) < tolerance:
            return True
    return False


def degeneracy(atoms, groups):
    """1 for a structure that is its own mirror image, 2 for a chiral one."""
    return 1 if same_structure(atoms, atoms, groups, mirror=True) else 2


def side_of_plane(atoms):
    """Signed distances of C2 and C5 from the plane of C1, O3, C4 and C6."""
    points = atoms.positions[list(TERMINAL)]
    centre = points.mean(0)
    normal = np.linalg.svd(points - centre)[2][2]
    return float((atoms.positions[1] - centre) @ normal), float((atoms.positions[4] - centre) @ normal)


def heavy_bonds(atoms):
    return {(i, j) for i, j in bonds(atoms) if atoms.numbers[i] > 1 and atoms.numbers[j] > 1}


def connects(atoms, mode, calculator, reactant, product):
    """Minima reached from the saddle point along +mode and -mode: reactant on one side, product on the other."""
    ends = []
    for sign in (1.0, -1.0):
        moved = atoms.copy()
        moved.positions += sign * 0.15 * mode
        moved.calc = calculator()
        BFGS(moved, logfile=None).run(fmax=1e-3, steps=3000)
        ends.append(heavy_bonds(moved))
    return sorted(map(sorted, ends)) == sorted(map(sorted, [reactant, product]))


def minimize(atoms, calculator, seed):
    atoms = atoms.copy()
    atoms.positions += np.random.default_rng(seed).normal(scale=0.02, size=atoms.positions.shape)
    atoms.calc = calculator()
    BFGS(atoms, logfile=None).run(fmax=FMAX, steps=2000)
    return atoms


def hessian(atoms, calculator, name):
    atoms.info["level_of_theory"] = LABEL
    return hessian_from_calculator(atoms, calculator(), delta=DELTA, nfree=NFREE, source=name)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--window", type=float, default=3.0, help="energy window for reactant conformers (kcal/mol)")
    options = parser.parse_args(argv)
    from sella import Sella

    def calculator():
        return build_calculator("xtb")

    start = read(os.path.join(START, "claisen_gs.xyz"))
    bond_set = bonds(start)
    groups = hydrogen_groups(start)
    n = len(start)

    # 1. reactant conformers
    found = []
    for seed, angles in enumerate(itertools.product(*DIHEDRALS.values())):
        guess = start.copy()
        for (a, b, c, d), angle in zip(DIHEDRALS, angles):
            guess.set_dihedral(a, b, c, d, angle, mask=moving_side(bond_set, n, b, c))
        atoms = minimize(guess, calculator, seed)
        if bonds(atoms) != bond_set:
            print("start %s: bonding changed on minimization, skipped" % (angles,))
            continue
        energy = atoms.get_potential_energy() * EV_TO_KCAL
        if any(
            abs(energy - e) < 0.05 and (same_structure(a, atoms, groups) or same_structure(a, atoms, groups, True))
            for a, e in found
        ):
            continue
        found.append((atoms, energy))
    found.sort(key=lambda item: item[1])
    lowest = found[0][1]
    kept = [(atoms, energy) for atoms, energy in found if energy - lowest <= options.window]
    print("%d distinct reactant conformers, %d within %.1f kcal/mol" % (len(found), len(kept), options.window))

    written = []
    for k, (atoms, energy) in enumerate(kept, 1):
        name = "gs_%d" % k
        data = hessian(atoms, calculator, name)
        problems = check_structure(data, saddle=False)
        if problems:
            # a start held near a saddle point: follow the imaginary mode down and minimize again
            _, mode = normal_modes(data)
            atoms.positions += 0.1 * mode
            atoms = minimize(atoms, calculator, 1000 + k)
            data = hessian(atoms, calculator, name)
            problems = check_structure(data, saddle=False)
        if problems:
            print("%s rejected: %s" % (name, "; ".join(problems)))
            return 1
        dihedrals = [atoms.get_dihedral(*d) for d in DIHEDRALS]
        g = degeneracy(atoms, groups)
        print(
            "%s: E %+.2f kcal/mol, dihedrals %s, degeneracy %d"
            % (name, energy - lowest, " ".join("%.0f" % v for v in dihedrals), g)
        )
        written.append((name, atoms, data, g))

    # 2. transition structures
    chair = read(os.path.join(START, "claisen_ts.xyz"))
    boat = chair.copy()
    points = boat.positions[list(TERMINAL)]
    centre = points.mean(0)
    normal = np.linalg.svd(points - centre)[2][2]
    allyl = moving_side(bond_set, n, 2, 3)  # C4, C5, C6 and their hydrogens (the O3-C4 bond is the breaking one)
    for i in np.flatnonzero(allyl):
        boat.positions[i] -= 2 * ((boat.positions[i] - centre) @ normal) * normal
    reactant = heavy_bonds(start)
    product = (reactant - {(2, 3)}) | {(0, 5)}  # C4-O3 broken, C1-C6 made
    for name, atoms, same_side in (("ts_chair", chair, False), ("ts_boat", boat, True)):
        atoms.calc = calculator()
        Sella(atoms, order=1, internal=True, logfile=None).run(fmax=FMAX, steps=2000)
        data = hessian(atoms, calculator, name)
        frequencies, mode = normal_modes(data)
        imaginary = frequencies[frequencies < -IMAGINARY]
        problems = [] if len(imaginary) == 1 else ["%d imaginary modes (one expected)" % len(imaginary)]
        if not problems and not connects(atoms, mode, calculator, reactant, product):
            problems.append("it does not connect allyl vinyl ether with 4-pentenal")
        c2, c5 = side_of_plane(atoms)
        if (c2 * c5 > 0) != same_side:
            problems.append("C2 and C5 are on %s sides: not a %s" % ("the same" if c2 * c5 > 0 else "opposite", name))
        print(
            "%s: E %+.2f kcal/mol above gs_1, imaginary %.1fi cm-1, C1-C6 %.3f, C4-O3 %.3f A, degeneracy %d"
            % (name, (data.energy * HARTREE_TO_KCAL_PER_MOL) - written[0][2].energy * HARTREE_TO_KCAL_PER_MOL,
               -frequencies[0], atoms.get_distance(0, 5), atoms.get_distance(2, 3), degeneracy(atoms, groups))
        )  # fmt: skip
        if problems:
            print("%s rejected: %s" % (name, "; ".join(problems)))
            return 1
        written.append((name, atoms, data, degeneracy(atoms, groups)))

    os.makedirs(options.out, exist_ok=True)
    for name, atoms, data, _ in written:
        write(os.path.join(options.out, name + ".xyz"), atoms)
        save_hessian_json(data, os.path.join(options.out, name + ".hessian.json"), atoms=atoms, calc_spec="xtb")
    with open(os.path.join(options.out, "degeneracies.txt"), "w") as handle:
        handle.write("# file, free energy (- : computed by Kinisot), degeneracy (2: its mirror image is not listed)\n")
        for name, _, _, g in written:
            handle.write("%s.hessian.json - %d\n" % (name, g))
    print("wrote %d structures to %s" % (len(written), options.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
