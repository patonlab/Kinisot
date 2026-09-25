# Input files

## Gaussian (supported)

Kinisot needs a **normally terminated frequency job** (`freq` keyword) for
every species. It reads three things from the output:

1. `NAtoms=` — the atom count.
2. `Atom N has atomic number Z and mass M` — one line per atom, printed by
   the frequency job. These masses respect `iso=` or `readisotopes`
   settings in the Gaussian input, which is why Kinisot refuses to
   substitute an atom whose mass is not the light isotope.
3. The **archive entry** at the end of the output (the block starting
   `1\1\GINC-...` and ending with `@`), specifically the level of theory and
   basis set fields and the lower-triangular force-constant matrix that
   follows `NImag=`.

Works: `opt freq` jobs (the last archive entry is the frequency job), `#p`,
`freq=noraman`, `freq=hpmodes`, restricted and unrestricted references,
Windows builds (`|` archive separators), archive entries wrapped across
lines. Does not work: single points, jobs that died before the archive was
written, `freq=readfc`/`geom=checkpoint` jobs whose masses were not printed,
and outputs edited by hand.

Atom numbering is the order of the atoms in the Gaussian input, starting at
1; the `--iso` labels use that numbering.

## ORCA (planned, implementation plan Phase 5)

ORCA writes the Hessian to a separate `name.hess` file next to `name.out`.
Kinisot will accept either file and locate the other. Note that ORCA's
default masses are standard atomic weights (C 12.011) rather than pure
isotopes; Kinisot will build both isotopologues from its own isotope
table so that Gaussian and ORCA inputs give the same numbers for the same
Hessian.

## ASE and machine-learned potentials (planned, Phase 8)

`ase.vibrations.VibrationsData` objects (JSON) and Hessians computed inside
Kinisot with any ASE calculator (MACE, UMA, ORB, SevenNet, AIMNet2, ...).
Finite-difference Hessians leave sizeable translation/rotation residuals,
so this backend will project them out before removing external modes.

## Results file

`Kinisot_output.dat` (or `--output FILE`) collects one block per run:
header with version and time, the species and their isotope labels, the
scaling factor found, the results table, and the modes kept/discarded for
every species. New runs are **appended**; `--overwrite` starts afresh.

`--json FILE` writes the full result of one run (the dictionary returned by
`IsotopeEffect.to_dict()`: every factor, frequency, mass and substitution)
and `--csv FILE` appends one summary row per run with a header line when the
file is new.
