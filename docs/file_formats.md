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

## ORCA (supported)

ORCA writes the Hessian to `name.hess` next to `name.out`; give Kinisot
either path and it finds the other by name (the `.hess` file alone is
enough for the numbers; the `.out` supplies the level of theory from the
`!` keyword line, parsed by GoodVibes). From the `.hess` file Kinisot reads
`$hessian` (through GoodVibes), `$atoms` (element symbols and the geometry
in the Hessian's frame, used to decide linearity) and
`$vibrational_frequencies` (for the self-check below).

**Masses.** ORCA lists standard atomic weights (C 12.011, H 1.008). For the
elements Kinisot can substitute, the light isotopologue is built from the
pure most-abundant-isotope masses Gaussian uses (¹²C 12.000, ¹H 1.00783,
¹⁶O 15.99491), so the same Hessian gives the same isotope effect from
either program; other elements keep ORCA's masses until the full isotope
table of Phase 7. This means Kinisot's unsubstituted frequencies differ
very slightly from ORCA's printed ones (well under the 1 cm⁻¹ self-check
tolerance for organic molecules).

Caveats: `NumFreq` Hessians are noisier than analytic ones (projection,
Phase 7, is recommended for them); `%freq scalfreq` scales ORCA's printed
frequencies but not `$hessian`, so Kinisot's factor is applied to the raw
frequencies as intended.

## Frequency self-check

For every unsubstituted species whose file lists the program's vibrational
frequencies, Kinisot compares them with the ones it obtains from the
Hessian and warns when they differ by more than 1 cm⁻¹ (the program
projects out translations and rotations, Kinisot does not, hence the
tolerance). The warning means the Hessian and the masses do not belong
together (wrong file pairing, edited output, unit problem).

## ASE and machine-learned potentials (supported, `pip install kinisot[ase]`)

Two forms of input:

1. **A `VibrationsData` JSON file**: `ase.vibrations.VibrationsData.todict()`
   encoded with `ase.io.jsonio` (geometry in Å, Hessian in eV/Å²).
   `kinisot.save_hessian_json()` writes it (from any HessianInput, e.g. to
   convert a Gaussian Hessian) and `--calc` writes it as its cache. The
   electronic energy (eV) may be stored under `atoms.info["energy"]` for the
   Skodje–Truhlar barrier, and `atoms.info["level_of_theory"]` for scaling.
2. **A geometry file plus `--calc`** (or `calculator=` in the API): anything
   `ase.io.read` accepts. The Hessian is computed by central finite
   differences (`--delta`, default 0.01 Å) or with the calculator's
   `get_hessian` when it has one, and cached as `<name>.hessian.json` next
   to the geometry (reused while the geometry is unchanged and the `--calc`
   specification matches). The geometry must be a stationary point of that
   calculator.

Masses come from Kinisot's isotope table (ASE's standard atomic weights
only identify elements). Projection of the external modes is on by default
for these inputs, no scaling factor is applied unless `-s` is given, and
imaginary modes below the cutoff are treated as real vibrations of the same
magnitude with a warning. Calculator names: `emt` (ASE's test potential),
`xtb[:method]` (GFN2-xTB through `tblite`, SCF accuracy tightened to 0.01),
`mace_mp[:model]`, `mace_off[:model]`, `mace_omol`, `orb`, `sevennet`,
`aimnet2`, or `module.path:callable`; each package must be installed
separately.

**Finite-difference Hessians need a tightly converged energy.** Noise in
the forces becomes noise in the Hessian and then in the isotope effect:
with tblite's default SCF accuracy, an EQE between water's two equivalent
hydrogens (exactly 1 by symmetry) comes out 1.0001 or 0.9996; with
`accuracy=0.01` it is 1 ± 3 × 10⁻⁷. Converge any calculator well beyond
what a geometry optimization needs, and use a symmetric EQE of this kind as
a quick check of a new setup.

## Frequency self-check

For every unsubstituted species whose file lists the program's vibrational
frequencies, Kinisot compares them with the ones it obtains from the
Hessian and warns when they differ by more than 1 cm⁻¹ (the program
projects out translations and rotations, Kinisot does not, hence the
tolerance). The warning means the Hessian and the masses do not belong
together (wrong file pairing, edited output, unit problem).

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
