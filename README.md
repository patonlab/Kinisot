![Kinisot Banner](https://github.com/patonlab/Kinisot/blob/master/kinisot_banner.png)

[![DOI](https://zenodo.org/badge/54840251.svg)](https://zenodo.org/badge/latestdoi/54840251)
[![PyPI version](https://badge.fury.io/py/kinisot.svg)](https://badge.fury.io/py/kinisot)
[![CI](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml/badge.svg)](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml)

**Kinisot** computes kinetic (KIE) and equilibrium (EQE) isotope effects
from Gaussian or ORCA frequency calculations, or from Hessians computed with
any ASE calculator, machine-learned potentials included. Give it the output files of a
reactant and a transition structure (or a product), say which atoms carry
the heavy isotope, and it re-mass-weights the Hessians, diagonalizes them,
and evaluates the Bigeleisen–Mayer equation with a Bell tunnelling
correction at any temperature. No new frequency job is needed for each
isotopologue. Kinisot is developed in the
[Paton group](https://patonlab.colostate.edu) at Colorado State University
and is a rewrite of the Fortran Kinisot by
[Henry Rzepa](https://en.wikipedia.org/wiki/Henry_Rzepa).

## Quick start

```
pip install kinisot            # or: conda install -c conda-forge kinisot
git clone https://github.com/patonlab/Kinisot && cd Kinisot/tests/data/gaussian
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 1 -t 393 -s 0.961
```

This computes the ¹³C KIE at carbon 1 of allyl vinyl ether for its Claisen
rearrangement at 393 K, with B3LYP/6-31G(d) frequencies scaled by 0.961:

```
  KINISOT.py v 2.5.0: 2026-09-25 12:45
  Species: claisen_gs.out isotopologue: 1
  Species: claisen_ts.out isotopologue: 1

                                                     Temp = 393.0K / Vib. scale factor = 0.961
                                                     V-ratio        ZPE        EXC       TRPF        KIE    1D-tunn   corr-KIE

o claisen_gs                                        --------------------------------------------------------------------------
o claisen_ts                                          463.9
o claisen_gs: iso @ 1                                        1.148e+00  1.024e+00  9.246e-01
o claisen_ts: iso @ 1                                 460.3  1.144e+00  1.021e+00  9.269e-01
                                                    --------------------------------------------------------------------------
  KIE @ 393.0 K                                    1.007880   1.003953   1.003421   0.997456   1.012744   1.001970   1.014739
                                                    --------------------------------------------------------------------------

  Vibrational modes (scaled, cm-1): kept in the partition function / discarded as external modes
  claisen_gs (light): 36 kept; discarded: -0.2 -0.0 0.1 5.4 16.5 23.1
  claisen_gs (iso @ 1): 36 kept; discarded: -0.2 -0.0 0.1 5.4 16.3 23.0
  claisen_ts (light): imaginary 463.9i; 35 kept; discarded: -9.7 -0.2 0.1 0.1 13.3 17.0
  claisen_ts (iso @ 1): imaginary 460.3i; 35 kept; discarded: -9.6 -0.2 0.1 0.1 13.2 16.9
```

The number to report is **corr-KIE = 1.015** (the semiclassical KIE of
1.013 times the tunnelling correction). The same block is appended to
`Kinisot_output.dat`.

## Reading the output

- **Species rows** show, for the reactant and the transition structure,
  the imaginary frequency of each isotopologue (cm⁻¹) and the three
  Bigeleisen–Mayer factors of that species as light/heavy ratios:
  the zero-point energy term (`ZPE`), the excitation term (`EXC`) and the
  Teller–Redlich product of frequencies (`TRPF`, heavy/light). Their
  product is the reduced isotopic partition function ratio (s/s′)f.
- **`KIE @ T`** is the ratio reactant/transition structure of every column:
  `V-ratio` = ν‡(light)/ν‡(heavy); `ZPE`, `EXC`, `TRPF` as above;
  `KIE` = V-ratio × ZPE × EXC × TRPF (semiclassical); `1D-tunn` = Bell
  infinite-parabola tunnelling correction; `corr-KIE` = KIE × 1D-tunn.
- **Vibrational modes**: how many modes entered each partition function
  and which six (five for a linear molecule, plus the reaction coordinate)
  were discarded as translations and rotations. On a converged geometry
  the discarded values are within a few tens of cm⁻¹ of zero; a genuine
  vibration in that list means the geometry needs tightening.
- For an EQE the line is labelled `EQE @ T`, `V-ratio` is empty and
  `1D-tunn` is 1.

The equations behind every column are in [docs/theory.md](docs/theory.md).

## Usage

```
kinisot --rct FILE [--rct FILE ...] (--ts FILE | --prd FILE) --iso ATOMS [--iso ATOMS ...]
        [-t K] [-s FACTOR] [--imag-cutoff CM-1] [--tunneling MODEL] [--barrier KCAL] [--project]
        [--reference ATOMS] [--calc SPEC] [-o FILE] [--overwrite] [-q] [--json FILE] [--csv FILE]
```

`python -m kinisot` is equivalent to `kinisot`.

**Three ways to describe a reaction**

| Case | Command | Labels |
| --- | --- | --- |
| KIE, one reactant | `kinisot --rct gs.out --ts ts.out --iso 5` | one `--iso` when the atom numbering is the same in both files |
| KIE, bimolecular | `kinisot --rct dienophile.out --rct diene.out --ts ts.out --iso 0 --iso 6 --iso 15` | one `--iso` per file, `--rct` files first, `0` for a file without a substituted atom |
| EQE | `kinisot --rct conf_a.out --prd conf_b.out --iso 24,25,26 --iso 28,29,30` | the same, with the product instead of a TS |

**Options**

| Flag | Meaning |
| --- | --- |
| `--iso ATOMS` | atom number(s) to replace with the heavy isotope, comma separated (`7,8` substitutes two hydrogens). Atom numbers follow the order of the atoms in the quantum-chemistry input. Kinisot checks that the atoms exist, can be substituted, still carry the light isotope, and that both sides of the reaction substitute the same elements. |
| `-t`, `--temperature` | temperature in K at which the partition functions are evaluated (default 298.15). It need not match the frequency job. A list (`273,298,323`) or range (`250:350:10`) gives one result line per temperature. |
| `-s`, `--scale` | vibrational scaling factor. Default: the ZPE factor of the [Truhlar database](https://comp.chem.umn.edu/freqscale/) (version 5, via GoodVibes) for the level of theory detected in the files, or 1.0 with a message if it is not listed. `--scale-type harm` or `fund` picks the harmonic or fundamental factor instead. |
| `--imag-cutoff` | a mode below −CUTOFF cm⁻¹ is the reaction coordinate (default 50). Reactants and products must have none. |
| `--tunneling` | `bell` (default), `wigner`, `skodje` (Skodje–Truhlar, needs the barrier: `--barrier KCAL` or the electronic energies in the files) or `none`. |
| `--project` | project translations and rotations out of the Hessian (Eckart) instead of discarding the six lowest modes. On by default for `--calc`/ASE inputs, off for Gaussian and ORCA (`--no-project` forces it off). |
| `--calc SPEC`, `--delta` | compute the Hessians of geometry inputs with an ASE calculator (`mace_mp:medium`, `emt`, `module:callable`, ...); finite-difference step in Å. |
| `--reference ATOMS` | a second isotopologue; the KIE is also reported divided by its KIE (natural-abundance NMR style). |
| `-o`, `--output` | results file (default `Kinisot_output.dat`); results are appended. `--overwrite` starts afresh, `-q` keeps the terminal quiet. |
| `--json FILE`, `--csv FILE` | also write the full result as JSON, or append one summary row to a CSV file. |
| `--version` | print the version. |

Invalid input (an atom number out of range, a reactant with an imaginary
frequency, a file that is not a completed frequency job, ...) stops the run
with a message and exit code 1.

**Isotopes.** A bare atom number substitutes the usual heavy label:
¹H → ²H, ¹²C → ¹³C, ¹⁴N → ¹⁵N, ¹⁶O → ¹⁸O, ³²S → ³⁴S, ³⁵Cl → ³⁷Cl,
⁷⁹Br → ⁸¹Br, ²⁸Si → ²⁹Si. Any isotope in the table can be asked for
explicitly: `--iso 3:17O`, `--iso 7:D` (or `7:T`), `--iso 5:14C`, or an
explicit mass `--iso 5:13.5`. Masses are AME 2020 values for every
naturally occurring isotope of 83 elements plus the common radioactive
labels; both isotopologues are built from this table (the light one from
the most abundant isotope, as Gaussian does), so Gaussian, ORCA and other
programs give the same numbers for the same Hessian. Note that Kinisot
≤ 2.3 substituted ¹⁷O for a bare oxygen index; 2.4 warns once when it sees
one.

## Input files

| Program | Status | What Kinisot reads |
| --- | --- | --- |
| Gaussian 09/16 | supported | a normally terminated `freq` job: atom masses and the force constants in the archive entry (`opt freq` jobs are fine) |
| ORCA 5/6 | supported | `name.out` plus the `name.hess` file ORCA writes next to it (give either path); level of theory from the `!` line |
| ASE / machine-learned potentials | supported (`pip install kinisot[ase]`) | a `VibrationsData` JSON file, or a geometry plus `--calc` (MACE, ORB, SevenNet, AIMNet2, any ASE calculator); the Hessian is computed and cached next to the geometry, external modes are projected out |

Details and pitfalls: [docs/file_formats.md](docs/file_formats.md).

## What Kinisot assumes

- Harmonic frequencies from the program's Hessian; the scaling factor is
  applied to all modes, including the imaginary one. Kinisot checks that it
  reproduces the frequencies the program printed and warns if not.
- Translations and rotations are removed by discarding the lowest 5/6
  eigenvalues; `--project` removes them by Eckart projection instead. On
  converged Gaussian geometries the two agree to better than 4 × 10⁻⁷ in
  the KIE; projection matters for noisy (finite-difference) Hessians.
- One imaginary mode per transition structure; a second one triggers a
  warning.
- Tunnelling by Bell's one-dimensional infinite-parabola model (default), the
  Wigner correction, the Skodje–Truhlar correction (which also uses the
  barrier height), or none; Bell is refused below the crossover
  temperature h c |ν‡| / 2π k.
- Rigid rotor / harmonic oscillator, no conformational averaging, no
  solvent corrections: compute those outside Kinisot if you need them.

## Python API

```python
from kinisot import compute_kie, parse_gaussian

r = compute_kie(rct="claisen_gs.out", ts="claisen_ts.out", iso="1", temperature=393.0, scale=0.961)
r.kie_tunnel                          # 1.014739 (the corr-KIE column)
r.kie, r.zpe, r.exc, r.trpf, r.imag_ratio, r.tunnel_corr
r.other.light.imaginary               # 463.9  (light TS reaction coordinate, cm-1)
r.reactant.rpfr, r.other.rpfr         # reduced partition function ratios (s/s')f
r.to_dict(); r.to_json()              # everything, JSON-serializable
r.summary_row()                       # one flat row, as written by --csv

gs, ts = parse_gaussian("claisen_gs.out"), parse_gaussian("claisen_ts.out")   # parse once,
[compute_kie(rct=gs, ts=ts, iso=a, temperature=393.0, scale=0.961).kie_tunnel  # scan positions
 for a in ("1", "2", "3", "4", "5", "6")]
```

`compute_kie(rct, ts=None, prd=None, iso=None, temperature=298.15, scale=1.0,
imag_cutoff=50.0, tunneling="bell", scale_type="zpe", project=None,
barrier=None, reference=None, calculator=None)` takes Gaussian or ORCA file paths or `HessianInput`
objects (one or a list per side), `scale=None` for automatic Truhlar
lookup, and `tunneling` in `"bell"`, `"wigner"`, `"none"`. It returns a
frozen `IsotopeEffect` whose `reactant` and `other` (transition structure
or product) sides hold the light and heavy isotopologues with their kept,
discarded and imaginary frequencies, masses and substitutions. Problems
raise `kinisot.KinisotInputError` / `KinisotParseError` (both `ValueError`
subclasses); non-fatal ones are `KinisotWarning`s and are listed in
`r.warnings`. See [examples/api_example.py](examples/api_example.py) and
[examples/examples.ipynb](examples/examples.ipynb). The 2.x function
`compute_isotope_effect()` still works but is deprecated.

Machine-readable output from the command line: `--json run.json` (the
full result of one run) and `--csv runs.csv` (one row per run, appended).

## Examples and documentation

- [examples/](examples/README.md): the Claisen rearrangement (¹³C, ¹⁸O/¹⁷O
  and ²H KIEs, projection, reference isotopologue, tunnelling models,
  temperature scans), a Diels–Alder reaction with one or two reactant files,
  the same KIE from ORCA and from ASE inputs (with the machine-learned
  potential workflow), and a conformational EQE, each with the commands, the
  expected numbers and the literature background.
- [docs/theory.md](docs/theory.md), [docs/file_formats.md](docs/file_formats.md),
  [docs/faq.md](docs/faq.md), [docs/comparison.md](docs/comparison.md)
  (PyQuiver, PyQuiverHS, Gaussian's `readisotopes`, GoodVibes).
- [benchmarks/](benchmarks/README.md): computed versus experimental KIEs,
  one directory per reaction, with a runner that regenerates the report.
- [CHANGELOG.md](CHANGELOG.md) and [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md)
  for what changed and what is coming.
- A [video guide](http://www.youtube.com/watch?v=r4x2gmkc0U8) to an older
  version, and Henry Rzepa's blog post on
  [computing KIE values](http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327).

## Citing Kinisot

Please cite the Zenodo record, DOI
[10.5281/zenodo.19272](http://dx.doi.org/10.5281/zenodo.19272)
([CITATION.cff](CITATION.cff) has the metadata; GitHub's "Cite this
repository" button formats it). Vibrational scaling factors come from
I. M. Alecu, J. Zheng, Y. Zhao, D. G. Truhlar, *J. Chem. Theory Comput.*
**2010**, *6*, 2872.

## Support and contributing

Questions and bug reports: the [issue tracker](https://github.com/patonlab/Kinisot/issues)
or [patonlab@colostate.edu](mailto:patonlab@colostate.edu). See
[SUPPORT.md](SUPPORT.md) for what to include and
[CONTRIBUTING.md](CONTRIBUTING.md) for how to run the tests and add a
backend. MIT licensed.
