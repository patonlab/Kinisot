![Kinisot Banner](https://github.com/patonlab/Kinisot/blob/master/kinisot_banner.png)

[![DOI](https://zenodo.org/badge/54840251.svg)](https://zenodo.org/badge/latestdoi/54840251)
[![PyPI version](https://badge.fury.io/py/kinisot.svg)](https://badge.fury.io/py/kinisot)
[![CI](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml/badge.svg)](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml)

**Kinisot** predicts kinetic and equilibrium isotope effects from the
frequency calculations you already have. Give it the Gaussian or ORCA
output for a reactant and a transition structure (or a product), say which
atoms carry the heavy isotope, and it returns the KIE, k(light)/k(heavy),
at any temperature. That works for ¹³C, ²H, ¹⁵N, ¹⁷O, ¹⁸O or any other
label, at every position, from one frequency job per structure. Kinisot
evaluates the Bigeleisen–Mayer equation with Bell's tunnelling correction,
the treatment used to test transition structures against natural-abundance
KIE measurements. It also handles conformer ensembles, transition
structures in series, and parallel pathways. It reads Hessians from xTB and
machine-learned potentials too, through ASE.

Kinisot is developed in the [Paton group](https://patonlab.colostate.edu)
at Colorado State University and is a rewrite of the Fortran Kinisot by
[Henry Rzepa](https://en.wikipedia.org/wiki/Henry_Rzepa).

## What you need

- **A frequency calculation on the reactant and on the transition
  structure**, at the same level of theory. That means a Gaussian `freq` (or
  `opt freq`) output, or an ORCA frequency job with its `.hess` file. The
  reactant must be a minimum, with no imaginary frequency. The transition
  structure must have exactly one, and it must describe the bond changes of
  the step you are testing.
- **Each labelled reactant.** For a bimolecular reaction, give each reactant
  (or a reactant complex). A reactant that carries no label can be left out,
  because it cancels.
- **The atom number of each labelled position, in each file.** Kinisot
  numbers atoms from 1, in the order of the input geometry, for every
  program; this is the numbering GaussView shows. Take the numbers from your
  input, not from ORCA's printed tables, which count from 0.
- **The temperature of the experiment.**

For a competition KIE, the transition structure is the one of the first
irreversible step. If no single step is irreversible, for example because
an intermediate can return to the reactant, see
[Several conformers or transition structures](#several-conformers-or-transition-structures).

## Quick start

Kinisot needs Python 3.9 or later. The example files come with the source:

```
pip install kinisot
git clone https://github.com/patonlab/Kinisot && cd Kinisot/tests/data/gaussian
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 1 -t 393 -s 0.961
```

This computes the ¹³C KIE at carbon 1 of allyl vinyl ether for its Claisen
rearrangement at 393 K, with B3LYP/6-31G(d) frequencies scaled by 0.961:

```
  KINISOT.py v 2.6.0: 2026-09-28 04:25
  Species: claisen_gs.out isotopologue: 1
  Species: claisen_ts.out isotopologue: 1

                                                     Temp = 393.0K / Vib. scale factor = 0.961
  Labelled atoms: claisen_gs C1 -> 13C; claisen_ts C1 -> 13C
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
`Kinisot_output.dat`. Its first line gives the Kinisot version, which
`kinisot --version` also prints.

## Reading the output

**The convention.** KIE = k(light)/k(heavy): above 1 is a normal isotope
effect, below 1 an inverse one. For an equilibrium isotope effect (EIE,
labelled `EQE` in the output) of reactant ⇌ product, the value is
K(light)/K(heavy). Above 1, the heavy isotope accumulates in the reactant,
the side with the stiffer vibrations.

**The number to report** is the last column, `corr-KIE`: the KIE with the
tunnelling correction. The other columns show where it comes from.

- **Labelled atoms** names each substituted atom in each file, with its
  element and the isotope it became (`claisen_ts C1 -> 13C`). Check it
  against your structures: an atom number that is off by one but lands on
  another carbon passes every other check. With `--reference`, a
  `Reference atoms` line does the same for the reference position. (This
  line is new in 2.6.)
- **Species rows** show, for the reactant and the transition structure,
  the imaginary frequency of each isotopologue (cm⁻¹) and three factors as
  light/heavy ratios:
  - `ZPE`, the zero-point energy term;
  - `EXC`, the excitation term (thermally populated vibrational levels);
  - `TRPF`, the Teller–Redlich product of frequencies (heavy/light), which
    stands in for the translational and rotational terms.

  Their product is the reduced isotopic partition function ratio (s/s′)f
  of that species.
- **`KIE @ T`** gives, for every column, the reactant's value divided by
  the transition structure's:
  - `V-ratio` = ν‡(light)/ν‡(heavy), the ratio of the imaginary
    frequencies;
  - `ZPE`, `EXC` and `TRPF` as above;
  - `KIE` = V-ratio × ZPE × EXC × TRPF, the semiclassical KIE;
  - `1D-tunn`, Bell's tunnelling correction;
  - `corr-KIE` = KIE × 1D-tunn.
- **Vibrational modes** lists how many modes entered each partition
  function, and the six lowest (five for a linear molecule) that were
  discarded as overall translations and rotations. On a well-converged
  geometry the discarded values are within a few tens of cm⁻¹ of zero. A
  real vibration in that list means the geometry needs tightening.
- For an EIE the line is labelled `EQE @ T`, `V-ratio` is empty and
  `1D-tunn` is 1.

The equations behind every column are in [docs/theory.md](docs/theory.md).

## Comparing with experiment

- **Use the temperature of the experiment** (`-t`). The frequency job does
  not need to have been run at that temperature.
- **Use the same reference.** Natural-abundance NMR measurements give KIEs
  relative to an internal standard assumed to have no isotope effect.
  `--reference ATOMS` computes that position as well and prints a
  `relative to` line. Compare that line with the measured values.
- **Average positions that give one signal.** Examples are the two ortho
  carbons of a phenyl ring, the three hydrogens of a methyl group, and the
  two oxygens of a nitro group. These are distinct in the transition
  structure but give one signal, so compute each placement of the label and
  average them. `kinisot.equivalent_positions()` in Python does this
  exactly. When the positions are equivalent in the reactant, the average
  is the harmonic mean of the separate KIEs, not their arithmetic mean.
- **Keep the tunnelling correction.** A one-dimensional tunnelling
  correction brings predicted heavy-atom KIEs to about the experimental
  uncertainty (Meyer, DelMonte, Singleton, J. Am. Chem. Soc. 1999, 121,
  10865). At the two bond-forming carbons of isoprene in its Diels–Alder
  reaction, Singleton and Thomas measured 1.022(3) and 1.017(2). B3LYP gives
  1.018 and 1.014 without tunnelling, and 1.022 and 1.017 with Bell's
  correction. For primary hydrogen KIEs with a large tunnelling
  contribution, any one-dimensional correction is only approximate.
- **Expect about 0.001 to 0.003 for heavy atoms.** At B3LYP/6-31G(d), with
  scaling and Bell tunnelling, the [benchmarks](benchmarks/README.md)
  reproduce measured ¹³C and ¹⁷O KIEs with mean absolute deviations of
  0.0009 (Claisen rearrangement), 0.0012 (Shi epoxidation) and 0.003
  (Diels–Alder, including ²H).
- **When prediction and experiment disagree, look at the transition
  structure first.** Bigeleisen–Mayer predictions for heavy atoms are
  accurate as long as the calculation has the right mechanism and
  transition-state geometry. Hirschi, Takeya, Hang and Singleton found a
  theory-independent relation between forming-bond distances and ¹³C KIEs
  (J. Am. Chem. Soc. 2009, 131, 2397), so calculations that reproduce the
  measured KIEs share nearly the same transition-state geometry. Then
  check the mechanism: other conformers, a step that is not rate-limiting
  alone, or competing pathways.
  [examples/mlip_claisen](examples/mlip_claisen/README.md) shows an extreme
  case: potentials whose saddle point is a different reaction.

## Usage

**Common cases**

| Case | Command | Labels |
| --- | --- | --- |
| KIE, one reactant | `kinisot --rct gs.out --ts ts.out --iso 5` | one `--iso` when the atom numbering is the same in both files |
| KIE, bimolecular | `kinisot --rct dienophile.out --rct diene.out --ts ts.out --iso 0 --iso 6 --iso 15` | one `--iso` per file, `--rct` files first, `0` for a file without a labelled atom |
| EIE | `kinisot --rct conf_a.out --prd conf_b.out --iso 24,25,26 --iso 28,29,30` | the same, with the product instead of a transition structure |
| Relative to a standard | `kinisot --rct gs.out --ts ts.out --iso 4 --reference 5` | the KIE of position 4 divided by that of position 5 |
| Conformers | `kinisot --rct gs_1.out gs_2.out --ts ts_*.out --iso 5` | several files after one flag are conformers of one species, with the same atom numbering |

**Everyday options**

| Flag | Meaning |
| --- | --- |
| `--iso ATOMS` | atom number(s) to label, comma separated (`7,8` labels two hydrogens). Numbers follow the order of the atoms in the input geometry, from 1. Kinisot checks that the atoms exist and that both sides of the reaction are labelled with the same isotopes. |
| `-t`, `--temperature` | temperature in K (default 298.15). A list (`273,298,323`) or a range (`250:350:10`) gives one result line per temperature. |
| `--reference ATOMS` | a reference position; the KIE is also reported divided by its KIE, as natural-abundance measurements are. |
| `-s`, `--scale` | vibrational scaling factor. By default Kinisot looks up the ZPE factor for the level of theory it detects in the files ([Truhlar database](https://comp.chem.umn.edu/freqscale/), version 5), or uses 1.0 with a message if the level is not listed. |
| `--tunneling` | `bell` (default), `wigner`, `skodje` (Skodje–Truhlar, which also uses the barrier height: `--barrier KCAL`, or the electronic energies in the files) or `none`. |
| `-o`, `--output` | results file (default `Kinisot_output.dat`); results are appended. `--overwrite` starts afresh, `-q` keeps the terminal quiet. |

**More options**

| Flag | Meaning |
| --- | --- |
| `--scale-type` | `harm` or `fund` picks the harmonic or fundamental scaling factor instead of the ZPE one. |
| `--imag-cutoff` | a mode below −CUTOFF cm⁻¹ is the reaction coordinate (default 50). Reactants and products must have none. |
| `--project` | remove overall translations and rotations exactly before computing frequencies, instead of discarding the six lowest. On by default for xTB and machine-learned-potential inputs. On converged Gaussian and ORCA geometries it changes KIEs by less than 10⁻⁵ (`--no-project` forces it off). |
| `--weights`, `--energies`, `--weight-uncertainty` | how conformers are weighted; see [Several conformers or transition structures](#several-conformers-or-transition-structures). |
| `--job FILE` | a JSON file with the structures, labels and settings, for transition structures in series, parallel pathways and several labelled positions in one run ([docs/job_files.md](docs/job_files.md)). The other flags then give only the outputs. |
| `--calc SPEC`, `--delta` | compute the Hessians of geometry files with an ASE calculator (`xtb`, `mace_mp:medium`, ...); finite-difference step in Å. |
| `--json FILE`, `--csv FILE` | also write the full result as JSON, or append one summary row to a CSV file. |
| `--version` | print the version. |

The full syntax (`python -m kinisot` is equivalent to `kinisot`):

```
kinisot --rct FILE [FILE ...] [--rct FILE ...] (--ts FILE [FILE ...] | --prd FILE [FILE ...])
        --iso ATOMS [--iso ATOMS ...] [-t K] [-s FACTOR] [--imag-cutoff CM-1] [--tunneling MODEL]
        [--barrier KCAL] [--project] [--reference ATOMS] [--calc SPEC]
        [--weights SCHEME] [--energies TABLE] [--weight-uncertainty KCAL]
        [-o FILE] [--overwrite] [-q] [--json FILE] [--csv FILE]
kinisot --job job.json [-o FILE] [--overwrite] [-q] [--json FILE] [--csv FILE]
```

Invalid input stops the run with a message and exit code 1: an atom number
out of range, a reactant with an imaginary frequency, a file that is not a
completed frequency job, and so on.

**Isotopes.** A bare atom number gives the usual heavy label:
¹H → ²H, ¹²C → ¹³C, ¹⁴N → ¹⁵N, ¹⁶O → ¹⁸O, ³²S → ³⁴S, ³⁵Cl → ³⁷Cl,
⁷⁹Br → ⁸¹Br, ²⁸Si → ²⁹Si. Any other isotope can be asked for explicitly:
`--iso 3:17O`, `--iso 7:D` (or `7:T`), `--iso 5:14C`, or an explicit mass
`--iso 5:13.5`. Masses are AME 2020 values for every naturally occurring
isotope of 83 elements plus the common radioactive labels. The light
isotopologue is built from the most abundant isotopes, as Gaussian does,
so Gaussian, ORCA and other programs give the same numbers for the same
Hessian. Kinisot 2.3 and earlier used ¹⁷O for a bare oxygen number; 2.4 and
later use ¹⁸O and print a note once.

## Several conformers or transition structures

One reactant and one transition structure are often enough. When they are
not, Kinisot combines several structures, weighting each by its free
energy (for the light isotopologue, which is exact).

- **Conformers.** Give several files after one flag, all with the same atom
  numbering: `kinisot --rct gs_*.out --ts ts_chair.out ts_boat.out --iso 5`.
  The result is not the Boltzmann average of the pairwise KIEs, but a
  ratio of population-weighted averages. The output lists each conformer's
  population and its own KIE. It also gives the KIE of the lowest
  conformers alone, and how much the result moves if each free energy moves
  by 0.5 kcal/mol. `--weights` chooses quasi-harmonic free energies
  (default), harmonic ones, your own (`--energies`), the lowest conformer,
  or equal weights. [examples/conformers](examples/conformers/README.md)
  works through allyl vinyl ether, where a secondary ²H KIE moves by 0.010
  between the lowest conformer and the ensemble.
- **Transition structures in series.** When an intermediate can go on or
  return, no single step commits the substrate, and the KIE is a weighted
  mean of the steps' KIEs. It depends on the partitioning of the
  intermediate (the commitment factor), which can come from free energies
  or from trajectories.
- **Parallel pathways.** When the substrate is consumed by several routes
  (for example two enantiomers through diastereomeric transition
  structures), each with its own labels, their KIEs combine according to
  each route's share, which can be a measured selectivity.

Series and pathways are set up in a JSON job file
([docs/job_files.md](docs/job_files.md)). The equations, with worked cases
from the Wittig reaction and a Rh-catalysed dynamic kinetic asymmetric
arylation, are in [docs/theory.md, section 7](docs/theory.md#7-conformer-ensembles).

## Input files

| Program | Status | What Kinisot reads |
| --- | --- | --- |
| Gaussian 09/16 | supported | a normally terminated `freq` job: atom masses and the force constants in the archive entry (`opt freq` jobs are fine) |
| ORCA 5/6 | supported | `name.out` plus the `name.hess` file ORCA writes next to it (give either path); level of theory from the `!` line |
| ASE: machine-learned potentials, GFN2-xTB, any ASE calculator | supported (`pip install kinisot[ase]`) | a `VibrationsData` JSON file, or a geometry plus `--calc` (`xtb`, `mace_mp`, `mace_off`, `orb`, `sevennet`, `aimnet2`, `module:callable`); the Hessian is computed and cached next to the geometry, external modes are projected out |

Atom numbers count from 1 in the order of the input geometry, whatever
the program. Details and pitfalls:
[docs/file_formats.md](docs/file_formats.md).

## Why the Bigeleisen–Mayer equation rather than free energies

A KIE can also be taken from the free energies a quantum chemistry program
prints for each isotopologue, as exp(ΔΔG‡/RT). With exact harmonic
frequencies the two routes give the same number, and in practice they
often agree to a few 10⁻⁴. Rzepa's Baeyer–Villiger KIEs came out 1.023
(free energies) against 1.0226 (Bigeleisen–Mayer) for ¹³C, and 0.928
against 0.92831 for ²H
([Henry Rzepa's blog, 2015](https://www.ch.ic.ac.uk/rzepa/blog/?p=14255)).
Against natural-abundance measurements, differences of 10⁻³ matter, and
there the Bigeleisen–Mayer route (used by QUIVER, PyQuiver and Kinisot) is
the more reliable one.

- **Soft vibrations cannot spoil it.** It uses only ratios of the two
  isotopologues' frequencies, and a low-frequency mode contributes a
  factor close to 1 however poorly it is computed. Through free energies
  the same mode's isotope shift enters in full. A 0.1 cm⁻¹ error in the
  isotope shift of one torsion of the Claisen reactant moves the ¹³C KIEs
  by up to 1.5 × 10⁻³ through free energies, and by less than 10⁻⁵ through
  Bigeleisen–Mayer.
- **Printed free energies are too coarse.** Gaussian prints six decimals of
  a hartree. At 393 K, an error of 10⁻⁶ hartree in one of the four free
  energies changes a KIE by 8 × 10⁻⁴.
- **Quasi-harmonic free energies distort it.** Raising soft modes to
  100 cm⁻¹, as quasi-harmonic thermochemistry does, moves the Claisen ¹⁸O
  KIE from 1.037 to 1.018 through free energies. The Bigeleisen–Mayer value
  moves by 10⁻⁴.
- **There are no translational terms to get wrong.** Leaving them out of a
  free-energy or enthalpy–entropy treatment costs 1.3% for
  Cl⁻ + CH₃Br.
- **It is less work.** One frequency calculation per structure gives every
  label at every temperature in seconds.

[docs/theory.md](docs/theory.md) (section 6) has the derivation and the
full comparison, and `scripts/compare_free_energy_route.py` reproduces it.

## What Kinisot assumes

- **Harmonic vibrations and rigid rotors.** There is no anharmonic
  correction. One scaling factor multiplies every frequency, the imaginary
  one included. Kinisot checks that it reproduces the frequencies the
  program printed and warns if not.
- **Transition-state theory.** Each transition structure has one imaginary
  frequency; a second one triggers a warning. Recrossing and dynamic effects
  are outside the model, except through a commitment factor you supply for
  steps in series.
- **One-dimensional tunnelling.** Bell's model is the default; Wigner and
  Skodje–Truhlar corrections are options. Bell's correction is not defined
  below the crossover temperature hc|ν‡|/2πk (229 K for a 1000i cm⁻¹ mode),
  and Kinisot refuses it there. Multidimensional tunnelling, which primary
  hydrogen KIEs can need, is outside Kinisot.
- **Solvent** enters only through the frequency calculation. An
  implicit-solvent frequency job is used as it is.
- **Overall translations and rotations** are removed by discarding the six
  lowest modes (five for a linear molecule), or by projecting them out with
  `--project`. For converged Gaussian geometries the two differ by less
  than 4 × 10⁻⁷ in the KIE. Projection matters for the noisier Hessians of
  finite-difference calculations.
- **Symmetry numbers are left out.** The reduced partition function ratio
  excludes them. If a substitution changes a symmetry number (CH₃ → CH₂D,
  for example), multiply by the ratio yourself.
- **Several structures.** Conformers given together are in fast
  equilibrium (Curtin–Hammett). Transition structures in series are
  combined at steady state, and parallel pathways at low conversion.

## Python API

For scripting and notebooks; everything above is also available from Python.

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

**Conformer ensembles.** A species given as a list of files, or as
`Conformers(files, free_energies=None, degeneracy=None)`, is an ensemble
of conformers with the same atom numbering:

```python
from kinisot import Conformers, compute_kie

r = compute_kie(rct="gs.out", ts=[["ts_chair.out", "ts_boat.out"]], iso="5", temperature=393.0)
r.kie_tunnel                          # sum_i x_i rho_i / sum_j y_j rho'_j
r.conformers                          # per conformer: free energy, population, its own KIE
r.kie_tunnel_lowest, r.n_effective, r.kie_tunnel_range
compute_kie(rct="gs.out", ts=Conformers(["a.out", "b.out"], free_energies=[0.0, 1.2]), iso="5", weights="user")
```

It returns an `EnsembleIsotopeEffect` with the same `kie`, `kie_tunnel`,
`kie_relative`, `to_dict()` and `summary_row()`, but no ZPE/EXC/TRPF
breakdown, which exists per conformer only. One conformer per species gives
the ordinary `IsotopeEffect`. `kinisot.equivalent_positions(results)` gives
the exact isotope effect for positions made equivalent by fast motion (the
hydrogens of a rotating methyl group) from one result per placement of the
label. [docs/theory.md, section 7](docs/theory.md#7-conformer-ensembles)
has the equations.

**Transition structures in series and parallel channels.**

```python
from kinisot import Series, channels, compute_kie, series_kie

# (atom numbers illustrative) steps none of which alone commits the substrate: one label per reactant, then per step
r = compute_kie(rct=["aldehyde.out", "ylide.out"], ts=Series(["ts_4.out", "ts_6.out"], commitment=128 / 76),
                iso=["8", "0", "20", "20"], temperature=340.15)
r.kie_tunnel, r.shares, r.commitment, r.commitment_for(1.033)
# parallel routes with their own files and labels; shares from a measured selectivity
r = channels([dict(rct="3.out", ts="ts_S.out", iso=["1", "40"]), dict(rct="3.out", ts="ts_R.out", iso=["3", "38"])],
             shares=[3.3, 1], temperature=313.15)
r.kie_tunnel, r.selectivity_for(1.0244)
series_kie([1.043, 1.015], free_energies=[25.9, 26.0], temperature=340.15)   # 1.0280, from KIEs you have
```

On the command line these come from a JSON job file, `kinisot --job
job.json` ([docs/job_files.md](docs/job_files.md)), which also runs several
isotopologues at once.

Machine-readable output from the command line: `--json run.json` (the
full result of one run) and `--csv runs.csv` (one row per run, appended).

## Examples and documentation

- [examples/](examples/README.md): the Claisen rearrangement (¹³C, ¹⁸O/¹⁷O
  and ²H KIEs, projection, reference isotopologue, tunnelling models,
  temperature scans), a Diels–Alder reaction with one or two reactant files,
  the same KIE from ORCA and from ASE inputs (with the machine-learned
  potential workflow), a conformational EQE, and a conformer ensemble (eight
  allyl vinyl ether conformers and the chair and boat Claisen transition
  structures), each with the commands, the expected numbers and the
  literature background.
- [docs/theory.md](docs/theory.md), [docs/file_formats.md](docs/file_formats.md),
  [docs/job_files.md](docs/job_files.md),
  [docs/faq.md](docs/faq.md), [docs/comparison.md](docs/comparison.md)
  (PyQuiver, PyQuiverHS, Gaussian's `readisotopes`, GoodVibes).
- [benchmarks/](benchmarks/README.md): computed versus experimental KIEs,
  one directory per reaction, with a runner that regenerates the report.
  At B3LYP/6-31G(d), Kinisot reproduces the measured KIEs of three reactions:
  - the Claisen rearrangement of allyl vinyl ether (five ¹³C and ¹⁷O KIEs,
    mean absolute deviation 0.0009);
  - the Diels–Alder reaction of isoprene with maleic anhydride (nine ¹³C
    and ²H KIEs, mean absolute deviation 0.003);
  - the Shi epoxidation of trans-β-methylstyrene (six ¹³C KIEs at 0 °C,
    mean absolute deviation 0.0012). Kinisot also matches the authors' own
    QUIVER predictions for all 18 published transition structures: 102 of
    108 values exactly at three decimals, the rest within 0.0008.

  Five more systems come from the SI of the PyQuiverHS paper (Grazioli et
  al., 2026), with their measured values: conformational KIEs of three
  biaryl and cyclophane ring flips, a gas-phase SN2 α-secondary KIE, and a
  CD₃ axial/equatorial EIE.
  - Six of the seven measured values are reproduced within 0.016.
  - The SN2 KIE is 0.08 too high: harmonic transition-state theory at
    HF/6-31+G(d) misses it, and PyQuiverHS gives the same value.
  - `tests/test_pyquiverhs.py` checks Kinisot against PyQuiverHS's own
    output for these files over 10–1000 K.
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
