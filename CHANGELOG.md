# Changelog

Notable changes to Kinisot. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [2.5.0] - 2026-09-25

The first release since 2.0.2. It contains everything developed as the
milestones 2.0.3, 2.1.0, 2.2.0, 2.3.0 and 2.4.0 below, none of which was
published separately; read those sections for the details. In short:

- **Correctness**: the atomic-mass-unit typo and duplicated constants
  removed, CODATA 2018 constants, AME 2020 isotope masses, exact
  scaling-factor lookup, linearity from the geometry, and every frequency
  checked against what the program printed.
- **Robustness**: validated input (atom numbers, elements, already
  substituted atoms, structures that are not minima or transition
  structures, mismatched substitutions), errors instead of crashes, a
  results file that is appended to rather than overwritten.
- **Inputs**: Gaussian, ORCA (`.out` + `.hess`) and ASE (`VibrationsData`
  JSON, or any geometry with `--calc` and an ASE calculator, machine-learned
  potentials included), mixable in one run.
- **Physics options**: explicit isotope syntax with a full table (`5:13C`,
  `3:17O`, `7:D`), a bare oxygen index now meaning ¹⁸O, Eckart projection of
  external modes (default for ASE Hessians), Bell, Wigner and
  Skodje–Truhlar tunnelling, a reference isotopologue, harmonic or
  fundamental Truhlar factors, temperature scans.
- **Interfaces**: the `kinisot` console script, `--json`/`--csv` output,
  and a Python API (`compute_kie` returning an `IsotopeEffect`); the 2.0
  functions remain as deprecated shims until 3.0.
- **Project**: `pyproject.toml`, GoodVibes ≥ 4.4 as a dependency, worked
  examples with verified outputs, documentation, a benchmark scaffold, CI
  on three platforms, automated PyPI and GitHub releases.

**Numerical changes relative to 2.0.2** (all documented in the milestone
sections): the constants fix (< 3 × 10⁻⁶ in KIEs), CODATA 2018
(7 × 10⁻⁹), AME 2020 isotope masses (≤ 1.2 × 10⁻⁶), the Truhlar v5
scaling table (a few factors), and ¹⁸O instead of ¹⁷O for a bare oxygen
index (`3:17O` restores the old value).

### Fixed in the release preparation

- The wheel now includes the `kinisot.backends` subpackage (package
  discovery was limited to the top-level package, so an installed copy could
  not be imported; the in-tree test runs did not notice). Verified by
  installing the wheel and running the suite from outside the repository on
  Python 3.9, 3.10, 3.12 and 3.13.

### Phase 8 details (ASE and machine-learned potentials)


- **ASE backend** (`pip install kinisot[ase]`): Hessians from
  `VibrationsData` JSON files, or computed from a geometry with any ASE
  calculator (`--calc mace_mp:medium`, `emt`, `orb`, `sevennet`, `aimnet2`,
  `module:callable`; `--delta` for the finite-difference step, analytic
  `get_hessian` used when available) and cached next to the geometry as
  `<name>.hessian.json`. `kinisot.hessian_from_calculator`,
  `hessian_for_geometry`, `save_hessian_json`, `build_calculator` in the API;
  `compute_kie(..., calculator=)`.
- Projection of external modes is now **on by default for ASE inputs** and
  off for Gaussian/ORCA (`project=None` means "decide by backend";
  `--project`/`--no-project` override). With projection, imaginary modes
  below the cutoff are treated as real vibrations of the same magnitude,
  with a warning.
- `examples/mlip_claisen/` and ASE JSON fixtures (`tests/data/ase/`); the
  Claisen KIE from the JSON form matches the Gaussian path to 10⁻¹⁰.
- `benchmarks/`: a runner and case-file format for comparing computed with
  experimental isotope effects, with the Claisen and Diels–Alder cases
  (experimental values to be entered from the papers).
- The release workflow now also creates a GitHub Release with the changelog
  section as notes; Dependabot watches the Actions versions; `.zenodo.json`
  carries the archive metadata.

## [2.4.0] - development milestone (included in 2.5.0)

Phase 7 of the [implementation plan](IMPLEMENTATION_PLAN.md): isotopes,
external modes, tunnelling.

### Added

- **Isotope table** (`kinisot/isotope_data.py`, generated from the
  `periodictable` package: AME 2020 masses for every naturally occurring
  isotope of 83 elements plus ³H, ¹¹C, ¹⁴C, ¹³N, ¹⁵O, ¹⁸F, ³²P, ³³P, ³⁵S,
  ³⁶Cl, ¹²⁵I, ¹³¹I). Both isotopologues are built from it, so Gaussian, ORCA
  and any other program give the same numbers for the same Hessian.
- **Explicit isotope syntax**: `--iso 5:13C`, `3:17O`, `7:D`, `7:T`, or an
  explicit mass `5:13.5`; a bare number means the default heavy label
  (H, C, N, O, S, Cl, Br, Si). Labels that do not match the atom's element
  are refused.
- **Eckart projection** of translations and rotations (`--project`,
  `project=True`); external modes are removed by value and the reaction
  coordinate is the most negative remaining mode. Needs the geometry, which
  the Gaussian (archive) and ORCA (`$atoms`) readers now provide.
- **Skodje–Truhlar tunnelling** (`--tunneling skodje`, barrier from
  `--barrier KCAL` or from the electronic energies now read from the files:
  Gaussian archive `HF=`, ORCA `FINAL SINGLE POINT ENERGY`).
- **Reference isotopologue** (`--reference ATOMS`, `reference=`): the KIE
  is also reported relative to a second substitution.
- **Temperature scans**: `-t 273,298,323` or `-t 250:350:10` give one
  result line per temperature (and one JSON entry / CSV row each).

### Changed

- **Bare oxygen index means ¹⁸O** (was ¹⁷O); a note is printed once per
  file in this release. `3:17O` gives the old behaviour.
- **Isotope masses** moved from the five-decimal values of Kinisot 1.x/2.x
  (¹H 1.00783, ²H 2.0141, ¹³C 13.00335, ¹⁷O 16.9991) to AME 2020. Largest
  change in a golden KIE: 1.2 × 10⁻⁶ relative (the H7,H8 Claisen case);
  a few bundled results move in the sixth decimal. Golden tests updated.
- Projection versus the lowest-six rule on the bundled Gaussian examples:
  the corrected KIEs differ by 5 × 10⁻⁸ (Claisen C4), 4 × 10⁻⁷ (Claisen
  H7,H8), 2 × 10⁻⁸ (Diels–Alder C15) and 3 × 10⁻⁸ (CD₃ EQE). Projection
  stays off by default for Gaussian/ORCA input through 2.x, as decided.
- The `substitute()` API returns AME masses and records the isotope in each
  `Substitution`; `parse_label()` returns `(index, symbol, mass number,
  mass)` tuples.

## [2.3.0] - development milestone (included in 2.5.0)

Phase 5 of the [implementation plan](IMPLEMENTATION_PLAN.md): GoodVibes
integration and ORCA support.

### Added

- **ORCA input**: give `name.out` or `name.hess` (ORCA writes the Hessian
  to the latter); the level of theory is read from the `!` line. The
  program is detected from the file, so Gaussian and ORCA files can be
  mixed. ORCA's standard atomic weights are replaced by pure-isotope masses
  for the elements Kinisot can substitute, so the same Hessian gives the
  same numbers from either program.
- **Frequency self-check**: for every unsubstituted species Kinisot compares
  the frequencies it obtains from the Hessian with the ones the program
  printed and warns when they differ by more than 1 cm⁻¹.
- `--scale-type zpe|harm|fund` (and `scale_type=` in the API) to apply the
  harmonic or fundamental Truhlar factor instead of the ZPE one.
- Geometry-based linearity test (`kinisot.hessian.linear_from_geometry`),
  used for ORCA and available to every backend; `HessianInput` gains
  `program_frequencies`.
- `examples/orca_claisen/` and ORCA-layout test fixtures (`tests/data/orca/`).

### Changed

- **GoodVibes ≥ 4.4 is a dependency.** Scaling factors now come from the
  Truhlar database version 5 shipped with GoodVibes, with its cross-program
  canonicalization (Gaussian `PBE1PBE` = ORCA `PBE0`, `6-31G*` = `6-31G(d)`,
  ...). Differences from the version 3b2 table Kinisot carried before:
  M06-2X/6-31+G(d,p) 0.967 → 0.968; `M06-2X/maug-cc-pVTZ` and
  `PW6B95/6-31+G(d,p)` are no longer listed (factor 1.0 with a message);
  four spellings (`MN15-L`, `MN12-L`, `MN12-SX`, `M06-L(DKH2)`) are aliased
  by Kinisot until GoodVibes carries them. The bundled examples are
  unaffected (B3LYP/6-31G(d) → 0.977 in both).
- `kinisot/vib_scale_factors.py` is gone; `kinisot.scaling.find_scaling_factor`
  takes an optional scale type.

## [2.2.0] - development milestone (included in 2.5.0)

Phase 4 of the [implementation plan](IMPLEMENTATION_PLAN.md): internal
refactor and Python API. Numerically identical to 2.1.0 apart from the
constants update listed under Changed.

### Added

- `kinisot.compute_kie()` returning a frozen `IsotopeEffect` (final factors,
  per-side `SideResult`s with the light and heavy isotopologues, their kept,
  discarded and imaginary frequencies, masses and substitutions, the scaling
  choice and any warnings) with `to_dict()`, `to_json()` and
  `summary_row()`. Inputs may be file paths or `HessianInput` objects.
- `--json FILE` (full result of a run) and `--csv FILE` (one row per run).
- `--tunneling bell|wigner|none` on the command line and `tunneling=` in the
  API. The Bell correction is refused below the crossover temperature,
  where it diverges, instead of returning a meaningless number.
- `HessianInput`, the program-independent interchange type every backend
  produces (Hessian, masses, atomic numbers, level of theory, linearity,
  positions), in `kinisot/hessian.py`; `kinisot/backends/gaussian.py`
  holds the Gaussian parser (now also reads the archive geometry).
- `examples/api_example.py` and an executed `examples/examples.ipynb`.

### Changed

- Module layout: physics in `kinisot/thermo.py`, scaling in
  `kinisot/scaling.py`, isotopes in `kinisot/isotopes.py`, the CLI in
  `kinisot/cli.py`, the API in `kinisot/api.py`. `kinisot.Kinisot` and
  `kinisot.Hess_to_Freq` remain as deprecation shims (removed in 3.0):
  `compute_isotope_effect()`, `calc_rpfr`, `get_frequency_scaling()` and
  the `calc_*_factor` functions warn with `DeprecationWarning` and forward
  to the new code.
- Physical constants updated from CODATA 2010 to CODATA 2018 (h, k, E_h,
  a_0, u). Largest change in any golden quantity: 7 × 10⁻⁹ relative; no
  printed digit of the bundled examples changes.
- The scaling-factor table is a plain dictionary of named tuples
  (`kinisot.vib_scale_factors.SCALING_FACTORS`) instead of a NumPy
  structured array; factors are now exact decimals rather than float32.

## [2.1.0] - development milestone (included in 2.5.0)

Phases 2 (robustness) and 3 (packaging and documentation) of the
[implementation plan](IMPLEMENTATION_PLAN.md). Every result line of the
bundled examples is unchanged (to the 6 printed decimals, modulo the 2.0.3
constants fix).

### Added

- `kinisot` console script (`python -m kinisot` still works) and a
  `pyproject.toml` (PEP 621) build; `setup.py`/`setup.cfg` are gone. The
  version is defined once, in `kinisot/__init__.py`, which now exports the
  public API (`compute_isotope_effect`, `parse_gaussian`, the exception
  classes, ...).
- Worked examples in `examples/` (Claisen rearrangement, Diels–Alder with
  one or two reactant files, conformational EQE) with the commands, the
  expected output and literature background; `examples/run_examples.sh`
  regenerates the outputs and `tests/test_examples.py` checks them.
- `docs/theory.md` (the equations as implemented), `docs/file_formats.md`,
  `docs/faq.md`, `docs/comparison.md`; a rewritten README with a quick
  start, a guide to every output column and the Python API.
- `CITATION.cff`, `CONTRIBUTING.md`, issue templates, `ruff` lint and
  format checks in CI, and a tag-driven PyPI publishing workflow
  (`.github/workflows/publish.yml`, trusted publishing).

- Input validation with clear messages: `--iso` atom numbers out of range,
  `0` combined with atom numbers, duplicated atoms, atoms of an element
  Kinisot cannot substitute, atoms that already carry a heavy isotope in the
  Gaussian input, a reactant or product with an imaginary frequency, a
  transition-structure side with more than one file carrying an imaginary
  frequency, different substituted elements on the two sides of the reaction
  (e.g. `13C` in the reactant but `2H` in the TS), and labels that request
  no substitution at all. Previously most of these were silently ignored.
- A warning when a transition structure has more than one imaginary
  frequency beyond the cutoff; a negative frequency can no longer reach the
  partition functions unnoticed (`nan` results).
- The results table now lists, for every species, the modes kept in the
  partition function and the modes discarded as external (translation and
  rotation), so a misassigned low-frequency mode is visible.
- `--output FILE` (`-o`), `--overwrite`, `--quiet` (`-q`), `--version`,
  `--imag-cutoff` (the old spelling `--cutoff` still works), and long forms
  `--temperature` and `--scale`.
- `kinisot.exceptions` (`KinisotError`, `KinisotParseError`,
  `KinisotInputError`, `KinisotWarning`); library code raises these instead
  of calling `sys.exit()`. Both error classes also subclass `ValueError`.
- `Hess_to_Freq.parse_gaussian()` returns a `FrequencyData` record (Hessian,
  masses, atomic numbers, level of theory, linearity) from one read of the
  file; `substitute()` reports which atoms were changed.
- Tests: frequencies of every bundled output are checked against the values
  Gaussian prints (kept modes to 0.05 cm-1, discarded modes against
  Gaussian's low frequencies); error paths; CLI behaviour including a
  `python -m kinisot` run. Coverage is reported in CI and must stay at or
  above 80 %.

### Changed

- The Gaussian outputs used by the tests and examples moved from
  `kinisot/examples/` to `tests/data/`; the wheel now weighs 23 kB instead
  of shipping 2.8 MB of test data.
- **Results file**: new results are appended to `Kinisot_output.dat`
  (or `--output`) instead of silently overwriting it, matching how the
  bundled example scripts collect a series of substitutions; use
  `--overwrite` to start afresh. The `Species: ... isotopologue: ...` lines
  are now written to the results file, not only to the terminal.
- **Per-species rows** of the results table now show each species' own
  Bigeleisen-Mayer factors, i.e. light/heavy ratios of the ZPE and
  excitation terms and heavy/light ratio of the frequency product (TRPF).
  Before 2.1.0 the TRPF column of the two rows was swapped (the reactant
  row showed the transition-structure quantity). The final `KIE @` line is
  unchanged and still equals the ratio of the two rows in every column.
- File names in the table are shown without directory and extension; a
  multi-file side is shown as `a + b`. The result line of an equilibrium
  isotope effect is labelled `EQE @` instead of `KIE @`.
- The level-of-theory check covers all files (not only the first two) and
  its warning is printed to the terminal.
- Unknown command-line flags are rejected (`parse_args`), `--ts` together
  with `--prd` is rejected, and a non-positive temperature or scaling
  factor is rejected. `-s` no longer treats `0` as "auto-detect".
- Windows Gaussian archives (`|` separators) are parsed; archive entries
  wrapped across lines are joined before parsing.
- The Bigeleisen-Mayer terms are evaluated with NumPy array operations.

### Removed

- The `try/except` import shim: Kinisot is run as `python -m kinisot`
  (or, from Phase 3, the `kinisot` console script); `python Kinisot.py`
  is no longer supported.
- `Logger.Fatal()`; the logger is now a context manager.

## [2.0.3] - development milestone (included in 2.5.0)

### Fixed

- **Physical constants**: removed a duplicated constants block that silently
  overrode the correct values with lower-precision ones, including a
  transcription error in the atomic mass unit (1.660468e-27 kg instead of
  CODATA 1.660539e-27 kg). Harmonic frequencies shift by ~2×10⁻⁵ relative;
  computed KIE/EQE values shift by <3×10⁻⁶ because the error largely cancels
  in the Bigeleisen-Mayer ratios. Results printed to 6 decimal places are
  unchanged for all bundled examples except in the last digit.
- **Crash on TS without an imaginary frequency**: passing a ground-state file
  as `--ts` crashed with `NameError`; it now reports a clear error message.
- **Vibrational scaling factor auto-detection**: the level of theory is now
  matched exactly against the Truhlar database (after normalizing case,
  hyphens, and Gaussian's R/U/RO spin prefix). Previously substring matching
  could apply a wrong factor (e.g. plain B3LYP factors to CAM-B3LYP jobs) and
  the last matching table row silently won.
- **Scaling-factor table**: the entry mislabeled `M06/maug-cc-pVTZ`
  (zpe_fac 0.971) is now correctly `M06-2X/maug-cc-pVTZ`, verified against
  the Truhlar database v3b2.
- **Linear molecule detection**: linearity is now determined from the
  rotational constants instead of the "prolate symmetric top" string, which
  non-linear prolate tops (e.g. CH3Cl) also print — those molecules would
  have had one contaminating near-zero mode included in the partition
  function products.

### Added

- Characterization test suite (23 tests) covering Claisen KIEs,
  single/multi-reactant Diels-Alder KIEs, the EQE (`--prd`) path,
  scaling-factor lookup, and linearity detection.
- GitHub Actions CI (Python 3.9-3.13 on Linux/macOS/Windows), replacing the
  defunct Travis setup.

## [2.0.2] - 2021

Last release before this changelog was introduced.
