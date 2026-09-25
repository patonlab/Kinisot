# Changelog

Notable changes to Kinisot. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [2.2.0] - Unreleased

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

## [2.1.0] - Unreleased

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

## [2.0.3] - Unreleased

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
