# Changelog

Notable changes to Kinisot. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **ORCA support**: Kinisot now reads ORCA frequency jobs (pass the `.out`;
  the companion `.hess` file is found automatically, or pass the `.hess`
  directly). Isotope substitution recognizes both Gaussian's isotopic masses
  and ORCA's average atomic masses. New `examples/orca/` n-pentane conformer
  EQE example and end-to-end tests.

### Changed

- **File parsing and scaling factors are now delegated to
  [GoodVibes](https://github.com/patonlab/goodvibes)** (new dependency,
  >= 4.4): the Cartesian Hessian and isotope-aware per-atom masses come from
  `goodvibes.io.parse_hessian`, level-of-theory detection and linearity come
  from GoodVibes' program-agnostic parsers, and the vendored Truhlar v3b2
  scaling-factor table is replaced by the **Truhlar v5 database** bundled
  with GoodVibes. Auto-detected scaling factors may differ where the v5
  database revised or renamed entries (e.g. `M06-2X/maug-cc-pVTZ` is now
  `M062X/maug-cc-pV(T+d)Z`). Explicitly supplied `-s` values are unaffected.
- Linearity is now determined from the point group detected by the QC
  program rather than parsed rotational constants.
- Requires Python >= 3.9 (matching GoodVibes).

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
