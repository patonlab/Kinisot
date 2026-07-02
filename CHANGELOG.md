# Changelog

Notable changes to Kinisot. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
