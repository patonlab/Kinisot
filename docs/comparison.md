# Comparison with other isotope-effect codes

Last checked 2026-09-25 (details and verification in [REVIEW.md](../REVIEW.md)).

## PyQuiver

[PyQuiver](https://github.com/ekwan/PyQuiver) (`pip install pyquiver-kie`,
Kwan group, Apache-2.0) is the closest peer: Kinisot is a rewrite of Rzepa's
Fortran Kinisot, PyQuiver a Python port of QUIVER, and both implement the
same Bigeleisen–Mayer treatment. Fed the same Hessians, the two agree to
2.3 × 10⁻⁶ or better on every Claisen and Diels–Alder KIE, uncorrected and
with Bell or Wigner tunnelling; `tests/test_pyquiver.py` asserts agreement
within 3 × 10⁻⁶ for all of these in CI.
The residual comes from PyQuiver's five-decimal isotope masses and its
atomic mass unit:

| Aspect | Kinisot 2.1 | PyQuiver |
| --- | --- | --- |
| Hessian source | Gaussian archive, ORCA `.hess` | Gaussian archive (`#p` output required), ORCA `.hess` |
| External modes | not projected, lowest 5/6 dropped | same |
| Imaginary-mode threshold | 50 cm⁻¹ | 50 cm⁻¹ |
| Scaling | one factor, all modes, after diagonalization | same |
| Bigeleisen–Mayer terms | logarithmic sums | direct products (algebraically identical) |
| Tunnelling | Bell | Wigner, Bell, Skodje–Truhlar |
| Constants | CODATA 2010 | atomic mass unit 1.660468 × 10⁻²⁷ kg (a transcription error, removed from Kinisot in 2.0.3) |
| Isotopes | ²H, ¹³C, ¹⁷O by atom number | ²H ³H ¹³C ¹⁴C ¹⁵N ¹⁷O ¹⁸O ¹⁸F by label |
| Reference isotopologue | no | yes |
| Several reactant files | yes | no (one ground-state file) |
| Scaling-factor lookup | automatic (Truhlar v5 via GoodVibes; ZPE, harmonic or fundamental) | manual |
| Output | text table | text, CSV, pandas |

## PyQuiverHS

A web front end to the PyQuiver family with Bigeleisen–Mayer and
enthalpy–entropy partition functions, Wigner and Bell corrections and
temperature scans (Grazioli et al., J. Phys. Org. Chem. 2026, 39, e70099,
[doi:10.1002/poc.70099](https://doi.org/10.1002/poc.70099)). Its
1,1,3,3-tetramethylcyclohexane EIE demonstration is the same system as
Kinisot's [eqe_cyclohexane](../examples/eqe_cyclohexane/README.md) example.

## Gaussian `freq=readisotopes`

Gaussian can print the thermochemistry of an isotopologue directly. This
projects external modes and uses Gaussian's own constants; the resulting
KIEs agree with Kinisot to the printed precision. Kinisot exists so that
every position and temperature does not need a new Gaussian job.

## GoodVibes

[GoodVibes](https://github.com/patonlab/GoodVibes) (same group) computes
quasi-harmonic thermochemistry, not isotope effects, but from version 4.4 it
parses Cartesian Hessians (Gaussian, ORCA) and program-independent levels
of theory and scaling factors. Kinisot depends on it since 2.3: the Truhlar
factors, the ORCA `$hessian` reader and ORCA level-of-theory detection come
from GoodVibes.
