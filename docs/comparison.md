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
[doi:10.1002/poc.70099](https://doi.org/10.1002/poc.70099); preprint
[doi:10.26434/chemrxiv.15004570](https://doi.org/10.26434/chemrxiv.15004570)).
Its 1,1,3,3-tetramethylcyclohexane EIE demonstration is the same system as
Kinisot's [eqe_cyclohexane](../examples/eqe_cyclohexane/README.md) example.

The Gaussian files and PyQuiverHS outputs of the preprint's SI are
Kinisot benchmark cases, with the authors' agreement:
- `dihydrophenanthrene`, `biaryl_diketone` and `metaparacyclophane`
  (conformational KIEs);
- `sn2_chloride_methyl_bromide`;
- `tetramethylcyclohexane_eie`.

`tests/test_pyquiverhs.py` runs Kinisot with the same labels, scaling
factors and imaginary-frequency thresholds, over their 10–1000 K scans:
- **Agreement.** Every Bigeleisen–Mayer term and total (uncorrected, Wigner
  and Bell), for all 11 isotopologues, agrees within 0.01 K / T, i.e.
  3 × 10⁻⁵ at 300 K. The difference grows as 1/T through the zero-point
  term, as with PyQuiver above.
- **Below the crossover temperature**, where Bell's formula does not apply,
  Kinisot stops with an error. PyQuiverHS prints a value there, negative
  for the 45i cm⁻¹ biaryl transition structure at 10 K.
- **Enthalpy–entropy terms.** PyQuiverHS's H-S terms contain vibrational and
  rotational entropy but no translational term. That cancels when the same
  atoms are on both sides. It does not when an unlabelled reactant is left
  out of the input, as the chloride is in the SN2 case: there the H-S KIE
  is 1.3% below the Bigeleisen–Mayer one (0.877 against 0.888 at 300 K;
  [theory.md, section 6](theory.md#6-bigeleisenmayer-versus-free-energy-differences)).

## Other programs

- **QUIVER**, PyQuiver's ancestor, first written by Keith Laidig around 1988
  (according to Rzepa's post below).
- **THERMISTP** (Saunders, Wolfsberg, Anet and Kronja, J. Am. Chem. Soc.
  2007, 129, 10276, [doi:10.1021/ja072375r](https://doi.org/10.1021/ja072375r)).
- **ISOEFF** (Anisimov and Paneth, J. Math. Chem. 1999, 26, 75,
  [doi:10.1023/a:1019173509273](https://doi.org/10.1023/a:1019173509273)).

These are the older Bigeleisen–Mayer programs. Rzepa's account of
reproducing a 20-year-old QUIVER calculation (fixed array sizes, a
preprocessor available only as an executable) is a reminder of why
Kinisot keeps its test inputs and outputs in the repository
([Henry Rzepa's blog, 2015](https://www.ch.ic.ac.uk/rzepa/blog/?p=14255)).
Rzepa's own program, which Kinisot descends from, was archived in 2015 as
the Zenodo record Kinisot is cited by
([doi:10.5281/zenodo.19272](https://doi.org/10.5281/zenodo.19272)).

## Gaussian `freq=readisotopes`

Gaussian can print the thermochemistry of an isotopologue directly, and a
KIE can be built from those free energies. That route differs from
Kinisot's Bigeleisen–Mayer evaluation in two ways. First, by the
Teller–Redlich product-rule violation of the computed frequencies: up to
8 × 10⁻⁴ for the Claisen example. Second, by the rounding of the printed
free energies: 8 × 10⁻⁴ per 10⁻⁶ hartree at 393 K. It is also far more
sensitive to errors in soft modes; see
[theory.md, section 6](theory.md#6-bigeleisenmayer-versus-free-energy-differences).
The `freq=readisotopes` route also needs a new Gaussian job for every
position and temperature; Kinisot reuses one Hessian per structure and
only changes the masses.

## GoodVibes

[GoodVibes](https://github.com/patonlab/GoodVibes) (same group) computes
quasi-harmonic thermochemistry, not isotope effects, but from version 4.4 it
parses Cartesian Hessians (Gaussian, ORCA) and program-independent levels
of theory and scaling factors. Kinisot depends on it since 2.3: the Truhlar
factors, the ORCA `$hessian` reader and ORCA level-of-theory detection come
from GoodVibes. Its quasi-harmonic free energies are meant for
thermochemistry and must not be used to compute isotope effects: raising
soft modes to 100 cm⁻¹ turns the Claisen ¹⁸O KIE of 1.037 into 1.018
([theory.md, section 6](theory.md#6-bigeleisenmayer-versus-free-energy-differences)).
