![Kinisot Banner](https://github.com/patonlab/Kinisot/blob/master/kinisot_banner.png)

[![DOI](https://zenodo.org/badge/54840251.svg)](https://zenodo.org/badge/latestdoi/54840251)
[![PyPI version](https://badge.fury.io/py/kinisot.svg)](https://badge.fury.io/py/kinisot)
[![CI](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml/badge.svg)](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml)

***
## Introduction

**Kinisot** is a Python program to compute kinetic (KIE) and equilibrium (EIE) isotope effects from two output files (Gaussian, or ORCA with its `.hess` file alongside): a ground state and either a transition state (for a KIE) or a product (for an EIE). File parsing and vibrational scaling factors are provided by [GoodVibes](https://github.com/patonlab/goodvibes). It is developed in the [Paton research group](https://patonlab.colostate.edu) at Colorado State University.

This is a Python version of [Kinisot](http://dx.doi.org/10.5281/zenodo.19272), inspired by the Fortran version originally written by [Henry Rzepa](https://en.wikipedia.org/wiki/Henry_Rzepa). This version does not require parameter files to run and allows easy manipulation of temperature and [vibrational scaling factors](http://t1.chem.umn.edu/freqscale/index.html). The level of theory and basis set are detected from in the output files and the program will attempt to assign the appropriate scaling factor based on data from the [Truhlar group](https://t1.chem.umn.edu/freqscale/index.html). Isotopic substitutions are to be specified by the command line, such that separate computations with Gaussian are not required. The program diagonalizes the mass-weighted Hessian matricies to obtain harmonic vibrational frequencies and Bigeleisen-Mayer Reduced Isotopic Partition Function Ratios. One difference with the Gaussian program itself is that the lowest five/six normal modes (translations and rotations) are not projected out, however, this is also the approach taken by [quiver](https://github.com/ekwan/quiver). Testing  this approach against calculation of the Reduced Isotopic Partition Function Ratios by hand (!) using the frequencies in the Gaussian output files led to agreeement up to 4DP. Three one-dimensional tunneling corrections are reported: the Bell infinite-parabola and Wigner models, and (when a barrier height is available from the input files' electronic energies, or supplied with `--barrier`) the Skodje-Truhlar truncated-parabola model. All three are unreliable below the crossover temperature; see the Usage notes.

Also see related discussions on [computing KIE values](http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327)

Isotopic substitutions available are <sup>2</sup>D and <sup>3</sup>T (for H), <sup>13</sup>C and <sup>14</sup>C, <sup>15</sup>N, and <sup>17</sup>O and <sup>18</sup>O. A plain atom number uses that element's default heavier isotope (<sup>2</sup>D, <sup>13</sup>C, <sup>15</sup>N, <sup>17</sup>O); a specific isotope is chosen with a suffix such as `--iso 5:18O`. The table lives in `isotopes.py`.

A video guide to using an older version this software is available at Youtube:

[![Kinisot Video Guide](http://img.youtube.com/vi/r4x2gmkc0U8/0.jpg)](http://www.youtube.com/watch?v=r4x2gmkc0U8)


## Installation

To install **Kinisot** with [conda](https://anaconda.org/conda-forge/kinisot):
```
conda install kinisot -c conda-forge
```
To install **Kinisot** with [pypi](https://pypi.org/project/kinisot/):
```
pip install kinisot
```

## Usage

```
kinisot --rct reactant_output --ts ts_output --iso 1,2,3 (-t temperature) (-s scalefactor)
```

(`python -m kinisot ...` works too.)

*	The two output files contain frequency calculations performed for the reactant and transition state at the same level of theory. Gaussian output files and ORCA outputs (with the `.hess` file next to the `.out`, or the `.hess` passed directly) are supported.
*	For an equilibrium isotope effect, pass `--prd product_output` instead of `--ts`.
*	The `--iso` flag is required and specifies comma-separated atom number(s) which are to be substituted for heavier isotopes. Plain atom numbers get the default heavy isotope of their element (²D, ¹³C, ¹⁵N, ¹⁷O); an explicit isotope can be chosen with a suffix, e.g. `--iso 5:18O,7:2D` (available: 2D, 3T, 13C, 14C, 15N, 17O, 18O). With several input files, repeat the flag once per file (`--iso 0` for a file with no substitution); a single `--iso` is applied to reactant and TS/product alike.
*	Several isotopologues can be computed in one run by separating them with `;` and (optionally) naming them: `--iso "C5=5;C4=4;HD=7,8"`. Adding `--ref C5` divides all other isotope effects by the one at C5, matching the internal-standard referencing of experimental KIE measurements.
*	The `-t` option specifies temperature (in Kelvin). N.B. This does not have to correspond to the temperature used in the underlying calculation since the Reduced Isotopic Partition Function Ratios are evalulated at the requested temperature. The default value is 298.15 K.
*	The `-s` option is a scaling factor for vibrational frequencies. Empirical ZPE-scaling factors from the [Truhlar group database](https://comp.chem.umn.edu/freqscale/) are applied automatically (via GoodVibes) based on detection of the level of theory and basis set in the output files. The default value when no scaling factor is available is 1 (no scale factor).
*	Bell infinite-parabola and Wigner one-dimensional tunneling corrections are reported alongside the uncorrected isotope effect. A Skodje-Truhlar correction (truncated parabola) is added automatically when a barrier height can be obtained from the input files' electronic energies, or set explicitly with `--barrier HARTREE`; each isotopologue uses its own ZPE-corrected barrier. All three 1-D corrections become unreliable in the deep-tunneling regime below the crossover temperature T_c = hν‡/2πk_B — the Bell value is reported as `nan` there, and Kinisot warns.
*	`--output PATH` sets the `.dat` file location, `--json PATH` additionally writes machine-readable results, and `-q`/`--quiet` suppresses terminal printing.

Programmatic use returns structured results:

```python
import kinisot
effect = kinisot.compute_isotope_effect(["gs.out"], ["ts.out"], None, ["5", "5"],
                                        temperature=393, freq_scale_factor=0.9614)
effect.kie, effect.kie_wigner, effect.kie_bell
```

See examples/ for more examples
