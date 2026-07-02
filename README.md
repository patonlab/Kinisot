![Kinisot Banner](https://github.com/patonlab/Kinisot/blob/master/kinisot_banner.png)

[![DOI](https://zenodo.org/badge/54840251.svg)](https://zenodo.org/badge/latestdoi/54840251)
[![PyPI version](https://badge.fury.io/py/kinisot.svg)](https://badge.fury.io/py/kinisot)
[![CI](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml/badge.svg)](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml)

***
## Introduction

**Kinisot** is a Python program to compute kinetic isotope effects from two output files (Gaussian, or ORCA with its `.hess` file alongside), one of which is a ground state and the other a transition state. File parsing and vibrational scaling factors are provided by [GoodVibes](https://github.com/patonlab/goodvibes). It is developed in the [Paton research group](https://patonlab.colostate.edu) at Colorado State University.

This is a Python version of [Kinisot](http://dx.doi.org/10.5281/zenodo.19272), inspired by the Fortran version originally written by [Henry Rzepa](https://en.wikipedia.org/wiki/Henry_Rzepa). This version does not require parameter files to run and allows easy manipulation of temperature and [vibrational scaling factors](http://t1.chem.umn.edu/freqscale/index.html). The level of theory and basis set are detected from in the output files and the program will attempt to assign the appropriate scaling factor based on data from the [Truhlar group](https://t1.chem.umn.edu/freqscale/index.html). Isotopic substitutions are to be specified by the command line, such that separate computations with Gaussian are not required. The program diagonalizes the mass-weighted Hessian matricies to obtain harmonic vibrational frequencies and Bigeleisen-Mayer Reduced Isotopic Partition Function Ratios. One difference with the Gaussian program itself is that the lowest five/six normal modes (translations and rotations) are not projected out, however, this is also the approach taken by [quiver](https://github.com/ekwan/quiver). Testing  this approach against calculation of the Reduced Isotopic Partition Function Ratios by hand (!) using the frequencies in the Gaussian output files led to agreeement up to 4DP. A one-dimensional tunneling correction is also included, which is the Bell infinite-parabola model.

Also see related discussions on [computing KIE values](http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327)

The current version is currently hard-coded to consider <sup>2</sup>D/<sup>1</sup>H, <sup>13</sup>C/<sup>12</sup>C and <sup>17</sup>O/<sup>16</sup>O isotopic replacements. This can be modified in Hess_to_Freq.py.

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
*	The `--iso` flag is required and specifies comma-separated atom number(s) which are to be substituted for heavier isotopes. With several input files, repeat the flag once per file (`--iso 0` for a file with no substitution); a single `--iso` is applied to reactant and TS/product alike.
*	The `-t` option specifies temperature (in Kelvin). N.B. This does not have to correspond to the temperature used in the underlying calculation since the Reduced Isotopic Partition Function Ratios are evalulated at the requested temperature. The default value is 298.15 K.
*	The `-s` option is a scaling factor for vibrational frequencies. Empirical ZPE-scaling factors from the [Truhlar group database](https://comp.chem.umn.edu/freqscale/) are applied automatically (via GoodVibes) based on detection of the level of theory and basis set in the output files. The default value when no scaling factor is available is 1 (no scale factor).

See examples/ for more examples
