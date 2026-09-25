![Kinisot Banner](https://github.com/patonlab/Kinisot/blob/master/kinisot_banner.png)

[![DOI](https://zenodo.org/badge/54840251.svg)](https://zenodo.org/badge/latestdoi/54840251)
[![PyPI version](https://badge.fury.io/py/kinisot.svg)](https://badge.fury.io/py/kinisot)
[![CI](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml/badge.svg)](https://github.com/patonlab/Kinisot/actions/workflows/ci.yml)

***
## Introduction

**Kinisot** is a Python program to compute kinetic isotope effects from two Gaussian output files, one of which is a ground state and the other a transition state. It is developed in the [Paton research group](https://patonlab.colostate.edu) at Colorado State University.

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
python -m kinisot --rct reactant.out --ts ts.out --iso 5 [-t 393] [-s 0.961]
python -m kinisot --rct reactant.out --prd product.out --iso 5          (equilibrium isotope effect)
python -m kinisot --rct diene.out --rct dienophile.out --ts ts.out --iso 6 --iso 0 --iso 15   (bimolecular)
```

*	The output files contain Gaussian frequency calculations performed for the reactant(s) and the transition structure (`--ts`, kinetic isotope effect) or product (`--prd`, equilibrium isotope effect) at the same level of theory.
*	`--iso` gives the atom number(s) to replace with the heavy isotope (<sup>2</sup>H, <sup>13</sup>C or <sup>17</sup>O), comma separated for several atoms (`--iso 7,8`). Give one `--iso` per file, in the order of the `--rct` files followed by the `--ts`/`--prd` file(s), or a single `--iso` when the atom numbering is the same in all files. Use `--iso 0` for a file with no substituted atom (e.g. the second reactant of a bimolecular reaction). Kinisot checks that the atom numbers exist, that the atoms can be substituted, and that both sides of the reaction substitute the same elements.
*	`-t` sets the temperature in Kelvin at which the reduced isotopic partition function ratios are evaluated (default 298.15 K). It does not have to match the temperature used in the Gaussian calculation.
*	`-s` sets the vibrational scaling factor. When it is omitted the level of theory is detected from the output files and the ZPE scaling factor from the [Truhlar group database](https://comp.chem.umn.edu/freqscale/) is applied; if the level is not in the database, or the files disagree, the factor is 1.0 and a message says so.
*	`--imag-cutoff` (default 50 cm<sup>-1</sup>): a mode below this value counts as the reaction coordinate of the transition structure. Reactants and products must not have one.
*	Results are printed to the terminal and appended to `Kinisot_output.dat` (change the file with `-o`, start afresh with `--overwrite`, silence the terminal with `-q`). Each block lists the frequencies of the reaction-coordinate mode, the Bigeleisen-Mayer factors of every species and, on the `KIE @` line: the ratio of imaginary frequencies (V-ratio), the ZPE, excitation (EXC) and Teller-Redlich product (TRPF) factors, the semiclassical KIE, the Bell tunnelling correction (1D-tunn) and the tunnelling-corrected KIE (corr-KIE). The modes kept in, and discarded from, each partition function are listed below the table.
*	Invalid input (an atom number out of range, a reactant with an imaginary frequency, different substitutions on the two sides, a file that is not a completed frequency job, ...) stops the run with a message and exit code 1.

See [kinisot/examples/gaussian](kinisot/examples/gaussian) for worked examples with their reference outputs.
