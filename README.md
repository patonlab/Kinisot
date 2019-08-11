![Kinisot Banner](kinisot_banner.png)

A Python program to compute kinetic isotope effects from two Gaussian output files, one of which is a ground state and the other a transition state. It is developed by [Robert Paton](https://orcid.org/0000-0002-0104-4166) at Colorado State University.

This is a Python version of [Kinisot](http://dx.doi.org/10.5281/zenodo.19272), inspired by the Fortran version originally written by [Henry Rzepa](https://en.wikipedia.org/wiki/Henry_Rzepa). This version does not require parameter files to run and allows easy manipulation of temperature and [vibrational scaling factors](http://t1.chem.umn.edu/freqscale/index.html). The level of theory and basis set are detected from in the output files and the program will attempt to assign the appropriate scaling factor based on data from the [Truhlar group](https://t1.chem.umn.edu/freqscale/index.html). Isotopic substitutions are to be specified by the command line, such that separate computations with Gaussian are not required. The program diagonalizes the mass-weighted Hessian matricies to obtain harmonic vibrational frequencies and Bigeleisen-Mayer Reduced Isotopic Partition Function Ratios. One difference with the Gaussian program itself is that the lowest five/six normal modes (translations and rotations) are not projected out, however, this is also the approach taken by [quiver](https://github.com/ekwan/quiver). Testing  this approach against calculation of the Reduced Isotopic Partition Function Ratios by hand (!) using the frequencies in the Gaussian output files led to agreeement up to 4DP. A one-dimensional tunneling correction is also included, which is the Bell infinite-parabola model.

Also see related discussions on [computing KIE values](http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327)

The current version is currently hard-coded to consider <sup>2</sup>D/<sup>1</sup>H, <sup>13</sup>C/<sup>12</sup>C and <sup>17</sup>O/<sup>16</sup>O isotopic replacements. This can be modified in Hess_to_Freq.py. 

A video guide to using this software is available at Youtube:

[![Kinisot Video Guide](http://img.youtube.com/vi/r4x2gmkc0U8/0.jpg)](http://www.youtube.com/watch?v=r4x2gmkc0U8)

#### Installation
1. Clone the repository https://github.com/bobbypaton/Kinisot.git
2. Add the directory (/pathto/GoodVibes/goodvibes) containing the python files to the PATH environment variable (optional). 
3. Run the script with your Gaussian output files. It has been tested with Python 2 and 3 on OSX

Alternatively `pip install kinisot` will install all classes

**Correct Usage**

```python
Kinisot.py reactant_output ts_output -iso "atom numbers" (-t temperature) (-s scalefactor)  
```
*	The two output files contain Gaussian frequency calculations performed for the reactant and transition state at the same level of theory. The temperature is unimportant.
*	The `-iso` flag is required and specifies a string of atom number(s) which are to be substituted for heavier isotopes. Multiple atom numbers require quotation marks and are separated by spaces.
*	The `-t` option specifies temperature (in Kelvin). N.B. This does not have to correspond to the temperature used in the Gaussian calculation since the Reduced Isotopic Partition Function Ratios are evalulated at the requested temperature. The default value is 298.15 K.
*	The `-s` option is a scaling factor for vibrational frequencies. Empirical scaling factors have been determined for several functional/basis set combinations, and these are applied automatically using values from the Truhlar group based on detection of the level of theory and basis set in the output files. The ZPE-scaling factors are selected if available. The default value when no scaling factor is available is 1 (no scale factor).

#### Example 1. Proton Transfer
Transfer of a proton between two chloride anions; the reactant and TS are both linear. The hydrogen atom is atom number 1 in both files (the numbering needs to be identical). In this example we consider the 2D/1H kie at the default temperature without vibrational scaling.

```python
python Kinisot.py Cl_H_Cl_RCT.out Cl_H_Cl_TS.out -iso "1" -s 1.0
```

The following output is produced:

```bash
                               Temp = 298.15K / Vib. scale factor = 1.0 
                               V-ratio         ZPE         EXC        TRPF         KIE     1D-tunn    corr-KIE 
o Cl_H_Cl_RCT                 --------------------------------------------------------------------------------
o Cl_H_Cl_TS                    1013.5 
o Cl_H_Cl_RCT: iso @ 1                   1.113e+01   1.129e+00   1.971e+00 
o Cl_H_Cl_TS: iso @ 1            722.0   3.494e+00   1.069e+00   2.765e+00 
                              -------------------------------------------------------------------------------- 
  TST-KIE @ 298.15 K          1.403746    3.185402    1.056428    0.712742    3.366858    2.156937    7.262101 
```

The first column (V-ratio) shows imaginary frequencies for the TS and its isotopomer and their ratio. The next three columns show the terms of the Reduced Isotopic Partition Function Ratios (RPFRs): the zero-point energy differences (ZPE), excitation factors (EXC) and enthalpic Teller–Redlich product factors (TRPF). The product of these four terms defines the KIE value, in this case 3.366858. This assumes classical nuclei (the Born-Oppenheimer approximation). A quantum mechanical tunnelling correction (1D-Tunn) is computed assuming an infinite parabolic barrier in one-dimension, and multiplication gives the corrected KIE (corr-KIE). The values obtained from manually computing the RPFRs are V-ratio = 1.4038; ZPE = 3.1854; EXC = 1.0564; TRPF =	0.7127; KIE=	3.3669 which all agree to 4DP with the above.

The level of theory used for the above calculations is HF/6-31G(d). If a scaling factor is not specified manually, Kinisot tries to match the level/basis set with a database of scaling factors. For the proton transfer example again:  

```python
python Kinisot.py Cl_H_Cl_RCT.out Cl_H_Cl_TS.out -iso "1"
```

To give:

```bash
  Found vibrational scaling factor 0.909 for RHF/6-31G(d) level of theory 
  REF: I. M. Alecu, J. Zheng, Y. Zhao, and D. G. Truhlar, J. Chem. Theory Comput. 6, 2872-2887 (2010). 

                               Temp = 298.15K / Vib. scale factor = 0.909 
                               V-ratio         ZPE         EXC        TRPF         KIE     1D-tunn    corr-KIE 
o Cl_H_Cl_RCT                 --------------------------------------------------------------------------------
o Cl_H_Cl_TS                    1013.5 
o Cl_H_Cl_RCT: iso @ 1                   8.938e+00   1.156e+00   1.971e+00  
o Cl_H_Cl_TS: iso @ 1            722.0   3.118e+00   1.088e+00   2.765e+00
                              -------------------------------------------------------------------------------- 
  TST-KIE @ 298.15 K          1.403746    2.866660    1.062490    0.712742    3.047348    2.156937    6.572937  
```

#### Example 2. Nucleophilic Substitution
The SN2 identity reaction between fluoromethane and fluoride: the output files of reactant and TS contain the results from B3LYP/aug-cc-pVTZ frequency calculations. The temperature defined in Gaussian is unimportant.

Suppose we want to obtain the KIE at a temperature of 273K from deuterating all three H atoms i.e. CD<sub>3</sub>F vs. CH<sub>3</sub>F: these atoms are numbered 2,3 and 4 in both structure files. This is performed with the following command:

```bash
python Kinisot.py CH3F_F_rc.out CH3F_F_ts.out -iso "2 3 4" -t 273
```

To give:

```bash  
  Found vibrational scaling factor 0.985 for RB3LYP/Aug-CC-pVTZ level of theory 
  REF: I. M. Alecu, unpublished (2011). 

                               Temp = 273.0K / Vib. scale factor = 0.985 
                               V-ratio         ZPE         EXC        TRPF         KIE     1D-tunn    corr-KIE 
o CH3F_F_rc                    --------------------------------------------------------------------------------
o CH3F_F_ts                       481.3 
o CH3F_F_rc: iso @ 2 3 4                  6.025e+04   1.253e+00   1.447e+01 
o CH3F_F_ts: iso @ 2 3 4          480.5   7.602e+04   1.112e+00   1.453e+01 
                               -------------------------------------------------------------------------------- 
   TST-KIE @ 273.0 K           1.001663    0.792526    1.126867    0.995606    0.890626    1.001003    0.891519 
                               --------------------------------------------------------------------------------
```

The secondary KIE is 0.8906 (or 0.8916 with tunnelling). 

[![DOI](https://zenodo.org/badge/16266/bobbypaton/Kinisot.svg)](https://zenodo.org/badge/latestdoi/16266/bobbypaton/Kinisot)
---
License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
