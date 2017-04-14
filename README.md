Kinisot
======

A Python program to compute kinetic isotope effects from two Gaussian output files, one of which is a ground state and the other a transition state.

This is a Python version of [Kinisot](http://dx.doi.org/10.5281/zenodo.19272), inspired by the Fortran version originally written by [Henry Rzepa](https://en.wikipedia.org/wiki/Henry_Rzepa). This version does not require parameter files to run and allows easy manipulation of temperature and [vibrational scaling factors](http://t1.chem.umn.edu/freqscale/index.html). The level of theory and basis set are detected from in the output files and the program will attempt to assign the appropriate scaling factor based on data from the [Truhlar group](https://t1.chem.umn.edu/freqscale/index.html). Isotopic substitutions are to be specified by the command line, such that separate computations with Gaussian are not required. The program diagonalizes the mass-weighted Hessian matricies to obtain harmonic vibrational frequencies and Bigeleisen-Mayer Reduced Isotopic Partition Function Ratios. One difference with the Gaussian program itself is that the lowest five/six normal modes (translations and rotations) are not projected out, however, this is also the approach taken by [quiver](https://github.com/ekwan/quiver). Testing  this approach against calculation of the Reduced Isotopic Partition Function Ratios by hand (!) using the frequencies in the Gaussian output files led to agreeement up to 4DP. A one-dimensional tunneling correction is also included, which is the Bell infinite-parabola model.

Also see related discussions on [computing KIE values](http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327)

The current version is currently hard-coded to consider <sup>2</sup>D/<sup>1</sup>H, <sup>13</sup>C/<sup>12</sup>C and <sup>17</sup>O/<sup>16</sup>O isotopic replacements. This can be modified in Hess_to_Freq.py. 

#### Installation
1. Clone the repository https://github.com/bobbypaton/Kinisot.git
2. Add the directory (/pathto/GoodVibes/goodvibes) containing the python files to the PATH environment variable (optional). 
3. Run the script with your Gaussian output files. It has been tested with Python 2 and 3 on OSX

Alternatively `pip install kinisot` will install all classes

**Correct Usage**

```python
Kinisot.py reactant_output ts_output -iso <"atom numbers"> (-t temperature) (-s scalefactor)  
```
*       The two output files contain Gaussian frequency calculations performed for the reactant and transition state at the same level of theory. The temperature is unimportant.
*       The `-iso` flag is required and specifies a string of atom number(s) which are to be substituted for heavier isotopes. Multiple atom numbers require quotation marks and are separated by spaces.
*       The `-t` option specifies temperature (in Kelvin). N.B. This does not have to correspond to the temperature used in the Gaussian calculation since the Reduced Isotopic Partition Function Ratios are evalulated at the requested temperature. The default value is 298.15 K.
*       The `-s` option is a scaling factor for vibrational frequencies. Empirical scaling factors have been determined for several functional/basis set combinations, and these are applied automatically using values from the Truhlar group based on detection of the level of theory and basis set in the output files. The ZPE-scaling factors are selected if available. The default value when no scaling factor is available is 1 (no scale factor). Automated scaling can also be surpressed by `-s 1.0`

#### Example 1. Proton Transfer
Transfer of a proton between two chloride anions; the reactant and TS are both linear. The hydrogen atom is atom number 1 in both files (the numbering needs to be identical). In this example we consider the 2D/1H kie at the default temperature without vibrational scaling.

```python
python Kinisot.py Cl_H_Cl_RCT.out Cl_H_Cl_TS.out -iso "1" -s 1.0
```

The following output is produced:

```python
                                                      Temp = 298.15K / Vib. scale factor = 1.0 
                                                      V-ratio         ZPE         EXC        TRPF         KIE     1D-tunn    corr-KIE 
                                                      -------------------------------------------------------------------------------- 
o kinisot/examples/Cl_H_Cl_RCT                     
o kinisot/examples/Cl_H_Cl_TS                           1013.5 
o kinisot/examples/Cl_H_Cl_RCT: iso @ 1                          1.113e+01   1.129e+00   1.971e+00 
o kinisot/examples/Cl_H_Cl_TS: iso @ 1                   722.0   3.494e+00   1.069e+00   2.765e+00 
                                                      -------------------------------------------------------------------------------- 
   TST-KIE @ 298.15 K                                 1.403746    3.185402    1.056428    0.712742    3.366858    2.156937    7.262101 
                                                      --------------------------------------------------------------------------------
```

The first column (V-ratio) shows imaginary frequencies for the TS and its isotopomer and their ratio. The next three columns show the terms of the Reduced Isotopic Partition Function Ratios: the zero-point energy differences (ZPE), excitation factors (EXC) and enthalpic Teller–Redlich product factors (TRPF). The product of these four terms defines the KIE value, in this case 3.366858. This assumes classical nuclei (the Born-Oppenheimer approximation). A quantum mechanical tunnelling correction (1D-Tunn) is computed assuming an infinite parabolic barrier in one-dimension, and multiplication gives the corrected KIE (corr-KIE). 

#### Example 2. Nucleophilic Substitution
he SN2 identity reaction between fluoromethane and Fluoride. We have previously optimized the structures of GS and TS structures with [Gaussian 09](http://www.gaussian.com): the output files must also contain the results from frequency calculations. The temperature defined in Gaussian is unimportant.

Suppose we want to obtain the KIE at a temperature of 273K from deuterating all three H atoms i.e. CD<sub>3</sub>F vs. CH<sub>3</sub>F: these atoms are numbered 2,3 and 4 in both structure files. This is performed with the following command:

```python
python Kinisot.py CH3F_F_rc.out CH3F_F_ts.out -t 273 -iso "2 3 4"
```

The following output results:


```python
   Temperature = 273.0 K;    Frequency scale factor (default) = 1.0

   Structure                                                       Im Freq         ZPE
o  Examples/CH3F_F_rc: label @ 0                                       N/A    0.039317
o  Examples/CH3F_F_ts: label @ 0                                481.395613    0.038712
o  Examples/CH3F_F_rc: label @ 2 3 4                                   N/A    0.029655
o  Examples/CH3F_F_ts: label @ 2 3 4                            480.596209    0.028845

                            V-ratio   ZPE-ratio         Exc     TR-Prod         KIE        Tunn    corr-KIE
   TST-KIE @ 273.0 K       1.001663    0.789683    1.127764    0.995606    0.888138    1.001003    0.889029
```

Individual terms within the Bigeleisen-Mayer equation are first specifed: the ratio of imaginary frequencies (V-ratio), zero-point energy differences (ZPE-ratio), entropic differences (Exc) and enthalpic Teller–Redlich differences (TR_Prod). The product of these four terms defines the KIE value. This first value assumes classical nuclei (the Born-Oppenheimer approximation). A quantum mechanical tunnelling correction (Tunn) is computed assuming an infinite parabolic barrier in one-dimension, and multiplication by this factor gives the corrected KIE (corr-KIE). This is unlikely to be accurate and should be used with caution.

#### Example 3. Claisen Rearrangement
Consider a \[3,3\]-sigmatropic Claisen rearrangement. Here we compute KIE values separately for various isotopic substitutions and compare with data obtained from [quiver](https://github.com/ekwan/quiver) reported by [Eugene Kwan]() to show that there is no inherent different between the two programs (at least in the first 3DP). The temperature is defined as 393K and the frequencies are scaled by a factor of 0.9614. Again there are two Gaussian output files for the GS and TS.

The first isotopologue is of atom 1 (carbon):

```python
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 1

   Temperature = 393.0 Kelvin    Frequency scale factor = 0.9614

   Structure                                                       Im Freq         ZPE
o  Examples/claisen_gs: label @ 0                                      N/A    0.114417
o  Examples/claisen_ts: label @ 0                               482.835702    0.112582
o  Examples/claisen_gs: label @ 1                                      N/A    0.114245
o  Examples/claisen_ts: label @ 1                               479.060524    0.112415

   V-ratio   ZPE-ratio         Exc     TR-Prod         KIE        Tunn    corr-KIE
   TST-KIE @ 393.0 K      1.007880    1.003956    1.003422    0.997456    1.012747    1.002143    1.014918
```

This was repeated for a series of isotopologues (truncated output):

```python
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 2
   TST-KIE @ 393.0 K      1.000281    0.999252    1.000204    1.002180    1.001916    1.000077    1.001994
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 3
   TST-KIE @ 393.0 K      1.007026    1.022302    1.005865    0.983895    1.018846    1.001913    1.020796
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 4
   TST-KIE @ 393.0 K      1.012716    1.036618    1.001967    0.978989    1.029763    1.003435    1.033300
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 5
   TST-KIE @ 393.0 K      1.000176    0.999077    1.000681    1.001962    1.001896    1.000048    1.001944
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 6
   TST-KIE @ 393.0 K      1.008093    1.004484    1.004208    0.997962    1.014802    1.002201    1.017035
python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso "7 8"
   TST-KIE @ 393.0 K      1.007258    0.885538    1.063932    1.006138    0.954816    1.001976    0.956703
```

Referencing these KIE values with respect to <sup>13</sup>C @ atom number 5 (i.e. this isotopomer is the  internal standard) requires division by this value (1.001896 or 1.001944):

| Isotopomer        | KIE           | corr-KIE  |
|:------------- |:-------------:|:-----:|
| <sup>13</sup>C @ posn 1	| 1.011	| 1.013|
| <sup>13</sup>C @ posn 2 	| 1.000	| 1.000|
| <sup>17</sup>O @ posn 3 	| 1.017	| 1.019|
| <sup>13</sup>C @ posn 4 	| 1.028	| 1.031|
| <sup>13</sup>C @ posn 5 	| 1.000	| 1.000|
| <sup>13</sup>C @ posn 6 	| 1.013	| 1.015|
| <sup>2</sup>H @ posns 7 & 8 	| 0.953	| 0.955|

[![DOI](https://zenodo.org/badge/16266/bobbypaton/Kinisot.svg)](https://zenodo.org/badge/latestdoi/16266/bobbypaton/Kinisot)
---
License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
