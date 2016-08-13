Kinisot
======

A Python program to compute kinetic isotope effects from two Gaussian output files (a ground state and transition state)

This is a Python version of [Kinisot](http://dx.doi.org/10.5281/zenodo.19272), from the Fortran version originally written by Henry Rzepa. This version gives identical numerical results as the original, although hopefully with greater flexibility in terms of specifying temperature, [vibrational scaling factors](http://t1.chem.umn.edu/freqscale/index.html). Isotopic substitutions can be specified when using the program, such that separate computations with Gaussian are not required. The program diagonalizes the mass-weighted Hessian matricies to obtain the vibrational frequencies and partition function. One difference with the Gaussian program itself is that the lowest six normal modes are not projected out, however, this is also the approach taken by [quiver](https://github.com/ekwan/quiver) which does not seem to have disastrous effects. Numerical results agree to 3DP with those obtained from quiver, which is not unsurprising since the underlying approach is the same. The tunneling correction is an infinite-parabola model.

Also see related discussions on [computing KIE values](http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327)


The current version is currently hard-coded to only consider <sup>2</sup>D/<sup>1</sup>H, <sup>13</sup>C/<sup>12</sup>C and <sup>17</sup>O/<sup>16</sup>O isotopic replacements. This can be modified in Hess_to_Freq.py

Example 1. Nucleophilic Substitution
------

Consider the SN2 identity reaction between fluoromethane and Fluoride. We have previously optimized the structures of GS and TS structures with [Gaussian 09](http://www.gaussian.com): the output files must also contain the results from frequency calculations. The temperature defined in Gaussian is unimportant.

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

Example 2. Claisen Rearrangement
------

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
