Kinisot
------

This is a faithful translation of Kinisot, written by Henry Rzepa in Fortran. This Python version is intended to give identical results, although with greater flexibility in terms of specifying temperature, vibrational scaling factors - and does not require the isotopomers to be calculated separately in e.g. Gaussian. The tunneling correction is an infinite-parabola model. The numerical results seem identical with quiver - https://github.com/ekwan/quiver - which is not unsurprising since the underlying euqations are the same. 

The original Fortran is available: Rzepa, Henry S.., "KINISOT. A basic program to calculate kinetic isotope effects using normal coordinate analysis of transition state and reactants.", 2015. http://dx.doi.org/10.5281/zenodo.19272

Also see related discussions on http://www.ch.imperial.ac.uk/rzepa/blog/?p=14327


This version is currently hard-coded to only consider 2D/1H, 13C/12C and 17O/16O replacements - in an inelegant fashion. This can be modified in Hess_to_Freq.py


Example 1: SN2 identity reaction of MeF with Fluoride, comparing CD3 with CH3 (these atoms are numbered 2,3 and 4 in both structure files). The temperature is set to 273K, and there is no scaling of force constants. Presently it is assumed that the isotopomer is the next most abundant (i.e. C13, H2, O17).

$ python Kinisot.py Kinisot_example/*out -t 273 - iso "2 3 4"

   Temperature = 273.0 Kelvin    Frequency scale factor (default) = 1.0

   Structure                                                       Im Freq         ZPE
o  Examples/CH3F_F_rc: label @ 0                                       N/A    0.039317
o  Examples/CH3F_F_ts: label @ 0                                481.395613    0.038712
o  Examples/CH3F_F_rc: label @ 2 3 4                                   N/A    0.029655
o  Examples/CH3F_F_ts: label @ 2 3 4                            480.596209    0.028845

                            V-ratio   ZPE-ratio         Exc     TR-Prod         KIE        Tunn    corr-KIE
   TST-KIE @ 273.0 K       1.001663    0.789683    1.127764    0.995606    0.888138    1.001003    0.889029


Example 2: Claisen rearrangement (comparison with data from quiver for different isotopologues). The temperature is set to 393K and the frequencies are scaled by 0.9614. The first five calculations consider a single atom substitution by its heavier isotope. The last considers substitution of two H-atoms (7 & 8) with deuterium.

$ python Kinisot.py Examples/claisen_*out -t 393 -s 0.9614 -iso 1

   Temperature = 393.0 Kelvin    Frequency scale factor = 0.9614

   Structure                                                       Im Freq         ZPE
o  Examples/claisen_gs: label @ 0                                      N/A    0.114417
o  Examples/claisen_ts: label @ 0                               482.835702    0.112582
o  Examples/claisen_gs: label @ 1                                      N/A    0.114245
o  Examples/claisen_ts: label @ 1                               479.060524    0.112415

                           V-ratio   ZPE-ratio         Exc     TR-Prod         KIE        Tunn    corr-KIE
   TST-KIE @ 393.0 K      1.007880    1.003956    1.003422    0.997456    1.012747    1.002143    1.014918

   ... -iso 2
   TST-KIE @ 393.0 K      1.000281    0.999252    1.000204    1.002180    1.001916    1.000077    1.001994
   ... -iso 3
   TST-KIE @ 393.0 K      1.007026    1.022302    1.005865    0.983895    1.018846    1.001913    1.020796
   ... -iso 4
   TST-KIE @ 393.0 K      1.012716    1.036618    1.001967    0.978989    1.029763    1.003435    1.033300
   ... -iso 5
   TST-KIE @ 393.0 K      1.000176    0.999077    1.000681    1.001962    1.001896    1.000048    1.001944
   ... -iso 6
   TST-KIE @ 393.0 K      1.008093    1.004484    1.004208    0.997962    1.014802    1.002201    1.017035
   ... -iso "7 8"
   TST-KIE @ 393.0 K      1.007258    0.885538    1.063932    1.006138    0.954816    1.001976    0.956703

	Classical and tunneling-corrected KIEs referenced with respect to 13C @ atom number 5 (i.e. divide through by this value):
	13C @ posn 1:	1.011	1.013
	13C @ posn 2: 	1.000	1.000
	17O @ posn 3: 	1.017	1.019
	13C @ posn 4: 	1.028	1.031
	13C @ posn 5: 	1.000	1.000
	13C @ posn 6: 	1.013	1.015
	2H @ posn 7&8: 	0.953	0.955
	These values appear to be indistinguishable (at least up to 3DP) to those produced by Quiver (https://github.com/ekwan/quiver) 
---
License: CC-BY
