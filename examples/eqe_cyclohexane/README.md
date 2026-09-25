# CD₃ axial/equatorial preference in 1,1,3,3-tetramethylcyclohexane

**Question.** Does a CD₃ group prefer the axial or the equatorial position
relative to CH₃? A geminal dimethyl group has one methyl axial and one
equatorial; swapping which one carries the deuterium is a conformational
equilibrium whose constant is an equilibrium isotope effect (EQE, also
called EIE). The same system is used as an EIE demonstration by Grazioli
et al. for the PyQuiverHS web tool (J. Phys. Org. Chem. 2026, 39, e70099,
[doi:10.1002/poc.70099](https://doi.org/10.1002/poc.70099)).

**Files.** One frequency calculation, `tetramethylcyclohexane.out`
(B3LYP/6-31G(d)), used as both reactant and product: the two isotopologues
differ only in which methyl group carries the three deuteriums (hydrogens
24–26 or 28–30).

**Commands** (unscaled frequencies, `-s 1`):

```
cd tests/data/gaussian
kinisot --rct tetramethylcyclohexane.out --prd tetramethylcyclohexane.out --iso 24,25,26 --iso 28,29,30 -s 1 -t 290
kinisot --rct tetramethylcyclohexane.out --prd tetramethylcyclohexane.out --iso 24,25,26 --iso 28,29,30 -s 1 -t 300
```

**Result lines** (from `expected_output.dat`):

```
                                  ZPE        EXC       TRPF        EQE    1D-tunn   corr-EQE
EQE @ 290.0 K                1.027782   1.025381   0.985373   1.038454   1.000000   1.038454
EQE @ 300.0 K                1.026844   1.025140   0.985373   1.037262   1.000000   1.037262
```

**Reading it.** There is no reaction coordinate, so the V-ratio column is
empty and the tunnelling correction is 1. The EQE is the equilibrium
constant for moving the CD₃ label from the first methyl group (hydrogens
24–26) to the second (28–30); a value of 1.038 at 290 K means the second
site is favoured by 3.8 % for CD₃ relative to CH₃. The effect comes almost
entirely from zero-point energy (ZPE 1.028) and the low-frequency
excitation term (EXC 1.025) and shrinks with temperature, as conformational
isotope effects do. The TRPF term is temperature independent.
