# Claisen rearrangement of allyl vinyl ether

**Question.** Which atoms change bonding in the concerted [3,3] transition
structure? Heavy-atom KIEs are largest at the atoms whose bonds break or
form, so a position-by-position ¹³C scan maps the transition structure.
This is the reaction of the combined experimental and computational study by
Meyer, DelMonte and Singleton (J. Am. Chem. Soc. 1999, 121, 10865,
[doi:10.1021/ja992372h](https://doi.org/10.1021/ja992372h)), who measured
¹³C, ²H and ¹⁷O KIEs by NMR at natural abundance and compared them with
Bigeleisen–Mayer predictions including a one-dimensional tunnelling
correction, which is the model Kinisot implements.

**Files.** `claisen_gs.out` (allyl vinyl ether) and `claisen_ts.out`
(chair transition structure), B3LYP/6-31G(d). Atoms 1–6 are the carbons and
the ether oxygen (C1, C2, O3, C4, C5, C6), atoms 7–14 the hydrogens; the
numbering is the same in both files, so one `--iso` per run is enough.

**Commands** (393 K, ZPE scaling factor 0.961). A bare atom number means
the default heavy isotope (¹³C, ²H, ¹⁸O); `3:17O` asks for oxygen-17
explicitly, which is what the NMR experiment measured:

```
cd tests/data/gaussian
for atoms in 1 2 3:18O 3:17O 4 5 6 7,8; do
    kinisot --rct claisen_gs.out --ts claisen_ts.out --iso $atoms -t 393 -s 0.961
done
```

**Result lines** (from `expected_output.dat`):

```
                          V-ratio        ZPE        EXC       TRPF        KIE    1D-tunn   corr-KIE
C1     KIE @ 393.0 K     1.007880   1.003953   1.003421   0.997456   1.012744   1.001970   1.014739
C2     KIE @ 393.0 K     1.000281   0.999252   1.000203   1.002180   1.001916   1.000071   1.001988
18O3   KIE @ 393.0 K     1.013529   1.042882   1.011560   0.969022   1.036088   1.003355   1.039564
17O3   KIE @ 393.0 K     1.007027   1.022288   1.005866   0.983895   1.018834   1.001759   1.020625
C4     KIE @ 393.0 K     1.012716   1.036594   1.001970   0.978989   1.029742   1.003157   1.032993
C5     KIE @ 393.0 K     1.000176   0.999077   1.000680   1.001962   1.001895   1.000044   1.001940
C6     KIE @ 393.0 K     1.008093   1.004481   1.004206   0.997962   1.014797   1.002022   1.016849
H7,H8  KIE @ 393.0 K     1.007258   0.885606   1.063924   1.006138   0.954882   1.001816   0.956616
```

**Reading it.** The ¹³C KIEs are largest at C4 (1.033) and at C1 and C6
(1.015, 1.017): the C4–O3 bond breaks and the C1–C6 bond forms in the
transition structure, while C2 and C5 (1.002) are spectators. The oxygen
KIE at O3 (1.040 for ¹⁸O, 1.021 for ¹⁷O, roughly the 2:1 ratio expected
from the mass differences) reflects the breaking C–O bond. Substituting both hydrogens 7 and
8 gives an inverse secondary ²H₂ KIE (0.957), the signature of a CH₂ group
that rehybridizes from sp² towards sp³ as a new bond forms. The
tunnelling correction (`1D-tunn`) changes heavy-atom KIEs by 0.2–0.3 %,
the size of the experimental uncertainty in the NMR measurements.

**More options, same C4 KIE.** `--project` removes translations and
rotations by Eckart projection instead of discarding the six lowest modes
(the KIE changes by 5 × 10⁻⁸ on this well-converged geometry);
`--reference 5` also reports the KIE relative to the C5 position, which is
how natural-abundance NMR experiments are referenced; `--tunneling skodje`
uses the Skodje–Truhlar correction with the barrier read from the energies
in the files (28.9 kcal/mol here, so it coincides with Bell); and a
temperature list or range gives one result line per temperature:

```
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393 -s 0.961 --project
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393 -s 0.961 --reference 5
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 393 -s 0.961 --tunneling skodje
kinisot --rct claisen_gs.out --ts claisen_ts.out --iso 4 -t 300:400:50 -s 0.961
```

```
  relative to iso @ 5 / 5:                                                    KIE = 1.001895   1.027795              1.030994
  KIE @ 300.0 K                                    1.012716   1.048208   1.001161   0.978989   1.040439   1.005638   1.046304
```
