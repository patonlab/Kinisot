# Diels–Alder reaction of isoprene with maleic anhydride

**Question.** Is the cycloaddition synchronous? ¹³C KIEs at the two diene
termini (and the two dienophile carbons) answer this: equal KIEs mean both
bonds form to the same extent in the transition structure, unequal KIEs an
asynchronous one. This is the reaction with which Singleton and Thomas
introduced the natural-abundance NMR method for measuring many small KIEs at
once (J. Am. Chem. Soc. 1995, 117, 9357,
[doi:10.1021/ja00141a030](https://doi.org/10.1021/ja00141a030)).

**Files.** `diene.out` (isoprene), `dienophile.out` (maleic anhydride),
`DATS_rct.out` (the pre-reaction complex of the two, with the transition
structure's atom numbering) and `DATS.out` (transition structure), all
B3LYP/6-31G(d).

**Commands** (298.15 K, ZPE scaling factor 0.963). Kinisot accepts the
reactants either as one file (the complex) or as two separate files. With
two files, give one `--iso` per file in the order of the `--rct` flags and
then the TS, using `0` for a file without a substituted atom. Atom 6 of the
diene is atom 15 of the transition structure:

```
cd tests/data/gaussian
kinisot --rct DATS_rct.out --ts DATS.out --iso 15 -s 0.963
kinisot --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 6 --iso 15 -s 0.963
```

and the other positions (see `../run_examples.sh` for the full list):

```
kinisot --rct dienophile.out --rct diene.out --ts DATS.out --iso 0 --iso 10 --iso 19 -s 0.963
kinisot --rct dienophile.out --rct diene.out --ts DATS.out --iso 1 --iso 0 --iso 1 -s 0.963
```

**Result lines** (from `expected_output.dat`):

```
                                          V-ratio        ZPE        EXC       TRPF        KIE    1D-tunn   corr-KIE
complex, TS C15        KIE @ 298.15 K   1.009798   1.003033   1.012514   0.992401   1.017743   1.003711   1.021520
diene C6 -> TS C15     KIE @ 298.15 K   1.009798   1.001688   0.986484   1.020846   1.018631   1.003711   1.022412
diene C10 -> TS C19    KIE @ 298.15 K   1.000073   0.995775   0.981489   1.023885   1.000760   1.000028   1.000787
diene C1 -> TS C10     KIE @ 298.15 K   1.000449   0.999119   0.996041   1.008008   1.003583   1.000172   1.003756
diene C2 -> TS C11     KIE @ 298.15 K   1.000209   0.994482   0.992761   1.014544   1.001851   1.000080   1.001931
diene C4 -> TS C13     KIE @ 298.15 K   1.007607   1.000009   0.986061   1.021625   1.015057   1.002890   1.017990
dienophile C1 -> TS C1 KIE @ 298.15 K   1.008974   1.004699   0.996274   1.008113   1.018131   1.003403   1.021597
dienophile C1 -> TS C2 KIE @ 298.15 K   1.010789   1.007034   0.996800   1.006756   1.021496   1.004081   1.025665
dienophile H5 -> TS H7 KIE @ 298.15 K   1.000107   0.995021   0.996275   1.006985   0.998345   1.000041   0.998385
dienophile H5 -> TS H5 KIE @ 298.15 K   1.000097   0.997604   0.996466   1.006177   1.000316   1.000037   1.000353
```

**Reading it.** The first two lines are the same ¹³C KIE computed from the
pre-reaction complex (one reactant file) and from the separated reactants
(two files): the individual ZPE/EXC/TRPF factors differ because the complex
and the free molecules have different low-frequency modes, but the KIE
agrees to 0.001. Kinisot multiplies the partition-function terms of all
files on a side, so a bimolecular reaction needs no complex. The two diene
termini (TS atoms 15 and 19) show KIEs of 1.022 and 1.001: the two new
bonds are formed to very different extents, i.e. the transition structure
is markedly asynchronous, in line with the experimental picture.
