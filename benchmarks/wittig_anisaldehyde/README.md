# Wittig reaction of anisaldehyde with a stabilized ylide: ¹³C KIEs

Chen, Nieves-Quiñones, Waas and Singleton (J. Am. Chem. Soc. 2014, 136,
13122) measured ¹³C KIEs for the reaction of anisaldehyde 1 with
Ph₃P=CHCOMe 2 at 67 °C. They are consistent with two transition structures
in series: C–C bond formation (4‡) and P–O bond formation (6‡), with a
commitment factor near 1 (see `case.json` and IMPLEMENTATION_PLAN.md,
Phase 10).

## Result

`series.json` puts 4‡ and 6‡ in series, weighted as the paper's
trajectories are: 128 go on to product for every 76 that return to the
starting materials (C_f = 128/76). It computes each KIE relative to the
experiment's standard, which is a meta carbon for the carbonyl carbon and
the ketone methyl carbon for the ylide carbons:

```
kinisot --job benchmarks/wittig_anisaldehyde/series.json
kinisot --job benchmarks/wittig_anisaldehyde/series_free_energies.json
```

| Position | Trajectories (C_f 1.684) | Free energies (C_f 0.822) | Measured |
| --- | --- | --- | --- |
| carbonyl carbon of 1 | 1.0326 | 1.0277 | 1.033 ± 0.002, 1.032 ± 0.002 |
| ylide CH carbon | 1.0111 | 1.0059 | 1.011 ± 0.002, 1.011 ± 0.004 |
| ylide ketone carbonyl carbon | 0.9979 | 0.9975 | 0.997 ± 0.002 |

The trajectory weighting reproduces all three measurements (mean absolute
deviation 0.0004); the free-energy weighting misses the two carbons that
change bonding by 0.005. That is the paper's conclusion. Most trajectories
pass the betaine without equilibrating, so the free energies of the two
transition structures do not decide the outcome.

- **Against the paper.** The absolute KIEs compare with Table 1 as
  follows:
  - trajectories: 1.0325 and 1.0118 (Table 1: 1.033 and 1.012);
  - free energies: 1.0276 and 1.0069 (Table 1: 1.028 and 1.008). The
    paper's own single-structure values give 1.0070 for the ylide CH
    carbon, so the printed 1.008 reflects the rounding of its inputs.
- **The free-energy C_f.** Kinisot weights each step by its rate
  including Bell tunnelling. C_f is therefore 0.822, not the
  exp[(G₄ − G₆)/RT] = 0.862 of the free energies alone. With 0.862 the
  KIEs would be 3 × 10⁻⁴ higher.
- **The meta standard.** The two meta carbons are equivalent in the
  experiment. Relative to the other one (anisaldehyde 2, TS 27), the
  carbonyl carbon's KIE is 5 × 10⁻⁴ lower.

The benchmark report lists the trajectory weighting as this case and the
free-energy weighting as `wittig_anisaldehyde_statistical`, an
alternative left out of the overall mean.

## Jobs

Gaussian frequency jobs at the SI geometries, at the paper's level:
M06-2X/6-31+G(d,p) with PCM for THF and default radii. The paper used
Gaussian 09, whose default grid is FineGrid, so the inputs ask for
`int=finegrid`.

| Input | Structure | Atoms | SI: E, ZPE (hartree) |
| --- | --- | --- | --- |
| `anisaldehyde_1.gjf` | anisaldehyde 1 | 18 | −459.931869643, 0.143666 |
| `ylide_2.gjf` | Ph₃P=CHCOMe 2 | 42 | −1227.87396966, 0.340415 |
| `ts_4.gjf` | 4‡, C–C bond formation ("TS1" in the SI) | 60 | −1687.79361546, 0.485845 |
| `ts_6.gjf` | 6‡, P–O bond formation ("TS1.5" in the SI) | 60 | −1687.79472605, 0.487529 |

Each title line carries these values. A job at the right geometry and level
reproduces its energy closely, and each transition structure has exactly
one imaginary mode. Gaussian 16 may differ from Gaussian 09 in its PCM or
integral defaults by a little; a difference of millihartrees points to a
different solvent model or radii.

The betaine 5 between the two transition structures is not needed: each
step's KIE runs from the starting materials 1 + 2.

## Status (2026-09-28)

Gaussian 16 (C.01) jobs for all four structures are in this directory.

| Log | E − E(SI) (hartree) | ZPE (SI) | Imaginary modes |
| --- | --- | --- | --- |
| `anisaldehyde_1.log` | −4 × 10⁻⁶ | 0.143664 (0.143666) | none |
| `ylide_2.log` | −2.5 × 10⁻⁵ | 0.340401 (0.340415) | 9i cm⁻¹, a phenyl torsion below the 50i cutoff |
| `ts_4.log` | −4.6 × 10⁻⁵ | 0.485865 (0.485845) | 285i cm⁻¹ |
| `ts_6.log` | −4.6 × 10⁻⁵ | 0.487545 (0.487529) | 111i cm⁻¹ |

The SI geometries are not exactly stationary under Gaussian 16, so the case
projects out translations and rotations (`project: true`). Without that,
low modes of the ylide and 4‡ are up to 42 cm⁻¹ from the printed ones;
the KIEs change by less than 10⁻⁶ either way.

(The first `ts_6.gjf` held the coordinates of the oxaphosphetane
`OP1direct`, not of 6‡. The SI prints 6‡ ("TS1.5") in a different
coordinate format from the other structures, and the extraction took the
next block. The corrected input has P···O 2.53 Å and C–C 1.59 Å, and its
job reproduces the SI.)

**Both transition structures reproduce SI Table S4** at 67 °C with a
scaling factor of 0.9614 and Bell tunnelling, as the SI states, within
5 × 10⁻⁴ at every position (`tests/test_wittig.py`):

| Position | anisaldehyde 1 / ylide 2 | 4‡ and 6‡ | Kinisot 4‡ (SI) | Kinisot 6‡ (SI) |
| --- | --- | --- | --- | --- |
| ylide CH carbon | ylide 1 | 4 | 1.0220 (1.022) | 0.9945 (0.994) |
| carbonyl carbon (CHO) | anisaldehyde 9 | 2 | 1.0427 (1.043) | 1.0152 (1.015) |
| ipso | anisaldehyde 4 | 5 | 0.9995 (0.999) | 0.9992 (0.999) |
| ortho | anisaldehyde 5, 3 | 24, 28 | 1.0011, 1.0007 (1.001, 1.001) | 1.0014, 1.0012 (1.001, 1.001) |
| meta | anisaldehyde 6, 2 | 25, 27 | 0.9998, 1.0002 (1.000, 1.000) | 1.0000, 1.0005 (1.000, 1.001) |
| para | anisaldehyde 1 | 26 | 1.0004 (1.000) | 1.0009 (1.001) |
| ketone carbonyl carbon | ylide 37 | 51 | 0.9985 (0.998) | 0.9987 (0.999) |
| ketone methyl carbon | ylide 38 | 52 | 0.9997 (1.000) | 1.0022 (1.002) |
| methoxy carbon | anisaldehyde 13 | 57 | 0.9995 (1.000) | 0.9997 (1.000) |
| aldehyde oxygen (¹⁸O) | anisaldehyde 17 | 1 | 1.0159 (1.016) | 1.0448 (1.045) |

Atom numbers count from 1 in each file. The mapping was found by matching
bond graphs and then geometries; 4‡ and 6‡ list their atoms in the same
order. For example, the ylide CH carbon is `--iso 0 --iso 1 --iso 4`
(anisaldehyde, ylide, transition structure).
