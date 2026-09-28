# Wittig reaction of anisaldehyde with a stabilized ylide: ¹³C KIEs

Chen, Nieves-Quiñones, Waas and Singleton (J. Am. Chem. Soc. 2014, 136,
13122) measured ¹³C KIEs for the reaction of anisaldehyde 1 with
Ph₃P=CHCOMe 2 at 67 °C. They are consistent with two transition structures
in series: C–C bond formation (4‡) and P–O bond formation (6‡), with a
commitment factor near 1 (see `case.json` and IMPLEMENTATION_PLAN.md,
Phase 10).

## Jobs to run

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

## What the SI predicts

SI Table S4 lists the M06-2X/6-31+G** KIEs at 67 °C (k12/k13) for each
transition structure; the jobs should reproduce them:

| Position | 4‡ | 6‡ |
| --- | --- | --- |
| ylide CH carbon | 1.022 | 0.994 |
| carbonyl carbon of 1 | 1.043 | 1.015 |
| ketone carbonyl carbon of 2 | 0.998 | 0.999 |
| ketone methyl carbon of 2 | 1.000 | 1.002 |
| aldehyde oxygen (¹⁶O/¹⁸O) | 1.016 | 1.045 |
