# Phosphetane-catalysed deoxygenation of nitrobenzene: ¹⁸O KIEs

Kang and Radosevich (Tetrahedron 2025, 186, 134892) measured ¹⁸O KIEs for
the P(III)/P(V)=O-catalysed reductive N-arylation of nitrobenzene at 120 °C:
1.033 ± 0.003 for one ¹⁸O and 1.066 ± 0.003 for two, both relative to
[¹⁵N]-nitrobenzene. Their SI gives the M06-2X/6-31G(d,p) structures, but no
frequencies, so the three frequency jobs here are still to be run.

## Jobs to run

| Input | Structure | Atoms |
| --- | --- | --- |
| `phno2.inp` | nitrobenzene | 14 |
| `ts2.inp` | TS2, the concerted [3+1] cheletropic transition structure (the paper's assignment) | 43 |
| `ts1b.inp` | TS1B, the monotopic frontside transition structure (the alternative) | 43 |

Each is an ORCA frequency job at the SI geometry, at the authors' level
(M06-2X/6-31G(d,p), gas phase; they used ORCA 6.0.1). Run them with
`orca phno2.inp > phno2.out` and so on, and keep each `.hess` file next to
its `.out`. If your ORCA version has no analytic Hessian for M06-2X, replace
`Freq` with `NumFreq`. A job at the right geometry has no imaginary mode for
nitrobenzene and exactly one for each transition structure.

Then set, in `case.json`, `"reactants": ["phno2.out"]` and
`"transition_structure": ["ts2.out"]`; in `../nitroarene_phosphetane_ts1b/case.json`
the paths are `../nitroarene_phosphetane/phno2.out` and
`../nitroarene_phosphetane/ts1b.out`. The isotope labels are already in both
files.

## Why two cases

In TS2 phosphorus bonds to both oxygens almost equally (P–O 2.271 and
2.261 Å). In TS1B it attacks one (2.223 Å) and not the other (2.727 Å).
Nitrobenzene's two oxygens are equivalent, so a singly labelled molecule
reacts through both TS1B isotopomers. Its KIE is then
2/(1/KIE_attacked + 1/KIE_spectator), not the KIE of either position. The
paper's TS1B prediction for one ¹⁸O (1.0468) may be the attacked position
alone, and the two cases will show how much the average changes it.
