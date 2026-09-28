# Phosphetane-catalysed deoxygenation of nitrobenzene: ¹⁸O KIEs

Kang and Radosevich (Tetrahedron 2025, 186, 134892) measured ¹⁸O KIEs for
the P(III)/P(V)=O-catalysed reductive N-arylation of nitrobenzene at 120 °C:
1.033 ± 0.003 for one ¹⁸O and 1.066 ± 0.003 for two, both relative to
[¹⁵N]-nitrobenzene. Their SI gives the M06-2X/6-31G(d,p) structures but no
frequencies. The ORCA 6.1.0 frequency jobs here are at those geometries.

| Files | Structure |
| --- | --- |
| `phno2.inp`, `.out`, `.hess` | nitrobenzene (no imaginary mode) |
| `ts2.inp`, `.out`, `.hess` | TS2, the concerted [3+1] cheletropic transition structure the paper assigns (326i cm⁻¹) |
| `ts1b.inp`, `.out`, `.hess` | TS1B, the monotopic frontside transition structure the paper rejects (365i cm⁻¹) |

The SI geometries are not exactly stationary at this ORCA version, so the
cases project out translation and rotation, as ORCA does for the frequencies
it prints. Kinisot then reproduces every printed frequency, and the KIEs move
by less than 4 × 10⁻⁶.

## Results

At 393 K with the paper's scale factor (0.9614), uncorrected:

| | Oxygen A | Oxygen B | One ¹⁸O (average) | Both ¹⁸O |
| --- | --- | --- | --- | --- |
| TS2 (O12, O13) | 1.0320 | 1.0325 | 1.0322 | 1.0661 |
| TS1B (O2 attacked, O3 spectator) | 1.0474 | 1.0140 | 1.0304 | 1.0621 |
| Paper, PyQuiver | | | TS2 1.0321, TS1B 1.0468 | TS2 1.0657, TS1B 1.0612 |
| Experiment | | | 1.033 ± 0.003 | 1.066 ± 0.003 |

With Bell tunnelling, TS2 gives 1.0346 and 1.0709, and TS1B 1.0319 and
1.0653. The singly labelled values here are the exact (harmonic) averages;
the benchmark runner averages the two positions geometrically and reports
1.0305 and, with Bell, 1.0321 for TS1B (the same to 10⁻⁵ for TS2).

- **The jobs reproduce the paper's predictions**, to within 0.0006 for
  TS2 and 0.0010 for TS1B.
- **The paper's TS1B value for one ¹⁸O is the attacked oxygen alone.**
  Nitrobenzene's two oxygens are equivalent, so a singly labelled molecule
  reacts through both TS1B isotopomers. Its KIE is then
  2/(1/KIE_attacked + 1/KIE_spectator) = 1.0304, within the measured value.
- **The KIEs do not separate the two mechanisms.** TS1B's doubly labelled
  KIE (1.0621) is close to the square of its singly labelled one, so it
  scales multiplicatively as TS2 does. Without tunnelling TS2 fits better;
  with Bell's correction TS1B does. TS1B is excluded by its energy
  (50.4 against 33.4 kcal/mol), not by the ¹⁸O KIEs.

`tests/test_nitroarene.py` checks these numbers.
