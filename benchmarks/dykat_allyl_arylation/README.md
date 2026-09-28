# Rh-catalysed arylation of racemic 3-chlorocyclohexene: ¹³C KIEs

Van Dijk et al. (Nat. Catal. 2021, 4, 284) measured ¹³C KIEs on the product,
3-phenylcyclohexene, at 40 °C. Both enantiomers of the allyl chloride 3 react
through their own anti-oxidative-addition transition structure and converge
on one η³-allyl intermediate, so each product carbon comes from a different
substrate carbon in each enantiomer's channel (see `case.json`).

## Result

`channels.json` computes each channel from its Rh-bound reactant complex
to its anti-oxidative-addition transition structure, as the paper did, at
313.15 K with Bell tunnelling, relative to C5, and combines them with the
paper's selectivity s = k(S)/k(R) = 3.3:

```
kinisot --job benchmarks/dykat_allyl_arylation/channels.json
```

| Product carbon | (R)-3 channel | (S)-3 channel | Combined | Figure 3d | Measured |
| --- | --- | --- | --- | --- | --- |
| C1 (bears Ph) | 1.0290 (C–Cl) | 0.9986 (alkene CH) | 1.0055 | 1.006 | 1.0097 ± 0.0044 |
| C2 | 0.9986 (central) | 0.9945 (central) | 0.9955 | 0.995 | 0.9952 ± 0.0016 |
| C3 | 0.9996 (alkene CH) | 1.0348 (C–Cl) | 1.0264 | 1.027 | 1.0244 ± 0.0027 |
| C4 | 1.0020 (CH₂ next to CH) | 1.0063 (CH₂ next to C–Cl) | 1.0053 | 1.005 | 1.0084 ± 0.0040 |
| C6 | 1.0058 (CH₂ next to C–Cl) | 1.0028 (CH₂ next to CH) | 1.0035 | 1.003 | 0.9983 ± 0.0048 |

Every combined KIE lies within 1.1 standard errors of the measurement
(mean absolute deviation 0.0030). The per-channel values reproduce SI
Tables 24 and 25, and the combined ones Figure 3d, within 6 × 10⁻⁴
(`tests/test_dykat.py`). The paper's three decimals and the reactant
complexes' coordinates, printed to three decimals in the SI, account for
the differences.

- **Unscaled frequencies.** The SI does not give a scaling factor. The
  Kinisot of 2021 applied 1.0 to any level of theory missing from its
  table, and neither ωB97X-D/6-31G(d) nor a `genecp` basis is in it; 1.0
  reproduces the tables. The
  Truhlar factor for ωB97X-D/6-31G(d), 0.975, moves the C–Cl KIEs up to
  9 × 10⁻⁴ away from them.
- **Atom numbers.** Each product carbon has the same atom number in both
  channels (C1 is atom 118, C2 117, C3 122, C4 121, C5 120, C6 119), so the
  `Labelled atoms` lines of the output show the substrate position it
  comes from in each.

## From free 3 instead of the reactant complex

In an intermolecular competition, free 3 is the reactant in the flask. If
its binding to Rh were a fast equilibrium ahead of oxidative addition, the
measured KIE would include the equilibrium isotope effect of that binding.
`channels_free3.json` computes the same channels from free 3:

| Product carbon | From the complex (paper) | From free 3 | Binding EIE, (R)-3 / (S)-3 | Measured |
| --- | --- | --- | --- | --- |
| C1 | 1.0055 | 1.0147 | 1.0000 / 1.0118 | 1.0097 ± 0.0044 |
| C2 | 0.9955 | 1.0060 | 1.0077 / 1.0114 | 0.9952 ± 0.0016 |
| C3 | 1.0264 | 1.0292 | 1.0105 / 1.0004 | 1.0244 ± 0.0027 |
| C4 | 1.0053 | 1.0041 | 0.9996 / 0.9987 | 1.0084 ± 0.0040 |
| C6 | 1.0035 | 1.0038 | 0.9987 / 1.0008 | 0.9983 ± 0.0048 |

- **The binding EIE** (free 3 → complex, per substrate position) is 1.008
  to 1.012 at the alkene carbons that bind to Rh, and near 1 elsewhere.
- **The measurement does not show it.** From free 3, C2 is 1.0060 against
  the measured 0.9952 ± 0.0016, and the mean absolute deviation doubles
  (0.0061 against 0.0030). The paper's model, from the complex, is the one
  the data support.
- **The d functions do not matter here.** Free 3 with Cartesian d functions
  (`allyl_chloride_3.log`) gives the same KIEs within 1 × 10⁻⁴.

The benchmark report lists this as `dykat_allyl_arylation_free3`, an
alternative left out of the overall mean.

## Jobs

Gaussian frequency jobs at the paper's geometries (Supplementary Data 1),
at its level: ωB97X-D/6-31G(d) with LANL2DZ plus an f function (exponent
1.35) on Rh, gas phase. The paper used Gaussian 09 D.01, whose default grid
is FineGrid, so the inputs ask for `int=finegrid`.

| Input | Structure | Atoms | SI Table 20: E, ZPE (hartree) |
| --- | --- | --- | --- |
| `allyl_chloride_3.gjf` | 3, pseudo-axial | 16 | −694.160144, 0.139852 |
| `allyl_chloride_3_5d.gjf` | the same, spherical d functions | 16 | none (for a consistent KIE from free 3) |
| `r3_reactant_complex.gjf` | (R)-3 bound to Rh–Ph | 132 | −3910.736886, 1.104755 |
| `s3_reactant_complex.gjf` | (S)-3 bound to Rh–Ph | 132 | −3910.739063, 1.105354 |
| `r3_anti_oa_ts.gjf` | (R)-3 anti-oxidative addition | 132 | −3910.712464, 1.105016 |
| `s3_anti_oa_ts.gjf` | (S)-3 anti-oxidative addition | 132 | −3910.702421, 1.104971 |

Each title line carries these values. The Gaussian 16 C.01 jobs (the
`.log` files) gave:

| Job | E (hartree) | ZPE | Imaginary mode |
| --- | --- | --- | --- |
| `allyl_chloride_3.log` | −694.160144 | 0.139852 | none |
| `allyl_chloride_3_5d.log` | −694.154046 | 0.139856 | none |
| `r3_reactant_complex.log` | −3910.736867 | 1.104806 | none |
| `s3_reactant_complex.log` | −3910.739041 | 1.105368 | none |
| `r3_anti_oa_ts.log` | −3910.712464 | 1.105016 | 176i |
| `s3_anti_oa_ts.log` | −3910.702421 | 1.104971 | 183i |

Free 3 and both transition structures reproduce Table 20 to the last
printed digit, which settles the (R)/(S) labels and the d functions below.
The reactant complexes are 2 × 10⁻⁵ hartree off, with residual forces of
10⁻³ hartree/bohr, because the SI prints their coordinates to three
decimals (the transition structures have five). The spherical-d job for
free 3 has no SI value. It lies 6 millihartree above the Cartesian one, as
a different basis should, with residual forces of 4 × 10⁻⁴ because the
geometry was optimized with Cartesian d functions.

- **Basis functions.** The four Rh inputs give 6-31G(d) through general
  basis input (`genecp`) without `5D` or `6D`, so Gaussian uses its default
  for general basis input: spherical d functions (1107 basis functions
  rather than 1179). The `6-31G(d)` keyword alone uses Cartesian ones, as
  in `allyl_chloride_3.gjf`. Both reproduce Table 20, so the paper mixed
  the two conventions.
- **Free 3.** For KIEs from free 3 with the same functions as the
  transition structures, `allyl_chloride_3_5d.gjf` repeats the job with
  spherical d functions. It turned out not to matter (see above).
