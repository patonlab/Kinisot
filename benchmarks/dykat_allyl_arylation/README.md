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
| `r3_reactant_complex.log` | −3910.736867 | 1.104806 | none |
| `s3_reactant_complex.log` | −3910.739041 | 1.105368 | none |
| `r3_anti_oa_ts.log` | −3910.712464 | 1.105016 | 176i |
| `s3_anti_oa_ts.log` | −3910.702421 | 1.104971 | 183i |

Free 3 and both transition structures reproduce Table 20 to the last
printed digit, which
settles the (R)/(S) labels and the d functions below. The reactant
complexes are 2 × 10⁻⁵ hartree off, with residual forces of 10⁻³
hartree/bohr, because the SI prints their coordinates to three decimals
(the transition structures have five). The job for `allyl_chloride_3_5d.gjf`
is still to be run.

- **Basis functions.** The four Rh inputs give 6-31G(d) through general
  basis input (`genecp`) without `5D` or `6D`, so Gaussian uses its default
  for general basis input: spherical d functions (1107 basis functions
  rather than 1179). The `6-31G(d)` keyword alone uses Cartesian ones, as
  in `allyl_chloride_3.gjf`. Both reproduce Table 20, so the paper mixed
  the two conventions.
- **Free 3.** A KIE from free 3 needs the same functions as the transition
  structures, so it needs the spherical-d job, `allyl_chloride_3_5d.gjf`.
  The paper computed each channel from the reactant complex. For an
  intermolecular KIE measured on the product, the reactant in the flask is
  free 3, and the difference between the two is the equilibrium isotope
  effect of binding.
