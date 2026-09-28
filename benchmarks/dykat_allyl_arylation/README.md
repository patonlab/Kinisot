# Rh-catalysed arylation of racemic 3-chlorocyclohexene: ¹³C KIEs

Van Dijk et al. (Nat. Catal. 2021, 4, 284) measured ¹³C KIEs on the product,
3-phenylcyclohexene, at 40 °C. Both enantiomers of the allyl chloride 3 react
through their own anti-oxidative-addition transition structure and converge
on one η³-allyl intermediate, so each product carbon comes from a different
substrate carbon in each enantiomer's channel (see `case.json`). The KIEs
here need the frequency jobs below; the combined KIE also needs the channel
treatment planned in IMPLEMENTATION_PLAN.md, Phase 10.

## Jobs to run

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

Each title line carries these values. A job at the right geometry and level
reproduces its energy to about 10⁻⁶ hartree, and each transition structure
has exactly one imaginary mode.

- **Basis functions.** The four Rh inputs give 6-31G(d) through general
  basis input (`genecp`) without `5D` or `6D`, so Gaussian uses its default
  for general basis input: spherical d functions. The `6-31G(d)` keyword
  alone would use Cartesian ones. The paper's jobs most likely used the same
  default, and the title-line energy shows whether they did. If a Rh job's
  energy is off by millihartrees, rerun it with `6D` added to the route.
- **Free 3 in both conventions.** `allyl_chloride_3.gjf` uses the keyword
  (Cartesian d, presumably as in the paper, so its energy can be checked).
  `allyl_chloride_3_5d.gjf` is the same 16-atom job with spherical d
  functions. A KIE from free 3 needs the same functions as the transition
  structure, so it uses whichever job matches the convention of the Rh
  jobs.
- **A transition structure matching the other one's energy.** Table 20
  labels them "(R)-TS complex" and "(S)-TS complex"; if the (R) input gives
  the (S) energy, tell me, and the labels will be swapped here.

Send back the `.log` files, including both jobs for 3. The paper computed each channel from the
reactant complex. The free allyl chloride is there as well, because for an
intermolecular KIE measured on the product the reactant in the flask is
free 3, and the difference is the equilibrium isotope effect of binding.
