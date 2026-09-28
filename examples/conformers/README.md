# Conformer ensembles: the Claisen rearrangement of allyl vinyl ether

Allyl vinyl ether has eight conformers within 2.2 kcal/mol, and its
rearrangement can pass through a chair or a boat transition structure.
Which of them should a KIE be computed from, and does it matter? This
example computes the KIEs of the lowest pair and of the whole ensemble at
393 K, relative to C5 as in the experiment of Meyer, DelMonte and Singleton
(J. Am. Chem. Soc. 1999, 121, 10865).

**In short:**

- **¹³C and ¹⁷O KIEs hardly depend on the conformer.** Every reactant
  conformer gives the same value to within 0.002, and the ensemble is within
  2 × 10⁻⁴ of the lowest pair.
- **The α-deuterium KIE at C4, the carbon whose bond to oxygen breaks,
  does.** It runs from 0.926 to 0.984 over the eight conformers, and the
  ensemble (0.968) differs from the lowest pair (0.978) by 0.010, more than
  a typical measurement error.
- **The boat does not contribute.** It lies 6.2 kcal/mol above the chair
  in free energy and carries 0.04 % of the rate (N_eff = 1.00).

The equations are in [docs/theory.md, section 7](../../docs/theory.md#7-conformer-ensembles).

## The structures

`python scripts/make_conformer_example.py --out examples/conformers` made
them with GFN2-xTB (`tblite`, through ASE) in about two minutes, and they
are committed here as Hessians in ASE's JSON form, so the example needs
neither `tblite` nor `sella` to run.

- **Reactant conformers.** Every combination of the three rotatable
  dihedrals (C1=C2–O3–C4 at 0° and 180°, C2–O3–C4–C5 at 60°, 180° and
  300°, O3–C4–C5=C6 at 0°, 120° and 240°) was minimized from the GFN2-xTB
  reactant of [mlip_claisen](../mlip_claisen/README.md). Of the 18 starts,
  10 distinct conformers remain once mirror images and renumbered
  hydrogens are recognized, and 8 lie within 3 kcal/mol.
- **Transition structures.** The chair is the GFN2-xTB saddle point of
  [mlip_claisen](../mlip_claisen/README.md). The boat started from it with
  the allyl fragment reflected through the plane of C1, O3, C4 and C6, and
  was refined with Sella, a saddle-point optimizer. Displaced along its imaginary mode either way
  and minimized, each saddle point gives allyl vinyl ether on one side and
  4-pentenal on the other. At this level the boat is asynchronous: C4–O3 is
  1.51 Å (1.59 in the chair) and C1–C6 formation leads.
- **One numbering.** The script rotates dihedrals and never renumbers atoms,
  so a label means the same atom in every file. Kinisot checks this: the
  same elements in the same order (an error otherwise), and the same bonds
  (a warning otherwise).

| Conformer | C1=C2–O3–C4 | C2–O3–C4–C5 | O3–C4–C5=C6 | ΔE | ΔG (qRRHO, 393 K) | g | Population |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gs_1 | 358° | 70° | 212° | 0.00 | 0.00 | 2 | 28.8 % |
| gs_2 | 357° | 75° | 8° | 0.11 | 0.23 | 2 | 21.3 % |
| gs_3 | 1° | 78° | 139° | 0.68 | 0.66 | 2 | 12.4 % |
| gs_4 | 360° | 179° | 138° | 0.77 | 0.56 | 2 | 14.1 % |
| gs_5 | 360° | 180° | 360° | 0.90 | 0.76 | 1 | 5.4 % |
| gs_6 | 187° | 73° | 3° | 1.63 | 1.15 | 2 | 6.6 % |
| gs_7 | 184° | 66° | 226° | 1.68 | 1.04 | 2 | 7.6 % |
| gs_8 | 192° | 75° | 134° | 2.14 | 1.57 | 2 | 3.8 % |
| ts_chair | | | | 21.46 | 0.00 | 2 | 99.96 % of the rate |
| ts_boat | | | | 27.26 | 6.16 | 2 | 0.04 % of the rate |

Energies in kcal/mol: ΔE is the electronic energy above gs_1 (for the
transition structures, the barrier from gs_1); ΔG is the free energy
Kinisot computes for the weights, relative within each species. qRRHO is
the quasi-harmonic free energy, Kinisot's default for weights, which treats
the lowest vibrations partly as free rotations (Grimme's method, as in
GoodVibes). RRHO is the plain harmonic free energy.

**Degeneracies.** A conformer without a mirror plane has a mirror image of
the same energy, which the list does not repeat. It therefore counts twice
(g = 2): all of these except gs_5, which is planar and its own mirror image
(g = 1). Both transition structures are chiral, so their factors of 2
cancel. `degeneracies.txt` gives them to Kinisot in the `--energies`
format, with `-` for free energies Kinisot computes. Here they change the
C4-d₂ KIE by 0.0022 and the others by 3 × 10⁻⁴ or less: gs_5, the one
conformer they single out, gives the most inverse C4-d₂ KIE.

## Commands

```
cd examples/conformers
# the lowest pair
kinisot --rct gs_1.hessian.json --ts ts_chair.hessian.json --iso 10,11 -t 393 -s 1 --reference 5
# every conformer: several files after one flag, with the degeneracies
kinisot --rct gs_1.hessian.json gs_2.hessian.json gs_3.hessian.json gs_4.hessian.json \
               gs_5.hessian.json gs_6.hessian.json gs_7.hessian.json gs_8.hessian.json \
        --ts ts_chair.hessian.json ts_boat.hessian.json \
        --iso 10,11 -t 393 -s 1 --reference 5 --energies degeneracies.txt
# every position at once, from a job file
kinisot --job ensemble.json
```

Hydrogens 10 and 11 are on C4. `-s 1` because GFN2-xTB has no scaling
factor, and external modes are projected out, the default for Hessians
from ASE.

The ensemble run prints the conformer table and then the ensemble KIE
(from `expected_output.dat`; the KIE columns of the table are absolute,
the last line is relative to C5):

```
  Conformers at 393.0 K:
                                                     dG     g   pop %    V-ratio        KIE   corr-KIE
o reactant 1, iso @ 10,11
    gs_1.hessian                                   0.00     2   28.76              0.981190   0.981525
    gs_2.hessian                                   0.23     2   21.30              0.969300   0.969631
    ...
    gs_5.hessian                                   0.76     1    5.44              0.929329   0.929647
    ...
o transition structure 1, iso @ 10,11
    ts_chair.hessian                               0.00     2   99.96     1.0012   0.971384   0.971716
    ts_boat.hessian                                6.16     2    0.04     1.0006   0.958826   0.958934

                                                          KIE    1D-tunn   corr-KIE     lowest        range (+/-0.5)  N_eff
  KIE (ensemble) @ 393.0 K                          0.971380   1.000342   0.971711   0.981530   0.969419-0.973724     1.00
  relative to iso @ 5 / 5:                          0.967498              0.967755
```

`ensemble.json` lists the same files, the degeneracies and eight
isotopologues ([docs/job_files.md](../../docs/job_files.md)).

## Results

All relative to C5, at 393 K, with Bell tunnelling.

| Position | Lowest pair | Over the 8 reactant conformers | Ensemble (qRRHO) | RRHO | Equal weights | Measured |
| --- | --- | --- | --- | --- | --- | --- |
| C1 ¹³C | 1.0220 | 1.0209–1.0226 | 1.0219 | 1.0218 | 1.0203 | 1.0135 |
| C2 ¹³C | 1.0017 | 1.0009–1.0020 | 1.0015 | 1.0014 | 0.9995 | 1.0005 |
| O3 ¹⁷O | 1.0099 | 1.0086–1.0104 | 1.0097 | 1.0095 | 1.0050 | 1.0190 |
| C4 ¹³C | 1.0165 | 1.0157–1.0170 | 1.0165 | 1.0164 | 1.0108 | 1.0340 |
| C6 ¹³C | 1.0238 | 1.0231–1.0243 | 1.0238 | 1.0238 | 1.0236 | 1.0150 |
| C1-d₂ | 0.9004 | 0.8837–0.9004 | 0.8959 | 0.8949 | 0.8772 | |
| C4-d₂ | 0.9777 | 0.9255–0.9844 | 0.9678 | 0.9685 | 0.9591 | |
| C6-d₂ | 0.8732 | 0.8732–0.8817 | 0.8760 | 0.8756 | 0.8745 | |

The measured values are the means of the two experiments in the paper's
Table 4 (as in `benchmarks/claisen`); no deuterium KIEs were measured.

**What the numbers show:**

- **Secondary deuterium KIEs are where conformers matter.** At C4 the two
  conformers with C2–O3–C4–C5 near 180° (gs_4 and gs_5) give 0.946 and
  0.926, the gauche ones 0.958–0.984. Heavy-atom KIEs change by 0.002 at
  most.
- **With one transition structure, the ensemble is the population-weighted
  mean of the pairwise KIEs.** With the chair alone, that mean and the
  ensemble agree to 10⁻¹⁵. When several transition structures share the
  rate, the transition-structure side becomes a rate-weighted harmonic mean
  instead ([theory.md, section 7](../../docs/theory.md#7-conformer-ensembles)).
- **qRRHO and RRHO weights agree** to 0.001 here. Equal weights count the
  boat as much as the chair and move even the ¹³C KIEs by up to 0.006,
  which shows that the weights on the transition-structure side matter.
- **The weights are uncertain too.** Moving each free energy by ±0.5
  kcal/mol in turn moves the C4-d₂ ensemble KIE between 0.9694 and 0.9737
  (absolute), a spread of 0.004: comparable with the conformer effect
  itself, so the free energies deserve as much care as the Hessians.
- **An ensemble does not fix the level of theory.** GFN2-xTB puts the C4
  ¹³C KIE at 1.0165 against the measured 1.034, and C6 at 1.024 against
  1.015. Its chair is tighter than B3LYP's: C4–O3 1.59 Å and C1–C6 1.94 Å,
  against 1.90 and 2.31 Å, so C–O cleavage is less advanced and C–C
  formation more. B3LYP/6-31G(d) reproduces the measurements
  ([claisen](../claisen/README.md), `benchmarks/REPORT.md`). This example
  shows the workflow; for real work, use DFT geometries and a conformer
  search such as CREST.

## From Python

```python
from kinisot import Conformers, compute_kie

reactant = Conformers(["gs_%d.hessian.json" % k for k in range(1, 9)], degeneracy=[2, 2, 2, 2, 1, 2, 2, 2])
r = compute_kie(rct=[reactant], ts=[["ts_chair.hessian.json", "ts_boat.hessian.json"]],
                iso="10,11", temperature=393, scale=1.0, reference="5")
r.kie_tunnel_relative                  # 0.9678
r.kie_tunnel_lowest, r.n_effective     # 0.9815 (absolute), 1.0007
[(c.name, c.free_energy, c.population, c.kie_tunnel) for c in r.conformers]
```
