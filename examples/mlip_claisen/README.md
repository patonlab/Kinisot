# ASE input: Hessians from any calculator, including machine-learned potentials

Kinisot reads Hessians saved in ASE's `VibrationsData` JSON form and can
compute them itself with any ASE calculator (`--calc`). The fixture files
in `tests/data/ase/` are the Claisen B3LYP/6-31G(d) Hessians written in
that form, so the numbers reproduce the [claisen](../claisen/README.md)
example and check the reader, the unit conversions (eV/Å² → Hartree/Bohr²)
and the projection path end to end.

```
cd tests/data/ase
kinisot --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso 4 -t 393 -s 0.961
kinisot --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso 4 -t 393 -s 0.961 --project
```

**Result lines** (from `expected_output.dat`):

```
                       V-ratio        ZPE        EXC       TRPF        KIE    1D-tunn   corr-KIE
C4     KIE @ 393.0 K  1.012716   1.036594   1.001970   0.978989   1.029742   1.003157   1.032993
C4 (P) KIE @ 393.0 K  1.012716   1.036593   1.001955   0.979004   1.029742   1.003157   1.032993
```

## A worked example with GFN2-xTB

`scripts/make_claisen_structures.py --calc xtb --out tests/data/xtb` starts
from the B3LYP geometries, minimizes the reactant and refines the transition
structure (with the Sella saddle point optimizer) using GFN2-xTB through its
ASE calculator (`tblite`), checks that the result is the Claisen transition
structure (see [below](#machine-learned-potentials-check-the-transition-structure-first)), and
computes the Hessians with the calculator `--calc xtb` builds (tblite's SCF
tightened to `accuracy=0.01`, see below) by central differences with a
0.005 Å step and four displacements per coordinate. The results are
committed in `tests/data/xtb/`, so the KIEs can be reproduced without
`tblite`:

```
cd tests/data/xtb
kinisot --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso 4 -t 393
cd ../../..
```

To recompute the Hessians with tblite, copy the geometries elsewhere first
(from the repository root).
Kinisot caches each Hessian next to its geometry as `<name>.hessian.json`,
which here would overwrite the committed files:

```
mkdir xtb_rerun && cp tests/data/xtb/*.xyz xtb_rerun && cd xtb_rerun
kinisot --rct claisen_gs.xyz --ts claisen_ts.xyz --iso 4 -t 393 --calc xtb --delta 0.005
```

The command line uses two displacements per coordinate, so the recomputed
Hessian differs slightly from the committed one. The corrected C4 KIE is
the same to the printed six decimals (1.020858). With the default 0.01 Å
step it is 1.020857.

The reactant is a minimum (lowest mode 58 cm⁻¹) and the transition
structure has one imaginary mode (494i cm⁻¹; B3LYP 483i unscaled). The
GFN2-xTB barrier is 20.6 kcal/mol against 28.9 for B3LYP. Corrected KIEs at
393 K (Bell tunnelling; no scaling factor exists for GFN2-xTB, so the B3LYP
values are shown both unscaled and with the 0.961 used in the
[claisen](../claisen/README.md) example):

| Position | B3LYP, s = 0.961 | B3LYP, unscaled | GFN2-xTB |
| --- | --- | --- | --- |
| C1 | 1.0147 | 1.0151 | 1.0259 |
| C2 | 1.0020 | 1.0020 | 1.0060 |
| O3 (¹⁸O) | 1.0396 | 1.0415 | 1.0273 |
| O3 (¹⁷O) | 1.0206 | 1.0216 | 1.0143 |
| C4 | 1.0330 | 1.0346 | 1.0209 |
| C5 | 1.0019 | 1.0020 | 1.0045 |
| C6 | 1.0168 | 1.0173 | 1.0282 |
| H7,H8 (²H₂) | 0.9566 | 0.9525 | 0.9010 |

GFN2-xTB describes a different transition structure: larger KIEs at the
bond-forming carbons C1 and C6, smaller ones at the bond-breaking C4–O3
pair, and a much more inverse secondary ²H₂ KIE at the hydrogens of C1 (the
vinyl CH₂ that forms the new bond). The geometries agree: the forming C1–C6
bond is 1.94 Å in the GFN2-xTB transition structure against 2.31 Å for
B3LYP, and the breaking C4–O3 bond 1.59 Å against 1.90 Å, a tighter
transition structure further along the C1–C6 bond formation and less far
along the C–O cleavage, consistent with its lower barrier. Which description matches experiment is the
question the benchmark suite (`benchmarks/claisen` and
`benchmarks/claisen_xtb`) will answer once the measured KIEs are entered.
Projection of the external modes changes these KIEs by less than
3 × 10⁻⁷ because the geometries are converged to 10⁻⁴ eV/Å; with looser
geometries it matters more.

**Converge the energy tightly for finite-difference Hessians.** With
tblite's default SCF accuracy the noise in the forces shifts isotope effects
by up to a few 10⁻⁴: for water, swapping ²H between the two equivalent
hydrogens (an EQE that must be exactly 1) gives 1.0001 or 0.9996 depending
on the stencil, against 1 ± 3 × 10⁻⁷ with `accuracy=0.01`, which is why
`--calc xtb` uses it. The same caution applies to any calculator: converge
the electronic structure well beyond what an optimization needs, and check
a symmetric EQE like this one when in doubt.

## With a machine-learned potential

Give geometries (anything ASE reads: `.xyz`, `.extxyz`, ...) of the
optimized reactant and transition structure and a calculator. The Hessian is
computed by central finite differences (or analytically when the calculator
offers `get_hessian`) and cached next to each geometry as
`<name>.hessian.json`, so a scan over positions and temperatures computes
it once. Projection of the external modes is on by default for these
inputs, because finite-difference Hessians leave translational and
rotational residuals of tens of cm⁻¹ that the lowest-six rule cannot
separate from real modes. No Truhlar scaling factor exists for a potential,
so the factor is 1.0 unless `-s` is given.

```
pip install "kinisot[ase]" mace-torch          # or orb-models, sevenn, aimnet2calc, ...
kinisot --rct claisen_gs.xyz --ts claisen_ts.xyz --iso 4 -t 393 --calc mace_mp:medium
kinisot --rct claisen_gs.xyz --ts claisen_ts.xyz --iso 4 -t 393 --calc mace_off:medium --delta 0.005
kinisot --rct rct.xyz --ts ts.xyz --iso 4 --calc my_package.calculators:make_calculator
```

`--calc` accepts `emt` (ASE's built-in test potential, not for chemistry),
`xtb[:method]` (GFN2-xTB by default, `xtb:GFN1-xTB` for GFN1; needs
`tblite`), `mace_mp[:model]`, `mace_off[:model]`, `mace_omol`, `orb[:model]`,
`sevennet`, `aimnet2`, or `module.path:callable` for anything else; the callable is
called without arguments (plus `model=` when given) and must return an ASE
calculator. The geometries must be stationary points **of the same
potential**: optimize the reactant and locate the transition structure with
the calculator you then pass to Kinisot, otherwise the Hessian has gradient
contamination and spurious imaginary modes. Kinisot warns when a
transition structure has more than one imaginary mode beyond the cutoff and
treats small imaginary modes below it as real vibrations of the same
magnitude (with a warning).

From Python:

```python
from ase.io import read
from mace.calculators import mace_mp
from kinisot import compute_kie, hessian_from_calculator

calc = mace_mp(model="medium", default_dtype="float64")
gs = hessian_from_calculator(read("claisen_gs.xyz"), calc)   # geometries optimized with the same calc
ts = hessian_from_calculator(read("claisen_ts.xyz"), calc)
r = compute_kie(rct=gs, ts=ts, iso="4", temperature=393.0)   # projection on by default for ASE Hessians
print(r.kie_tunnel, r.other.light.imaginary)
```

## Machine-learned potentials: check the transition structure first

`scripts/make_claisen_structures.py` runs the GFN2-xTB workflow above with
any `--calc`. It was run with four MACE foundation models (mace-torch
0.3.16, float64, analytic Hessians). With none of them does the saddle point
search, started from the B3LYP transition structure, find the concerted
Claisen transition structure. The script rejects every one of these saddle
points, so the repository has no MACE structures or KIEs for this reaction:

| Potential | ΔE at the B3LYP geometries (kcal/mol) | Saddle point reached | C1–C6 / C4–O3 (Å) | Imaginary modes (cm⁻¹) | Barrier (kcal/mol) |
| --- | --- | --- | --- | --- | --- |
| B3LYP/6-31G(d) | 28.9 | concerted [3,3] shift | 2.31 / 1.90 | 483i | 28.9 |
| GFN2-xTB | 28.4 | concerted [3,3] shift | 1.94 / 1.59 | 494i | 20.6 |
| MACE-OFF23 small | 60.4 | C4–O3 cleavage, no C1–C6 bond | 3.81 / 2.38 | 229i | 50.4 |
| MACE-OFF23 medium | 70.7 | C4–O3 cleavage, no C1–C6 bond | 3.56 / 2.16 | 158i | 59.3 |
| MACE-OFF23 large | 65.6 | none (not converged, two imaginary modes) | 2.03 / 1.47 | 337i, 117i | – |
| MACE-MP-0 medium | 25.1 | C1–C6 ring closure, C4–O3 intact | 2.42 / 1.46 | 355i | 9.2 |

"ΔE at the B3LYP geometries" is the energy of the B3LYP transition structure
above the B3LYP reactant, both evaluated with the potential without
re-optimizing. MACE-OFF23 was trained on near-equilibrium organic molecules
(SPICE), and it places the pericyclic region 30 to 40 kcal/mol too high,
so the saddle search slides into C–O dissociation instead. MACE-MP-0 was
trained on inorganic materials (the Materials Project). It happens to give
a reasonable ΔE at the B3LYP geometries, but its surface has no concerted
saddle point nearby. Sella reaches a C1–C6 ring closure with the C4–O3 bond
intact, 9 kcal/mol above a reactant that is not a minimum either (a 42i cm⁻¹
mode).

These saddle points are real stationary points of the potentials, and KIEs
computed from them would look plausible. They would still describe a
different reaction. So before trusting isotope effects from a potential:

- **Validate the transition structure.** The script requires exactly one
  imaginary mode, both partial bonds inside generous windows (C1–C6
  1.75–2.90 Å, C4–O3 1.50–2.50 Å), and the two stretches of those bonds
  dominating the imaginary mode. It exits with status 1 and writes nothing
  otherwise (`--keep-invalid` writes the structures for inspection). For
  your own reaction, look at the imaginary mode and connect the transition
  structure to its reactant and product (IRC or NEB).
- **Try single points first.** The potential's energies at the DFT
  reactant and transition structure are cheap to compute. A ΔE far from the
  DFT barrier, like MACE-OFF23's, means the potential does not describe that
  region. A close ΔE, like MACE-MP-0's, is necessary but not sufficient.
- **Prefer potentials trained on reactive data.** A general-purpose
  foundation model cannot be assumed to know transition-structure regions.
  Fine-tuning on reaction-path data for the reaction class, or GFN2-xTB or
  DFT, is the safer route. `mace_omol`, `orb`, `sevennet`, `aimnet2` and UMA
  were not tested here.

To reproduce the rejections:

```
pip install "kinisot[ase]" sella mace-torch
python scripts/make_claisen_structures.py --calc mace_mp:medium --out mace_mp0 --label MACE-MP-0-medium
python scripts/make_claisen_structures.py --calc mace_off:medium --out mace_off23_medium
```

`--keep-invalid` writes the rejected structures anyway. The MACE-MP-0 ones
are kept in `tests/data/mace_mp0_rejected/` as the regression test for
these checks.

MACE-MP-0 is MIT licensed. The MACE-OFF23 weights are distributed under the
Academic Software License, which does not allow commercial use. Kinisot
does not ship or download any weights itself. `mace-torch` fetches them on
first use.
