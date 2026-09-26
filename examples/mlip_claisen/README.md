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

`scripts/make_xtb_claisen.py` starts from the B3LYP geometries, minimizes
the reactant and refines the transition structure (with the Sella saddle
point optimizer) using GFN2-xTB through its ASE calculator (`tblite`), and
computes the Hessians exactly as `--calc xtb` would (central differences,
0.005 Å, with tblite's SCF tightened to `accuracy=0.01`, see below). The
results are committed in `tests/data/xtb/`, so the KIEs can be reproduced
without `tblite`:

```
cd tests/data/xtb
kinisot --rct claisen_gs.hessian.json --ts claisen_ts.hessian.json --iso 4 -t 393
kinisot --rct claisen_gs.xyz --ts claisen_ts.xyz --iso 4 -t 393 --calc xtb      # recomputes with tblite
```

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
`tblite`), `mace_mp[:model]`, `mace_off[:model]`, `mace_omol`, `orb`,
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

Machine-learned potential weights could not be downloaded where this
example was prepared, so the potential-optimized structures in the
repository are the GFN2-xTB ones above; the same script works with any
calculator (replace `calculator()` in `scripts/make_xtb_claisen.py`).
