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
`mace_mp[:model]`, `mace_off[:model]`, `mace_omol`, `orb`, `sevennet`,
`aimnet2`, or `module.path:callable` for anything else; the callable is
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

The DFT-quality Claisen structures needed to reproduce the numbers above
with an actual potential are not part of the repository (see the
implementation plan, Phase 9, for the validation suite that will hold them).
