# Kinisot Review: Usability, Correctness, and ORCA/ASE Support

**Date:** 2026-09-25
**Version reviewed:** 2.0.3 (branch head 826faa1, all Phase 0–1 fixes applied)
**Companions:** [AUDIT.md](AUDIT.md) (2026-07-02 code audit) and
[IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) (revised alongside this review).

## 1. Summary

1. **The numbers are right.** Kinisot reproduces Gaussian's printed harmonic
   frequencies to within 0.03 cm⁻¹ on all eight bundled outputs, the six
   modes it discards coincide with Gaussian's "Low frequencies" line, and its
   Bigeleisen–Mayer and Bell implementations are algebraically identical to
   PyQuiver's (Section 3). Where the two codes differ it is PyQuiver that
   still carries the atomic-mass-unit typo Kinisot removed in 2.0.3.
2. **The plan's central assumption is stale, in Kinisot's favour.** Phase 5
   assumed a PR to GoodVibes would be needed to expose the Cartesian
   Hessian. GoodVibes 4.4.0 already ships `HessianData`, `parse_hessian()`
   (Gaussian archive and ORCA `.hess`), `QCData.per_atom_masses`, an
   ORCA-aware `level_of_theory()`, and a Truhlar v5 scaling table with a
   cross-program canonicalizer. An ORCA path through that code was
   prototyped here and reproduces the Claisen golden KIE to 1×10⁻⁹
   (Section 2.3). ORCA support is an adapter, not an upstream project.
3. **The largest gap is user-facing, not numerical.** The output table
   cannot be interpreted without reading the source; the `'0'` label, the
   multi-reactant convention, and the EQE mode are undocumented; only three
   isotopes exist (and the oxygen default is ¹⁷O, not the commonly measured
   ¹⁸O); `--iso` silently ignores atoms it cannot substitute; there is no
   Python API, no machine-readable output, and no console script.
4. **"Sort eigenvalues and drop the lowest 5/6" does not survive contact
   with MLIP Hessians.** It is fine for tightly converged Gaussian and ORCA
   jobs (Section 2.1) but finite-difference Hessians from machine-learned
   potentials routinely leave translation/rotation residuals of 10–30 cm⁻¹,
   above genuine low-frequency modes, and transition structures can carry
   extra small imaginary modes. Eckart projection has to ship with, and
   default on for, the ASE backend (Section 2.4, Section 6).

## 2. What was verified

Environment: Python 3.11, NumPy 2.4, GoodVibes 4.4.0, ASE 3.29, pytest 9.1.
The test suite passes (23/23). The checks below were run with throwaway
scripts; the plan turns them into permanent tests (Phase 0, item 4).

### 2.1 Frequencies versus Gaussian (all bundled outputs)

Kinisot's kept vibrational modes were compared with the `Frequencies --`
lines printed by Gaussian, and the six discarded modes with Gaussian's
`Low frequencies ---` line (which Gaussian prints *before* projection).

| File | Modes (Gaussian / Kinisot) | max \|Δν\| / cm⁻¹ | Six discarded modes (Kinisot) | Gaussian low modes |
| --- | --- | --- | --- | --- |
| claisen_gs.out | 36 / 36 | 0.009 | −0.2 −0.1 0.1 5.6 17.2 24.0 | −0.0 0.0 0.0 5.6 17.2 24.0 |
| claisen_ts.out | 36 / 36 | 0.007 | −10.1 −0.2 0.1 0.1 13.9 17.7 | −10.1 0.0 0.0 0.0 13.9 (after −482.7) |
| claisen_prd.out | 36 / 36 | 0.007 | −6.4 −2.7 −0.1 0.2 0.2 2.7 | −6.4 −2.7 0.0 0.0 0.0 2.7 |
| DATS.out | 60 / 60 | 0.007 | −0.1 0.1 0.2 1.6 4.9 6.6 | 0.0 0.0 0.0 1.6 4.9 (after −444.6) |
| DATS_rct.out | 60 / 60 | 0.022 | −6.7 −0.2 −0.1 0.1 2.1 4.0 | −6.7 0.0 0.0 0.0 2.1 4.0 |
| diene.out | 33 / 33 | 0.007 | −0.2 0.0 0.1 1.5 5.4 6.0 | 0.0 0.0 0.0 1.5 5.4 6.0 |
| dienophile.out | 21 / 21 | 0.007 | −0.1 0.1 0.1 3.3 5.6 5.8 | 0.0 0.0 0.0 3.3 5.6 5.8 |
| tetramethylcyclohexane.out | 84 / 84 | 0.007 | −4.9 −0.1 0.2 0.2 6.2 11.6 | −4.9 0.0 0.0 0.0 6.2 11.6 |

The residual 0.007 cm⁻¹ is the difference between Gaussian's internal
constants and Kinisot's CODATA 2010 values, plus Gaussian's projection. The
mode counts match in every case, so "drop the lowest six (seven with an
imaginary mode)" is validated for these inputs. Note how close the call
already is for `claisen_ts.out`: the sixth discarded mode is 17.7 cm⁻¹ and
the first kept mode is 170.7 cm⁻¹; a floppier molecule with a genuine 15 cm⁻¹
torsion would be misclassified.

### 2.2 Hessian parsing versus GoodVibes

`goodvibes.io.parse_hessian('claisen_gs.out')` mass-weighted with its own
`masses` differs from Kinisot's `read_hess()` matrix by at most 2.8×10⁻¹⁷.
The two parsers read the same archive block; GoodVibes additionally validates
the triangular length against `NAtoms=`, handles the `|` archive separator
written by Windows builds, and raises `ValueError` instead of calling
`sys.exit()`.

### 2.3 ORCA path, end to end

The Claisen Gaussian Hessians were written out in ORCA `.hess` layout
(`$hessian` in column blocks of five, `$atoms` with per-atom masses), read
back with `goodvibes.io._parse_orca_hess`, mass-weighted, and pushed through
Kinisot's unchanged partition-function code:

```
round-trip hessian max diff 0.0     masses equal True
ORCA-path claisen C5 @393K: KIE = 1.001895135   (golden 1.001895135)
```

So the only Kinisot-side work for ORCA is (a) routing `read_hess` through
`parse_hessian`, (b) taking the level of theory from GoodVibes' ORCA-aware
parser, and (c) determining linearity from something other than a Gaussian
`Rotational constants (GHZ):` line. Section 5 details this.

### 2.4 ASE path, end to end

- `ase.vibrations.VibrationsData.from_2d(atoms, H × Hartree/Bohr²)` built from
  the Gaussian Hessian reproduces Kinisot's frequencies to 1.6×10⁻⁶ cm⁻¹ (the
  residual is CODATA 2010 versus 2014 for the atomic mass unit).
- `VibrationsData.with_new_masses()` gives the ¹³C isotopologue to the same
  precision, so ASE already provides the isotopologue re-weighting primitive.
- A finite-difference Hessian of ethane from ASE's built-in EMT calculator
  (a stand-in for any MLIP; EMT is not a chemical potential and the imaginary
  modes are its artefact, which is the point) gives, unprojected:

  ```
  lowest 8:  -135.4 -135.4 -98.0 -98.0 -16.2 -14.0 -10.6 -0.0
  ```

  Kinisot's rule would discard −135.4, −135.4, −98.0, −98.0, −16.2, −14.0 and
  keep −10.6, 0.0, 0.1, 0.6, 3.9 as "vibrations". After Eckart projection of
  the six external modes the same Hessian gives

  ```
  lowest 8:  -135.4 -135.4 -98.0 -98.0 -16.2  -0.0  -0.0  -0.0
  ```

  and the six zeros can be removed by value, not by position.

### 2.5 Scaling factors versus GoodVibes (Truhlar v3b2 versus v5)

- All bundled examples resolve to the same factor in both tables
  (B3LYP/6-31G(d) → 0.977).
- One value changes: M06-2X/6-31+G(d,p) ZPE factor 0.967 → 0.968.
- Seven of Kinisot's 177 rows do not resolve through GoodVibes'
  `canonicalize_level()`:
  - `MN15-L/MG3S`, `MN15-L/maug-cc-pVTZ`, `MN12-L/MG3S`,
    `MN12-SX/6-311++G(d,p)`: GoodVibes keys are `MN15L/…`, `MN12L/…`,
    `MN12SX/…`; three `FUNCTIONAL_ALIASES` entries fix this (small GoodVibes
    PR).
  - `M06-L(DKH2)/aug-cc-pwcVTZ-DK`: hyphen handling in the basis; also an
    alias fix.
  - `M06-2X/maug-cc-pVTZ`: v5 lists `maug-cc-pV(T+D)Z` instead.
  - `PW6B95/6-31+G(d,p)`: absent from v5.

## 3. Correctness versus PyQuiver

PyQuiver (`pip install pyquiver-kie`, Apache-2.0, Kwan group) is the closest
peer. Its `pyquiver/kie.py`, `quiver.py`, `constants.py`, `tunneling.py`
and `weights.dat` were read for this comparison.

| Aspect | Kinisot 2.0.3 | PyQuiver (master) | Verdict |
| --- | --- | --- | --- |
| Hessian source | Gaussian archive block | Gaussian archive block; ORCA `.hess` | PyQuiver ahead |
| Mass weighting | H/√(mᵢmⱼ) | identical | same |
| Translations/rotations | not projected; drop lowest 5/6 | not projected; `DROP_NUM_LINEAR = 5` (+1 non-linear) | same approach, same weakness |
| Linearity | rotational constants | bond-vector parallelism, `LINEARITY_THRESHOLD = 1e-6` | both program-specific; see Section 5 |
| Imaginary mode | lowest mode, if \|ν\| > 50 cm⁻¹ | `f < -imag_threshold` (default 50) | same |
| Scaling | one factor, all modes incl. imaginary, after diagonalization | identical | same |
| Product / ZPE / excitation terms | log-sum form, `exp` at the end | direct products | algebraically identical |
| Tunnelling | Bell infinite parabola | Wigner, Bell, Skodje–Truhlar (needs barrier) | PyQuiver ahead |
| Bell formula | `(ν_L/ν_H) · sin(u_H/2)/sin(u_L/2)`, applied on top of the ν ratio | `raw_ratio · (u_H/u_D) · sin(u_D/2)/sin(u_H/2)` | identical |
| Constants | CODATA 2010, `amu = 1.660538921e-27` | `amu = 1.660468E-27` (typo), 4-s.f. `h`, `c`, `kB` | Kinisot correct; PyQuiver KIEs shift by < 1×10⁻⁵ |
| Isotopes | ²H, ¹³C, ¹⁷O only, chosen by mass match | `weights.dat`: ¹H ²D ³T, ¹²C ¹³C ¹⁴C, ¹⁴N ¹⁵N, ¹⁶O ¹⁷O ¹⁸O, ¹⁸F ¹⁹F | PyQuiver ahead |
| Specifying substitutions | atom index only, on the CLI | `(atom, atom, "13C")` tuples, config file / dict | PyQuiver clearer |
| Reference isotopologue | none | `apply_reference()` divides by a reference KIE | PyQuiver ahead (needed for Singleton-style relative KIEs) |
| Multiple reactant files | yes (`--rct a --rct b`, `'0'` label) | one ground-state file | Kinisot ahead |
| EQE / EIE | explicit `--prd` | auto when no imaginary mode | same capability |
| Scaling-factor lookup | automatic (Truhlar) | manual | Kinisot ahead |
| Output | fixed-width text + `.dat` | table, dict, CSV, `to_dataframe()` | PyQuiver ahead |
| Python API | `compute_isotope_effect()` returning an 8-tuple | `KIE_Calculation`, `Config.from_dict()`, batch API | PyQuiver ahead |

Conclusion: on a Gaussian input the two programs implement the same
algorithm and should agree to the constant-induced 10⁻⁵ level. A permanent
cross-validation test against PyQuiver on the bundled examples (tolerance
2×10⁻⁵) is cheap insurance and is added to the plan (Phase 0, item 5).

### 3.1 Issues found in this review (not in AUDIT.md)

1. **Silent no-op substitutions.** `--iso 99` on the 14-atom Claisen
   reactant changes nothing and prints nothing; so does `--iso` on any atom
   that is not ¹H, ¹²C or ¹⁶O to five decimals (N, F, S, Cl, ²H already
   substituted in Gaussian, …). Both must be errors.
2. **Oxygen default is ¹⁷O.** The label used experimentally is ¹⁸O.
   Keep ¹⁷O reachable through the explicit isotope syntax, but the plain
   numeric label should mean the standard heavy label (²H, ¹³C, ¹⁵N, ¹⁸O).
3. **Additional imaginary modes are absorbed silently.** A second imaginary
   mode with \|ν\| > cutoff becomes one of the "six discarded modes", which
   pushes a genuine external-mode residual into the vibrational product. If
   a negative frequency does survive the drop, `np.log()` yields `nan` with
   no message. Count imaginary modes and warn or fail.
4. **`--cutoff` is not what its name says.** It only decides whether the
   lowest mode is the reaction coordinate; it is not a low-frequency
   treatment of the partition function. Rename to `--imag-cutoff` (keep the
   old spelling as an alias) and say so in `--help`.
5. **Per-species rows print the TRPF column swapped.** The reactant row shows
   `exp(PF_TS,light − PF_TS,heavy)` and the TS row shows the reactant
   quantity. The layout was chosen so that row 1 / row 2 equals the final
   line in every column, but it is undocumented and mislabels the physics.
   Print each species' own RPFR components and let the final line be the
   ratio (this is also what PyQuiver prints).
6. **`Species:` lines go to stdout only.** They never reach
   `Kinisot_output.dat`, so the file does not record which isotopologue it
   describes (compare the bundled `Kinisot_output.dat`, which has no label
   lines, with `claisen_kinisot.dat`, which was captured from stdout).
7. **Same-numbering assumption is unvalidated.** With one `--iso` the label
   is duplicated for reactant and TS. If atom 5 is carbon in one file and
   hydrogen in the other, two different isotopes are swapped without
   warning. Check that substituted atoms carry the same element (or the same
   light mass) across species.
8. **Scaling-factor type is fixed to ZPE.** The Truhlar table also carries
   harmonic and fundamental factors; make the choice explicit
   (`--scale-type zpe|harm|fund`, default `zpe` to preserve behaviour) and
   document why ZPE is the default for isotope effects.
9. **Program masses are trusted blindly.** Gaussian uses pure most-abundant
   isotopes (C = 12.000), but ORCA and ASE default to standard atomic
   weights (C = 12.011, H = 1.008). Mixing conventions changes the *light*
   isotopologue and would make the same Hessian give different KIEs from
   different programs. Kinisot should build both isotopologues from its own
   isotope table (Section 6.3).
10. **Constants are CODATA 2010.** Harmless (10⁻⁸ level) but update to 2018
    and cite them when the constants move into `thermo.py`.

## 4. Usability recommendations

### 4.1 README

The current README is a paragraph of history plus a one-line usage string.
Replace it with:

1. **A five-line quick start** that installs, runs the Claisen example, and
   shows the resulting table, followed by *one paragraph explaining every
   column*: V-ratio (ν‡_light/ν‡_heavy), ZPE, EXC, TRPF (the three
   Bigeleisen–Mayer factors as ratios reactant/TS), KIE (their product
   including V-ratio, "semiclassical"), 1D-tunn (Bell correction factor),
   corr-KIE (KIE × 1D-tunn, the number to report).
2. **Input conventions**, each with a runnable example: one reactant + TS
   (KIE); two reactant files with the `'0'` label (bimolecular KIE);
   reactant + product (`--prd`, EQE); labels when atom numbering differs
   between files; multiple simultaneous substitutions (`--iso 7,8`).
3. **Isotopes supported** and the explicit isotope syntax once Phase 7 lands.
4. **What is and is not done**: harmonic, no projection of external modes
   (until Phase 7), one imaginary mode, Bell tunnelling only, temperature
   is applied at evaluation time, scaling factors come from Truhlar's
   database and which factor type is used.
5. **Program support matrix** (Gaussian, ORCA, ASE/MLIP) with the exact
   files needed for each.
6. **Python API** example returning a result object.
7. **How to cite** (Zenodo DOI, and a `CITATION.cff` so GitHub shows the
   "Cite this repository" button), **how to get help**, and links to
   `docs/theory.md` and the PyQuiver comparison.

### 4.2 Output

- Print each species' own RPFR components (fix 3.1 item 5), write the
  species/isotopologue header lines into the `.dat` as well as stdout
  (item 6), and print the imaginary frequencies of both isotopologues on
  labelled lines.
- Add `--output PATH` (default `Kinisot_output.dat`, never silently
  overwrite without `--overwrite`), `--json` and `--csv` for
  machine-readable results, and `--quiet`.
- Print the level of theory, the scaling factor *and its type and source*,
  and the number and values of the modes discarded as external modes, so a
  user can spot the 17.7 cm⁻¹ case in Section 2.1.

### 4.3 CLI

- Console entry point `kinisot` (keep `python -m kinisot`).
- `--version`; `parse_args()` so typos are errors; refuse `--ts` together
  with `--prd`; validate atom indices and element identity (3.1 items 1, 7).
- Explicit isotope syntax `--iso 5:13C,7:2H` alongside the numeric form.
- `--tunneling none|bell|wigner` (Skodje–Truhlar later, it needs energies).
- `--project/--no-project`, `--imag-cutoff`, `--scale-type`.

### 4.4 Python API

`compute_isotope_effect()` returns an eight-tuple whose first element is a
list of `calc_rpfr` instances. Replace with

```python
from kinisot import compute_kie
r = compute_kie(rct=["claisen_gs.out"], ts=["claisen_ts.out"], iso="5", T=393, scale=0.961)
r.kie, r.kie_tunnel, r.zpe, r.exc, r.trpf, r.imag_ratio, r.species[0].frequencies
```

with a frozen dataclass result, a `to_dict()`, and a documented notebook.
Keep `compute_isotope_effect` as a deprecated shim for one minor release.

### 4.5 Examples

- Move the Gaussian outputs to `tests/data/` and keep a curated top-level
  `examples/` with one directory per case (`claisen/`, `diels_alder/`,
  `eqe_cyclohexane/`, later `orca_claisen/`, `mlip_claisen/`), each with a
  `README.md` stating the chemical question, the command, the expected
  table, and the literature value where one exists.
- Replace the `.sh` files with a single `run_examples.sh` and a
  `examples.ipynb` that uses the Python API and plots KIE versus temperature.

### 4.6 Documentation

A `docs/` folder (MkDocs later if wanted) with `theory.md` (the exact
equations as implemented, sign conventions, what "TRPF" means, the Bell
formula, the scaling convention), `file_formats.md` (what each program must
produce: Gaussian `freq` with archive; ORCA `.out` + `.hess`; ASE
`VibrationsData` JSON), `faq.md` (why my KIE differs from Gaussian's
`iso=` route, why the ¹⁷O default changed, how to treat multiple
conformers), and `comparison.md` (Section 3 of this review, kept current).

### 4.7 Project hygiene

`CITATION.cff`, `CONTRIBUTING.md` (how to run tests, how to add a backend),
issue templates, `pyproject.toml`, a `kinisot[ase]` extra, ruff in CI, and
tag-driven publishing. All already in the plan; the ordering is moved
earlier because the documentation items block adoption more than the
refactor does.

## 5. ORCA support: design

**Inputs.** The ORCA output (`name.out`) and the Hessian file ORCA writes
next to it (`name.hess`). Users pass either file; the other is located by
stub, exactly as `goodvibes.io.parse_hessian()` already does.

**What GoodVibes 4.4 provides (verified in this environment).**

| Need | GoodVibes 4.4 API | Notes |
| --- | --- | --- |
| Cartesian Hessian, Eh/Bohr² | `parse_hessian(file) -> HessianData` | Gaussian archive or ORCA `$hessian`; symmetrized; dimension checked |
| Per-atom masses | `HessianData.masses` | ORCA `$atoms` respects mass overrides in the input |
| Level of theory | `level_of_theory(file)` | ORCA-aware (`_parse_orca_lot`), strips R/U |
| Scaling factor | `vib_scale_factors.scaling_data_dict[canonicalize_level(lot)]` | Truhlar v5; resolves ORCA `PBE0` versus Gaussian `PBE1PBE`, `M06-2X` versus `M062X` |
| Linear molecule | `parse_qcdata(file).linear_mol` / `.roconst` | program-specific parsing, see below |
| Program detection | `_detect_program(lines)` | Gaussian, ORCA, NWChem, xtb, Q-Chem, ASE extxyz |

**Kinisot changes.**

1. `read_hess(file, iso)` becomes: `hd = parse_hessian(file)`; build the
   mass vector from Kinisot's isotope table (Section 6.3) using
   `hd.masses` to identify the element of each atom; apply substitutions;
   mass-weight. The Gaussian path is numerically unchanged (Section 2.2).
2. `level_of_theory()` and `find_scaling_factor()` delegate to GoodVibes;
   Kinisot's `vib_scale_factors.py` is deleted after the Section 2.5 alias
   PR lands upstream (or with a temporary local alias shim).
3. `is_linear()` stops reading Gaussian text. Determine linearity from the
   geometry (principal moments of inertia, one ≈ 0) or, once projection
   exists, from the rank of the external-mode projector (5 versus 6). This
   is program-independent and also serves ASE.
4. Frequency self-check: ORCA prints `VIBRATIONAL FREQUENCIES` (projected);
   after Phase 7 projection Kinisot can compare its unsubstituted frequencies
   against the program's and warn if they differ by more than ~1 cm⁻¹,
   which catches unit or mass-convention mistakes for any backend.
5. Tests: two small ORCA fixture pairs (a reactant/TS pair and a linear
   molecule), each `.out` trimmed to the sections the parsers read so the
   files stay under 100 kB; golden values produced by converting the ORCA
   `.hess` to Gaussian-archive form and running the existing path (the
   inverse of the Section 2.3 prototype).

**ORCA-specific caveats to document.** `NumFreq` Hessians are noisier than
analytic ones (projection recommended); `%freq scalfreq` scales ORCA's
printed frequencies but not `$hessian`, so Kinisot's factor is applied to
raw frequencies as intended; ORCA 6 keeps the `.hess` layout but adds
sections Kinisot ignores; ORCA's default masses are standard atomic weights
(3.1 item 9).

## 6. ASE / MLIP support: design

### 6.1 Two entry points

1. **Consume an existing ASE result.** `kinisot --rct rct_vib.json --ts
   ts_vib.json …` where each file is `VibrationsData.todict()` written as
   JSON (or the `vib.*.json` cache directory of `ase.vibrations.Vibrations`).
   Detection: JSON with `hessian` and `atoms` keys.
2. **Compute the Hessian inside Kinisot.**

   ```python
   from kinisot.ase import hessian_from_calculator
   hd_rct = hessian_from_calculator(atoms_rct, calc, delta=0.01, nfree=2)
   hd_ts  = hessian_from_calculator(atoms_ts,  calc)
   compute_kie(rct=[hd_rct], ts=[hd_ts], iso="5:13C", T=393, project=True)
   ```

   `hessian_from_calculator` runs `ase.vibrations.Vibrations` (central
   differences) or, when the calculator exposes an analytic Hessian
   (`calc.get_hessian`, as MACE and AIMNet2 do), uses that. Every current
   MLIP (MACE, UMA/FAIRChem, ORB, SevenNet, MatterSim, AIMNet2, TorchANI)
   is an ASE calculator, so this single adapter covers all of them, and a
   CLI form `kinisot --rct rct.xyz --ts ts.xyz --calc mace_mp:medium …`
   can follow once the API is stable.

### 6.2 Interchange type

Kinisot's internal container (a frozen dataclass; GoodVibes' `HessianData`
plus the fields Kinisot needs):

```
hessian   (3N,3N) Eh/Bohr²        symbols  [N]         positions (N,3) Bohr
masses    [N] amu (light isotopologue, from Kinisot's table)
program   'Gaussian'|'Orca'|'ase'  source  path        level_of_theory  str|None
```

Positions are required for projection; `level_of_theory` is `None` for
MLIPs (no Truhlar factor; scale 1.0 unless `-s` is given, with a printed
note). Units from ASE: multiply the eV/Å² Hessian by `Bohr²/Hartree`
(verified to 1.6×10⁻⁶ cm⁻¹ in Section 2.4).

### 6.3 Masses

Neither ASE (`ase.data.atomic_masses`: C = 12.011) nor ORCA use pure
isotopic masses. Kinisot must own an isotope table (most-abundant isotope as
the light default; explicit labels ²H ³T ¹³C ¹⁴C ¹⁵N ¹⁷O ¹⁸O ¹⁸F ³⁴S ³⁷Cl …
from NIST/AME2020) and build *both* isotopologues from it, using program
masses only to identify elements and to honour explicit user mass overrides.
This keeps Gaussian results unchanged and makes ORCA/ASE results
program-independent.

### 6.4 External modes

Projection is mandatory here (Section 2.4). Implement Eckart projection in
mass-weighted Cartesians (three translations, three or two rotations built
from the centre-of-mass geometry, Gram–Schmidt, `P = 1 − DDᵀ`, diagonalize
`P H P`), then drop the modes with \|ν\| < 1 cm⁻¹ by value. Identify the
reaction coordinate as the largest-magnitude imaginary mode *after*
projection and warn about any others. Make projection available for all
backends (`--project`), default on for ASE, and report the KIE difference
with and without projection on the bundled Gaussian examples in the
changelog so users know what to expect (expected 10⁻⁵–10⁻⁴ relative).

### 6.5 Packaging and tests

`pip install kinisot[ase]` (ASE pulls SciPy and Matplotlib, so keep it
optional). Tests use ASE's built-in EMT calculator, which needs no model
download; an optional MACE test is skipped unless `mace-torch` is installed.
Ship `examples/mlip_claisen/` reproducing the DFT Claisen KIEs with
MACE-MP-0 or UMA and stating the discrepancy honestly.

## 7. Changes made to the implementation plan

- Phases 0–1 marked done; new Phase 0 items for the cross-checks in
  Section 2 (Gaussian-printed frequencies, GoodVibes parity, PyQuiver
  parity) so they run in CI.
- Documentation and UX moved up to Phase 3 (was the tail of Phase 5).
- Phase 5 (GoodVibes) rewritten: no upstream Hessian PR is needed; the only
  upstream ask is the scaling-factor aliases in Section 2.5. ORCA support is
  delivered in this phase.
- New Phase 7 items: isotope table with explicit syntax and validation
  (3.1 items 1, 2, 7, 9), Eckart projection, imaginary-mode accounting,
  tunnelling options, reference isotopologue, scale-type selection.
- New Phase 8: ASE/MLIP backend (Section 6).
- Version targets and the sequencing table updated accordingly.
