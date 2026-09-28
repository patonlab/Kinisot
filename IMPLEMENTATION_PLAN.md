# Kinisot Strategic Implementation Plan

Step-by-step plan addressing the findings in [AUDIT.md](AUDIT.md)
(2026-07-02) and [REVIEW.md](REVIEW.md) (2026-09-25). Ordered so that each
phase leaves the repo in a working, releasable state; later phases build on
earlier ones.

**Revision history**

- 2026-07-02: initial plan from the code audit.
- 2026-07-02: Phase 5 redirected to reuse GoodVibes rather than rewrite.
- 2026-09-25: Phases 0–1 marked done. Documentation/UX promoted to Phase 3.
  Phase 5 rewritten around GoodVibes 4.4.0, which already ships the Hessian
  parsers the previous plan assumed would need an upstream PR; ORCA support
  moves into that phase. New Phase 7 (isotope table, projection, tunnelling
  options) and Phase 8 (ASE / machine-learned interatomic potentials).
  Evidence for every change is in REVIEW.md §2.
- 2026-09-27: Phase 10 (conformer ensembles) planned; extended the same day to
  transition structures in series and parallel channels (Dale et al. 2021).
- 2026-09-28: Phase 10 decisions confirmed (recommended options);
  conformer ensembles, transition structures in series, parallel channels
  and the JSON job file implemented.

**Guiding decisions**

1. *Reuse GoodVibes 4.4+* (same lab) for file parsing, Hessian extraction,
   level-of-theory detection and scaling factors. GoodVibes depends only on
   NumPy, `rich` and `pymsym`, so it is an acceptable hard dependency.
2. *Kinisot owns the physics*: isotope table, mass weighting, external-mode
   projection, Bigeleisen–Mayer, tunnelling. Nothing in GoodVibes is
   isotope-aware, and it should stay that way.
3. *Golden numbers move only deliberately.* Every numerical change is a
   separate commit that updates the characterization tests and the
   CHANGELOG with the size of the shift.
4. *Documentation is a release blocker*, not a tail item: a release that
   adds a feature without a README section and an example is not done.

---

## Phase 0 — Safety net ✅ (done, extended)

Done (commits 36f850f…776812e): characterization tests on the bundled
examples, `pytest.approx`, GitHub Actions on Python 3.9–3.13 × three OSes,
Travis removed.

Still to add (they are the checks run by hand in REVIEW.md §2 and belong in
CI before the parser is replaced):

4. ✅ `tests/test_frequencies.py`: for every bundled Gaussian output, Kinisot's
   kept modes must match the `Frequencies --` lines to 0.05 cm⁻¹ and the
   discarded modes must match `Low frequencies ---`. Catches any regression
   in mass weighting, unit conversion, or mode dropping, for any backend.
5. ✅ (in `tests/test_orca.py` and `tests/test_pyquiver.py`) mass-weighted
   Hessian from `goodvibes.io.parse_hessian` equals Kinisot's parser to
   1e-12; PyQuiver, given the same Hessians as `quiver.System` objects (no
   `#p` output needed), agrees to 2.3e-6 on every Claisen and Diels–Alder
   KIE (uncorrected, Bell, Wigner). The original text asked for an optional
   PyQuiver (`pyquiver-kie`) comparison on the Claisen and Diels–Alder KIEs
   to 2e-5 relative (skipped when PyQuiver is not installed; PyQuiver still
   carries the 1.660468e-27 amu typo, hence the loose tolerance).

**Exit criteria:** both tests green in CI (item 4 done; item 5 waits for
the GoodVibes dependency in Phase 5).

## Phase 1 — Confirmed bug fixes ✅ (v2.0.3, a development milestone released in 2.5.0)

All six items done (constants block, `NameError`, exact scaling-factor
match, table row, linearity from rotational constants, CHANGELOG). v2.0.3
was never published on its own; the fixes reached users in 2.5.0.

## Phase 2 — Robustness ✅ (v2.1.0, implemented 2026-09-25)

All eight items below are implemented (85 tests, 94 % coverage, every
result line of the bundled examples unchanged). Two deliberate deviations
from the text as first written:

- Item 2, mixed elements: there is **no** `--allow-mixed-elements`
  override. An isotope effect that substitutes different elements on the
  two sides is not a defined quantity, so an override would only produce
  meaningless numbers; the check is a plain error. (Different atom
  numbering between files is not affected: only the *set* of substituted
  elements must match.)
- Item 6, results file: instead of refusing to overwrite without
  `--overwrite`, new results are **appended** to the results file, which is
  how the bundled example scripts and their reference `.dat` files were
  produced (one block per substitution). `--overwrite` starts a fresh file.
  Nothing is ever lost silently, and existing scripts keep working.

1. **Exceptions, not exits**: `KinisotError` hierarchy
   (`KinisotParseError`, `KinisotInputError`) in `kinisot/exceptions.py`;
   `main()` is the only place that calls `sys.exit`.
2. **Input validation** (REVIEW §3.1 items 1, 3, 7):
   - `--iso` index out of range → error naming the file and atom count.
   - `--iso` on an atom whose mass is not in the substitution table → error
     listing the atom's element and the isotopes Kinisot knows.
   - Substituted atoms must have the same element in every species that
     shares a label position; otherwise error (override with
     `--allow-mixed-elements` for genuinely different numbering).
   - Count imaginary modes with \|ν\| > cutoff; warn on more than one in a
     TS, error on any in a reactant/product, and never let a negative
     frequency reach `log()`.
   - `level_of_theory` mismatch between files becomes a visible warning on
     stdout, not only in the `.dat`.
3. **Context managers**: `with open()` everywhere; `Logger` becomes a
   context manager and is always closed.
4. **Imports**: delete the try/except shim and `import *`; explicit relative
   imports.
5. **Numerical cleanups**: `np.log(math.exp(x))` → `x`; delete dead
   `hess_mat` and the rounded-frequency block; vectorize the three factor
   functions (`np.log`, `np.log1p(-np.exp(-u))`).
6. **CLI fixes**: `-s` default `None`; reject `--ts` with `--prd`;
   `parse_args()`; `--version`; rename `--cutoff` → `--imag-cutoff` (old
   name kept as a hidden alias); document the `'0'` label in `--help`;
   `--output PATH`, refuse to overwrite without `--overwrite`; `--quiet`.
7. **Output fixes** (REVIEW §3.1 items 5, 6): per-species rows show each
   species' own ZPE/EXC/TRPF ratios; species and isotopologue header lines
   are written to the `.dat`; discarded external modes and both imaginary
   frequencies are printed with labels.
8. Tests for every error path and a CLI smoke test via `python -m kinisot`
   on the examples asserting on the `.dat`.

**Exit criteria:** no `sys.exit`/bare `except`/unclosed file in library
code; error paths tested; coverage ≥ 80 % reported in CI.

## Phase 3 — Packaging and documentation ✅ (v2.1.0, implemented 2026-09-25)

All items below are implemented; the wheel is 23 kB, `kinisot --version`
works from a `pip install .`, ruff is clean and gates CI, and the README
quick start is the Claisen block that `tests/test_examples.py` replays.
Phase 6 item 2 (tag-driven PyPI publishing) was pulled forward into this
phase at the maintainer's request.

Packaging:

1. `pyproject.toml` (PEP 621); delete `setup.py`/`setup.cfg`; version
   single-sourced from `kinisot/__init__.py`; `requires-python >= 3.9`;
   `dependencies = ["numpy>=1.22"]` (GoodVibes is added in Phase 5).
2. Console entry point `kinisot = kinisot.cli:main`; `__main__.py` reduced to
   two lines; `__init__.py` exports `__version__` and the public API.
3. Move `kinisot/examples/` → `tests/data/`; the wheel must not ship
   megabytes of Gaussian output (target < 100 kB).
4. `ruff` (lint + format) in CI; one mechanical re-indentation commit.

Documentation (REVIEW §4; each item is a checklist entry for the release):

5. **README rewrite**: quick start with the Claisen example and its output;
   a paragraph defining every output column; the three input conventions
   (KIE, bimolecular KIE with `'0'`, EQE with `--prd`) each with a runnable
   command; supported isotopes and the default heavy label per element;
   assumptions and limitations (harmonic, external modes, one imaginary
   mode, Bell tunnelling, scaling-factor type); program support matrix;
   Python API snippet; how to cite; where to get help.
6. Top-level `examples/` with one directory per case, each with a
   `README.md` (question, command, expected table, literature value) and
   a single `run_examples.sh`; the Gaussian outputs live in `tests/data/`
   and are symlinked or referenced, not duplicated.
7. `docs/theory.md` (equations exactly as implemented, sign conventions,
   Bell formula, scaling convention), `docs/file_formats.md`,
   `docs/faq.md`, `docs/comparison.md` (REVIEW §3 kept current).
8. `CITATION.cff` (Zenodo DOI), `CONTRIBUTING.md`, issue templates.

**Exit criteria:** `pip install .` from a clean checkout via
`pyproject.toml` only; `kinisot --version` works; wheel < 100 kB; lint
clean; README quick start reproduces `claisen_kinisot.dat` verbatim.

## Phase 4 — Internal refactor and Python API ✅ (v2.2.0, implemented 2026-09-25)

Delivered as planned; `tunneling="bell"|"wigner"|"none"` (Phase 7 item 4,
except Skodje–Truhlar) came along because the API signature needed it. The
`project=` argument waits for the projection code (Phase 7).

Golden tests from Phase 0 guarantee numerical invariance.

1. Module split:
   - `kinisot/hessian.py` — `HessianInput` frozen dataclass (hessian in
     Eh/Bohr², symbols, positions in Bohr, light masses, program, source,
     level_of_theory) and mass weighting. This is the interchange type
     every backend produces (REVIEW §6.2).
   - `kinisot/backends/gaussian.py` — the current parser, returning
     `HessianInput` (replaced by the GoodVibes adapter in Phase 5, kept as
     the reference implementation for the parity test).
   - `kinisot/thermo.py` — constants (CODATA 2018, cited) and the
     Bigeleisen–Mayer / tunnelling math on frequency arrays.
   - `kinisot/cli.py` — argument parsing and all formatting.
   - `kinisot/Kinisot.py` — thin deprecation shim for one minor release.
2. Public API: `compute_kie(rct, ts=None, prd=None, iso=..., T=298.15,
   scale=None, imag_cutoff=50.0, tunneling="bell", project=False)` returning
   a frozen `IsotopeEffect` dataclass (`kie`, `kie_tunnel`, `zpe`, `exc`,
   `trpf`, `imag_ratio`, `tunnel_corr`, per-species `SpeciesResult` with
   frequencies, discarded modes, imaginary frequency, masses) with
   `to_dict()`; inputs may be paths or `HessianInput` objects.
   `compute_isotope_effect` becomes a deprecated wrapper.
3. `--json`/`--csv` output produced from `to_dict()`.
4. Replace the NumPy structured array in `vib_scale_factors.py` with a dict
   (deleted entirely in Phase 5).

**Exit criteria:** golden tests unchanged to 1e-9 relative; API documented
in README and `examples/examples.ipynb`.

## Phase 5 — GoodVibes integration and ORCA support ✅ (v2.3.0, implemented 2026-09-25)

Delivered with two deviations from the text below:

- Item 4: the in-repo Gaussian parser is **kept** (as `backends/gaussian.py`)
  rather than replaced by GoodVibes' `parse_hessian`. GoodVibes' `HessianData`
  carries the Hessian and masses only; Kinisot also needs atomic numbers, the
  archive-frame geometry (for projection), the level of theory and the
  printed frequencies, and its own parser gives more specific error messages.
  A parity test (`tests/test_orca.py::test_goodvibes_hessian_parity`) pins
  the two parsers to 1e-12; GoodVibes is used for the ORCA `$hessian` block,
  ORCA level-of-theory detection and the scaling database.
- Item 6: no real ORCA output was available, so the fixtures in
  `tests/data/orca/` are the Claisen Gaussian Hessians rewritten in ORCA's
  `.hess` layout with a minimal `.out` (documented there). They exercise the
  reader end to end against the Gaussian goldens; **real ORCA fixtures (a
  small reactant/TS pair and a linear molecule) are still wanted**, and so is
  a `#p` Gaussian pair for the PyQuiver parity test, which PyQuiver
  requires and which the bundled outputs lack (the test skips until then).


GoodVibes 4.4.0 (verified 2026-09-25) provides `goodvibes.io.parse_hessian`
→ `HessianData(hessian, masses, program, source)` for Gaussian archives and
ORCA `.hess` files, an ORCA-aware `level_of_theory()`, program detection,
and `vib_scale_factors.scaling_data_dict` keyed by `canonicalize_level()`
(Truhlar v5). The ORCA route was prototyped end to end and reproduces the
Claisen golden KIE to 1e-9 (REVIEW §2.3). No upstream Hessian work is
needed.

1. **Dependency**: `goodvibes>=4.4` in `pyproject.toml` and the conda
   recipe.
2. **Upstream PR to GoodVibes (prepared 2026-09-26: branch
   `claude/scaling-factor-aliases` on patonlab/GoodVibes, commit 20a216b,
   full GoodVibes suite passing; pull request to be opened by a maintainer)**:
   `FUNCTIONAL_ALIASES` for
   `MN15-L`, `MN12-L`, `MN12-SX`, and hyphen handling for
   `M06-L(DKH2)/aug-cc-pwcVTZ-DK`, so the seven Kinisot rows in REVIEW §2.5
   resolve. Until merged, keep a five-line alias shim in Kinisot.
3. **Scaling factors via GoodVibes**: delete `kinisot/vib_scale_factors.py`
   and `find_scaling_factor`; add `--scale-type zpe|harm|fund`
   (default `zpe`, unchanged behaviour). Document the v3b2 → v5 changes
   (M06-2X/6-31+G(d,p) 0.967 → 0.968; `M06-2X/maug-cc-pVTZ` and
   `PW6B95/6-31+G(d,p)` no longer resolve) in the CHANGELOG; update the
   scaling tests.
4. **Hessian via GoodVibes**: the Gaussian backend becomes an adapter over
   `parse_hessian()`; the Phase 0 parity test pins equality with the
   in-repo parser, which is then removed.
5. **Linearity without program text**: from principal moments of inertia
   of the parsed geometry (GoodVibes `QCData.cartesians`) or, after Phase 7,
   from the rank of the external-mode projector. Remove `is_linear()`'s
   Gaussian-specific parsing.
6. **ORCA backend**: `kinisot --rct rct.out --ts ts.out` with the `.hess`
   files found by stub, or `.hess` paths given directly. Level of theory and
   scaling factor via GoodVibes. Fixtures: one reactant/TS pair and one
   linear molecule, `.out` trimmed to the sections read (< 100 kB each);
   goldens generated by converting `.hess` → Gaussian-archive form and
   running the Gaussian path (the inverse of the REVIEW §2.3 prototype).
   Document ORCA caveats (NumFreq noise, `scalfreq`, standard atomic
   weights) in `docs/file_formats.md`, and add `examples/orca_claisen/`.
7. **Frequency self-check**: when the program's printed frequencies are
   available (`QCData.frequency_wn`), compare with Kinisot's unsubstituted,
   projected frequencies and warn above 1 cm⁻¹. Catches unit and mass
   mistakes for every backend.

**Exit criteria:** Gaussian goldens unchanged (except documented scaling
changes); ORCA goldens in CI; `docs/file_formats.md` lists exactly which
files each program must produce.

## Phase 6 — Release hygiene ✅ (ongoing; automation in place 2026-09-25)

Items 1–4 are in place: the changelog, the tag-driven PyPI workflow (which
now also creates the GitHub Release with the changelog section as notes),
Dependabot for the Actions versions, and `.zenodo.json` metadata. PyPI
trusted publishing and the Zenodo connection are configured (confirmed by
the maintainer 2026-09-25), so a pushed `v*` tag is a complete release.

1. CHANGELOG.md in Keep a Changelog format (started in 2.0.3).
2. ✅ Tag-driven publishing (`.github/workflows/publish.yml`, 2026-09-25):
   on a `v*` tag the workflow checks the tag against `__version__`, runs
   the tests, builds sdist/wheel and uploads with PyPI trusted publishing
   (environment `pypi`). One-time PyPI/GitHub setup is in CONTRIBUTING.md.
   Kinisot is not on conda-forge: `recipe/meta.yaml` is a template for a
   staged-recipes submission (CONTRIBUTING.md). Once a feedstock exists,
   the conda-forge bot picks up each PyPI release.
3. Dependabot for Actions versions.
4. Zenodo integration so each tag gets a DOI; `CITATION.cff` updated by the
   release workflow.

## Phase 7 — Isotopes, external modes, tunnelling ✅ (v2.4.0, implemented 2026-09-25)

All six items delivered. Measured on the bundled Gaussian examples,
projection changes the corrected KIEs by at most 4 × 10⁻⁷, so the default
stays off for Gaussian/ORCA through 2.x as decided, and the ASE backend
(Phase 8) will turn it on. Isotope masses come from `periodictable`
(AME 2020) through `scripts/make_isotope_data.py`; the NIST download that
ASE offers is blocked from this environment, so the generated table is
committed. Item 5 delivers a single `--reference` isotopologue rather than
a batch of isotopologues per run (the API loop and `--csv` cover batches).

1. **Isotope table** (`kinisot/isotopes.py`): most-abundant isotope masses
   and labelled isotopes (²H ³T, ¹³C ¹⁴C, ¹⁵N, ¹⁷O ¹⁸O, ¹⁸F, ³³S ³⁴S,
   ³⁷Cl, ⁸¹Br, …) from AME2020/NIST with the source cited. Both
   isotopologues are built from this table; program masses only identify
   elements (REVIEW §6.3). Gaussian numbers are unchanged because Gaussian
   already uses these masses.
2. **Explicit substitution syntax**: `--iso 5:13C,7:2H` alongside the
   numeric form; the numeric form means the *standard heavy label*
   (²H, ¹³C, ¹⁵N, ¹⁸O). This changes the oxygen default from ¹⁷O to ¹⁸O:
   a documented behaviour change with a CHANGELOG entry and a one-release
   warning whenever a bare oxygen index is used.
3. **Eckart projection** (`--project/--no-project`; default off for
   Gaussian/ORCA in this release, on for ASE in Phase 8): project
   translations and rotations from the mass-weighted Hessian, drop external
   modes by value (\|ν\| < 1 cm⁻¹), pick the reaction coordinate as the
   largest-magnitude imaginary mode after projection. Report the
   with/without-projection KIE differences on the bundled examples in the
   CHANGELOG. Revisit the default for all backends in v3.0.
4. **Tunnelling models**: ✅ `--tunneling none|bell|wigner` (Phase 4); Skodje–Truhlar
   as `--tunneling skodje --barrier <kcal/mol>` (or energies parsed from the
   files via GoodVibes `QCData.scf_energy`) once the API exists. Print all
   available corrections in the JSON output.
5. **Reference isotopologue**: `--reference <label>` divides every KIE by
   the KIE of a reference substitution (Singleton-style relative KIEs);
   multiple `--iso` sets in one run produce one table.
6. **Temperature scans**: `-t 273,298,323` or `-t 250:350:10` producing
   one row per temperature (the notebook example plots this).

## Phase 8 — ASE and machine-learned potentials ✅ (v2.5.0, implemented 2026-09-25)

Delivered as designed, with two notes. The CLI form (`--calc`) shipped in
the same release as the API rather than one release later, because the
cache next to the geometry makes it the natural way to scan positions. No
machine-learned potential could be exercised in this environment (the model
packages and weights are not installable here), so the tests use ASE's
built-in EMT potential on water, ammonia and N₂, and the MACE test is skipped
unless `mace-torch` is installed; the `examples/mlip_claisen` numbers come
from the Gaussian Hessians written in ASE's JSON form. **Update 2026-09-26:**
the Claisen reactant and transition structure were re-optimized with
GFN2-xTB through the ASE backend (`scripts/make_claisen_structures.py`, Sella
for the saddle point) and compared with B3LYP in `examples/mlip_claisen`;
this exposed that finite-difference Hessians need a tight SCF, now the
`--calc xtb` default. **Update 2026-09-26 (MACE):** the same script was run
with MACE-OFF23 (small, medium, large) and MACE-MP-0 (medium), with weights
from GitHub releases. The analytic MACE Hessian (`get_hessian`) matches
finite differences to 1.5 × 10⁻⁷. None of the four has the concerted
Claisen transition structure. MACE-OFF23 puts the pericyclic region 30 to
40 kcal/mol too high and its saddle points are C–O cleavage. MACE-MP-0's is
a C1–C6 ring closure with C–O intact. The script now validates every
structure (one imaginary mode, partial-bond windows, bonds dominating the
imaginary mode) and refuses these. No MACE KIEs are reported; the
diagnostics are in the example. **Still wanted:** a potential trained on
reactive data that passes the checks. `mace_omol`, `orb`, `sevennet`,
`aimnet2` and UMA have not been tried.

Design in REVIEW §6; prototypes verified 2026-09-25 (unit conversion to
1.6e-6 cm⁻¹, `with_new_masses` isotopologues, EMT finite-difference failure
mode without projection).

1. **Extra**: `pip install kinisot[ase]` (`ase>=3.23`); everything in
   `kinisot/backends/ase.py` imports ASE lazily.
2. **Consume ASE results**: accept `VibrationsData.todict()` JSON files and
   `ase.vibrations.Vibrations` cache directories as `--rct/--ts/--prd`
   inputs; convert eV/Å² → Eh/Bohr², positions Å → Bohr, and replace ASE's
   standard atomic weights with the Phase 7 isotope masses.
3. **Compute Hessians**: `kinisot.ase.hessian_from_calculator(atoms, calc,
   delta=0.01, nfree=2, analytic=True)` using `Vibrations` central
   differences, or the calculator's analytic Hessian when it has one
   (MACE, AIMNet2). Any ASE calculator works, which covers MACE,
   UMA/FAIRChem, ORB, SevenNet, MatterSim, TorchANI.
4. **Projection on by default** for this backend; require at most one
   imaginary mode above the cutoff after projection and print the others.
5. **Scaling**: no Truhlar factor for MLIPs; default 1.0 with a printed
   note; `-s` honoured.
6. **CLI form**: `kinisot --rct rct.xyz --ts ts.xyz --calc mace_mp:medium
   --iso 5:13C` once the Python API has been used for a release.
7. **Tests** with ASE's EMT calculator (no downloads); a MACE test skipped
   unless `mace-torch` is installed.
8. **Example**: `examples/mlip_claisen/` reproducing the DFT Claisen KIEs
   with an MLIP and stating the discrepancy, plus a notebook cell showing
   the projected versus unprojected difference.

**Exit criteria:** `compute_kie` accepts `HessianInput` objects from all
three backends; the ASE path reproduces the Gaussian golden KIE when fed the
Gaussian Hessian; documentation lists the MLIPs tried and their results.

## Phase 9 — Experimental validation suite (v2.5.x, scaffold in place 2026-09-25)

Scaffolded: `benchmarks/` holds `run.py` (computes every case through the
API and writes `REPORT.md`, optionally `report.json`), the case-file format
(`benchmarks/README.md`) and two cases built from the structures in the
repository (Claisen, Diels–Alder) with **experimental values left null**:
they must be entered from the papers, not from memory, and the report says
so until they are. Items 1 (further reactions) and 2 (their frequency jobs)
need a maintainer; item 3 is done apart from the plot and the CI job, which
wait for the first measured values. **Update 2026-09-26:** the first
measured values are in. Two Baeyer–Villiger cases (intermolecular addition
and intramolecular migration) take them from Singleton and Szymanski, JACS
1999, 121, 9455, Figure 1, and from the Crow et al. preprint, ChemRxiv 2026,
SI Tables S3a, S3b and S8b, both read from the PDFs. Their structures are
still needed; Rzepa's blog (post 14112) links transition structures, but
those files are no longer available (September 2026). The runner now
handles cases without structures, replicate measurements, per-KIE reference
positions and several sources. The Diels–Alder values came from the 1995
paper (Figure 1b): nine positions, mean absolute deviation 0.003, the first
measured comparison. The Claisen values came from Meyer et al. 1999
(Table 4): mean absolute deviation 0.0009 at B3LYP and 0.0089 at GFN2-xTB. **Update 2026-09-27:**
the Shi epoxidation case (Singleton, Wang, JACS 2005) is computed. It uses
Gaussian 16 frequency jobs at the SI geometries of the alkene and
transition structure 10, run by the maintainer; they reproduce the SI's
energies and zero-point energies. Kinisot matches the paper's six QUIVER
predictions to the three decimals given, and deviates from experiment by
0.0012 on average. Overall: 25 measured positions, mean absolute deviation
0.0034. The transition structure's Hessian needs about 13 hours with PySCF
on this environment's four cores, against 10 minutes for Gaussian on a
16-core node.
**Update 2026-09-27 (PyQuiverHS SI):** five more cases come from the SI of
Grazioli et al.: three conformational KIEs, a gas-phase SN2 α-secondary KIE
and a CD₃ axial/equatorial EQE. They are included with the authors'
agreement, together with PyQuiverHS's outputs for the same files, and
`tests/test_pyquiverhs.py` cross-checks Kinisot against those outputs.
- **Against experiment:** six of the seven measured values within 0.016;
  the SN2 KIE is 0.08 too high (harmonic TST at HF/6-31+G(d)).
- **Overall:** 11 cases, 32 measured positions, mean absolute deviation
  0.0060.
- **Findings along the way:**
  - PyQuiverHS's enthalpy–entropy terms omit translation (docs/theory.md,
    section 6);
  - a false self-check warning for multi-block Gaussian logs;
  - the EQE direction was stated backwards in the docs.

**Goal (requested 2026-09-25):** a `benchmarks/` directory that compares
Kinisot's predictions with published experimental KIEs, primarily the
natural-abundance NMR measurements of the Singleton group, so that every
release can state how well it reproduces experiment and so that users have
a template for their own comparisons.

1. **Curate the set.** Start with reactions whose structures are already in
   hand: the Claisen rearrangement (Meyer, DelMonte, Singleton, JACS 1999,
   121, 10865: ¹³C, ²H and ¹⁷O KIEs at 393 K) and the isoprene + maleic
   anhydride Diels–Alder reaction (Singleton, Thomas, JACS 1995, 117, 9357).
   Add 6–10 further cases spanning KIE types: an SN2 reaction, an
   epoxidation or dihydroxylation, an ene reaction, a hydride transfer with
   a large primary ²H KIE, and at least one EQE. Each case needs the
   experimental values with uncertainties, temperature, and the reference,
   verified against the paper (not transcribed from memory).
2. **Compute the structures** at a documented level of theory (new Gaussian
   or ORCA frequency jobs; these cannot be produced inside this repository
   and need cluster time). Store the trimmed outputs under
   `benchmarks/<case>/` with the input files, and a `case.yaml` describing
   atom mapping, temperature, scaling and the experimental numbers.
3. **Runner and report**: `benchmarks/run.py` computes every case through
   the API and writes `benchmarks/REPORT.md` with computed vs experimental
   tables (semiclassical, Bell, Wigner), mean absolute deviation per case
   and overall, and a plot of predicted vs measured. Run it in CI as an
   optional job so drift is caught.
4. **Use it**: the report feeds the README ("agreement with experiment"),
   decides the Phase 7 projection default, and provides the reference
   numbers for validating the ORCA (Phase 5) and MLIP (Phase 8) backends
   against the Gaussian results.

**Exit criteria:** ≥ 8 reactions, every experimental number traceable to a
DOI, REPORT.md regenerated by one command.

## Phase 10 — Conformer ensembles (planned 2026-09-27; implemented 2026-09-28)

**Status (2026-09-28).** Conformer ensembles are implemented
(`kinisot/ensemble.py`, `Conformers`, `EnsembleIsotopeEffect`,
`equivalent_positions`; the CLI takes several files per flag with
`--weights`, `--energies` and `--weight-uncertainty`; theory.md section
7). Tests: `tests/test_ensemble.py`, and `tests/test_shi_ensemble.py`,
where `compute_kie` over all 18 Shi transition structures reproduces the
prototype to 2 × 10⁻⁶ at every site. The benchmark runner averages
equivalent positions exactly, which moved two reported values by
1 × 10⁻⁴ (Diels–Alder H3, nitroarene TS1B). Transition structures in
series (`Series`) and parallel channels (`channels`) followed the same day,
with `series_kie`/`channel_kie`, `kinisot --job` (docs/job_files.md) and
theory.md sections 7a and 7b. Tests: `tests/test_pathways.py` (the Wittig
and DyKAT published combinations, a three-step rate matrix to 7 × 10⁻⁸,
channels from one reactant equal to the ensemble to 10⁻¹³) and
`tests/test_jobs.py`. The worked example (`examples/conformers/`, made by
`scripts/make_conformer_example.py`) followed: eight GFN2-xTB reactant
conformers and the chair and boat transition structures, where the C4-d₂
KIE moves by 0.010 from the lowest pair to the ensemble. The benchmark
runner computes a case from a job file (`job` in case.json), so series and
channels enter the report. DyKAT is done that way (2026-09-28): the Gaussian
16 jobs at the paper's geometries reproduce SI Tables 24 and 25 and Figure 3d
within 6 × 10⁻⁴, unscaled as the paper's Kinisot was, and all five combined
KIEs lie within 1.1 standard errors of experiment (`tests/test_dykat.py`).
The Wittig series followed the same day: 4‡ and 6‡ reproduce SI Table S4
within 5 × 10⁻⁴, and weighted by the paper's trajectories the three ¹³C
KIEs match experiment (mean absolute deviation 0.0004), while the
free-energy weighting misses by 0.005 (`tests/test_wittig.py`). Computed
from free 3 instead of the Rh-bound complex, the DyKAT KIEs include the
equilibrium isotope effect of binding (1.008–1.012 at the alkene carbons)
and miss the central carbon by 0.011, so the measurement supports the
paper's model from the complex.

**Goal (requested 2026-09-27):** KIEs and EQEs from several conformers of
each reactant, transition structure and product, given with the same atom
numbering, instead of from one reactant and one transition structure.

**Theory** (a new section of docs/theory.md, before the constants). With
transition-state theory and reactant conformers in equilibrium
(Curtin–Hammett), isotopologue X reacts with

    k^X = (k_B T/h) Σ_j κ_j^X Q‡_j^X exp(−E‡_j/k_B T) / Σ_i Q_i^X exp(−E_i/k_B T)

where i runs over reactant conformers and j over transition-structure
conformers. All conformers of a species have the same atoms, so the
Teller–Redlich mass factor cancels conformer by conformer, and

    KIE = Σ_i x_i ρ_i / Σ_j y_j ρ‡_j

- ρ_i is (s/s')f of reactant conformer i, today's `SideResult.rpfr`.
- ρ‡_j = (s/s')f‡_j · (ν‡_H/ν‡_L)_j · (κ_H/κ_L)_j. With this definition
  ρ_i/ρ‡_j is exactly the pairwise `kie_tunnel` Kinisot computes now.
- x_i is the Boltzmann population of reactant conformer i for the light
  isotopologue, ∝ g_i exp(−G_i/RT), with g_i the conformer's degeneracy.
- y_j is the share of the light isotopologue's rate through transition
  structure j, ∝ g_j κ_j^L exp(−G‡_j/RT).

Four consequences shape the design:

1. **Only populations within each ensemble matter.** The free-energy gap
   between the reactant and transition-structure ensembles cancels.
2. **The result is a ratio of means.** The ensemble KIE is a ratio of
   population-weighted means of partition-function ratios. It is neither
   the Boltzmann average of the pairwise KIEs nor the KIE of the lowest
   pair. A minor transition structure with a different KIE counts in
   proportion to its share of the rate.
3. **Light-isotopologue weights are exact.** No isotope-dependent weights
   are needed, because Σ_i Q_i^H e^{−E_i/kT} / Σ_i Q_i^L e^{−E_i/kT} =
   Σ_i x_i^L (Q_i^H/Q_i^L).
4. **Equivalent positions are permuted conformers.** Positions made
   equivalent by fast motion are one case of the ensemble: the three
   hydrogens of a rotating methyl group, or the two ortho or meta carbons
   of a flipping phenyl ring (`iso_average` and `reference_average` in the
   benchmarks). They are conformers with permuted labels and equal weights.
   - The site KIE is therefore the mean ρ over the labels divided by the
     mean ρ‡.
   - The benchmark runner took the geometric mean of the pairwise KIEs
     until 2026-09-28. On the Diels–Alder methyl reference, the two differ
     by 3.9 × 10⁻⁵.

Two variants follow from the same formula:

- **Semiclassical KIE.** Set κ = 1 in both y and ρ‡.
- **EQE.** Use product conformers in place of transition structures, with
  κ = 1 and no ν‡ ratio.

When the reactant side has several species (a bimolecular reaction), their
ensemble ratios multiply.

Dale, Leach and Lloyd-Jones (J. Am. Chem. Soc. 2021, 143, 21079) propose
the same treatment for flexible systems: Boltzmann-weight the reduced
isotopic partition-function ratios of all low-lying reactant and
transition-structure conformers, then compute a single KIE from the
weighted ratios.

**Transition structures in series (added 2026-09-27).** A reaction may pass
through several transition structures in sequence (R ⇌ I₁ ⇌ I₂ … → P),
none of which alone commits the substrate. Its observed KIE is then a
weighted mean of the KIEs of the individual steps (Dale et al. 2021, eq 25
and SI eqs S69–S81; Paneth, J. Am. Chem. Soc. 1985, 107, 7070). For an
unbranched sequence whose last step is irreversible, the steady-state
approximation makes inverse rate constants add like resistances in series:

    1/k_obs = Σ_n 1/k_n,   k_n = (k_B T/h) exp(−ΔG‡_n/RT)

with ΔG‡_n the free energy of transition structure n above the starting
material. Hence

    KIE_obs = Σ_n w_n KIE_n,   w_n ∝ exp(+ΔG‡_n/RT), for the light isotopologue

- **KIE_n is an ordinary KIE.** It runs from the starting material to
  transition structure n, which Kinisot computes today, so the
  intermediates need no frequencies. For a later step it is the product
  of the equilibrium isotope effects of the steps before it and that
  step's own KIE, as Dale et al. note.
- **Two steps give their eq 25:** KIE_obs = (KIE_TS2 + C_f·KIE_TS1)/(1 + C_f),
  with the commitment factor C_f = k₂/k₋₁ = exp[(ΔG‡₁ − ΔG‡₂)/RT]. A large
  C_f (TS2 below TS1) gives KIE_TS1; a small one gives KIE_TS2.
- **The weights are exact.** As for conformers, weights from the light
  isotopologue alone are exact. A three-step check against the slowest
  eigenvalue of the full rate matrix agrees to 10⁻⁶.
- **It is the opposite average.** The series form is an arithmetic mean
  dominated by the highest transition structure. The parallel form above is
  a harmonic mean dominated by the lowest.
- **The free-energy gap no longer cancels.** The result depends directly on
  ΔΔG‡ between the transition structures, and it matters most when they
  are within about 5 kJ/mol (Dale et al.). The output therefore reports:
  - both limiting KIEs;
  - KIE_obs as a function of ΔΔG‡ (or C_f);
  - the C_f that reproduces a measured value.

  C_f can also be given directly, for example from trajectories.

  The Wittig analysis in Dale et al. (Case Study III) is of this kind:
  the experiment fits C_f ≈ 1, while each transition structure alone
  misses it.
- **Each step can be a conformer ensemble.** Its KIE_n is then the
  ensemble KIE, and its weight uses the ensemble's effective free energy,
  −RT ln Σ_j exp(−G_j/RT).

**Parallel channels with their own labels or reactants (added
2026-09-27).** The formula above assumes one reactant ensemble in
equilibrium and the same label in every structure. Two kinds of
measurement need more.

- **A different label in each channel.** In the Rh-catalysed arylation of
  racemic 3-chlorocyclohexene (van Dijk et al., Nat. Catal. 2021, 4, 284;
  Case Study IV of Dale et al.), the two enantiomers react through
  diastereomeric transition structures, TS-2R and TS-2S. Both converge on
  one η³-allyl intermediate, so the product carbon that is measured came
  from C-1 of one enantiomer and from C-3 of the other.
- **Reactants that do not interconvert.** The two enantiomers are separate
  species in fixed proportion, not conformers in equilibrium.

At low conversion both reduce to

    1/KIE = Σ_c y_c / KIE_c

- KIE_c is the pairwise KIE of channel c, with its own reactant,
  transition structure and label.
- y_c is the channel's share of the light isotopologue's rate. For two
  channels with selectivity s = k_S/k_R this is eq S111 of Dale et al.,
  KIE = (1 + s)·KIE_R·KIE_S/(KIE_S + s·KIE_R). The measured KIEs fit
  s = 3.3.
- **The shares** come from computed free energies, including the
  concentrations of reactants that do not interconvert. They can also be
  given directly, since s is often better known from experiment.
- **Equivalent positions are a special case.** They are one reactant with
  permuted labels (consequence 4). The singly ¹⁸O-labelled nitrobenzene of
  Kang and Radosevich (benchmarks/nitroarene_phosphetane_ts1b) is another:
  its two oxygens are equivalent, and a transition structure that attacks
  one of them makes the observed KIE an average of the attacked and
  spectator positions.
- **Only low conversion.** Shares that drift with conversion, as the faster
  enantiomer is depleted, are left out. The arylation KIEs were measured
  at F_T ≤ 0.17, where that drift is small.

**Populations.** The isotope ratios stay harmonic Bigeleisen–Mayer, because
quasi-harmonic treatments break the product rule (theory.md section 6). The
weights, however, are free energies of the light isotopologue only, and for
floppy conformers a quasi-harmonic treatment of those is appropriate. The
weighting options are:

- `rrho`: electronic energy plus the harmonic free energy from Kinisot's
  scaled frequencies, plus the rotational term from the moments of inertia.
  Translation cancels within an ensemble.
- `qrrho`: Grimme's entropy interpolation at 100 cm⁻¹, as in GoodVibes.
  GoodVibes 4.4 exposes `calc_rrho_entropy`, `calc_freerot_entropy`,
  `calc_damp` and `calc_rotational_entropy`, which take frequencies, so
  Kinisot needs no new thermochemistry code (Guiding decision 1).
- `user`: relative free energies from a table (file, ΔG or G, degeneracy).
  Typical sources are coupled-cluster single points with GoodVibes
  corrections, or a CREST/CENSO ensemble.
- `lowest` and `equal`: for diagnostics.

Degeneracies enter g. Examples are mirror-image conformers in an achiral
environment and symmetry-equivalent rotamers. Rotational symmetry numbers
stay the user's responsibility, as in section 3 of the theory; they could
optionally be detected with `pymsym`, which is already a GoodVibes
dependency.

1. **Theory**: the section above, with the derivation and the four
   consequences.
2. **Input model** (`kinisot/ensemble.py`): a `Conformers(files,
   free_energies=None, degeneracy=None)` wrapper, accepted wherever
   `compute_kie` accepts a species. A nested list (`rct=[["gs1.out",
   "gs2.out"], "b.out"]`) is shorthand for it. Validation:
   - **Same atoms in the same order in every conformer of a species.** A
     mismatch is an error naming the first differing atom, because a label
     must mean the same atom in every conformer.
   - **Same bond graph** (covalent radii). A difference is a warning,
     because a renumbered conformer passes the element check but moves the
     label to a different atom.
   - **Duplicates** found by energy and RMSD after alignment, mirror images
     included. They are a warning, because a duplicate doubles its weight.
   - **One level of theory and one scaling factor** for every file.

   Two further wrappers:
   - `Series([ts_1, ts_2, ...], free_energies=None)`, accepted as `ts`.
     Each member is a file or a `Conformers`.
   - `channels([dict(rct=..., ts=..., iso=...), ...], shares=None)`, a
     function returning the combined result, with the pairwise results
     kept.
3. **Evaluation**:
   - The light isotopologue of every conformer is diagonalized once per
     temperature. That pass gives the weights, ν‡ and κ_L, and is reused
     across labels and for `--reference`.
   - Each label then needs only the heavy diagonalizations.
   - Errors name the conformer's file, for example a transition-structure
     conformer without exactly one imaginary mode, or a minimum with one.
4. **Result**: an `EnsembleIsotopeEffect` with the attributes the CLI, the
   CSV writer and the benchmark runner already read (`kie`, `kie_tunnel`,
   `kie_relative`, `kie_tunnel_relative`, `warnings`, `to_dict`,
   `summary_row`). It adds:
   - **A conformer table**: source, role, ΔG, degeneracy, population or
     rate share, ρ, ν‡, κ_L and κ_H, and the pairwise KIE against the
     reactant ensemble.
   - **Diagnostics**:
     - the KIE of the lowest pair;
     - the effective number of transition structures, 1/Σ y_j²;
     - the range of the KIE when each ΔG moves by ±δ (default
       0.5 kcal/mol).
   - **No ensemble ZPE/EXC/TRPF decomposition.** That breakdown exists per
     conformer only, and the result says so rather than inventing one.
5. **Dispatch in `compute_kie`**:
   - One conformer per species returns today's `IsotopeEffect` unchanged,
     so no golden number moves.
   - Any ensemble returns an `EnsembleIsotopeEffect`.
   - `reference=` and temperature scans apply per label and per temperature.
6. **Tunnelling**: Bell and Wigner apply per transition-structure conformer.
   Skodje–Truhlar needs a barrier per conformer. By default it is measured
   from the lowest reactant conformer, and `--barrier` overrides it for all.
7. **CLI**: `--rct`, `--ts` and `--prd` take one or more files
   (`nargs="+"`). Each flag is one species, and the files after it are its
   conformers. So `--rct a.out --rct b.out` keeps meaning two species, and
   `--ts ts_*.out` is one transition-structure ensemble. New options:
   - `--weights rrho|qrrho|user|lowest|equal`;
   - `--energies table.csv`;
   - `--weight-uncertainty 0.5`.

   The text output prints the conformer table, then `KIE (ensemble) @ T`
   lines. The JSON and CSV outputs carry the table. Series and channels
   need more structure than flags carry, so on the command line they come
   from a JSON job file (`kinisot --job job.json`).
8. **Benchmarks**: `case.json` accepts a list of conformer files per
   species. `iso_average` and `reference_average` move from the geometric
   mean to the exact equal-weight ensemble. That shifts the Diels–Alder
   methyl reference by 4 × 10⁻⁵, in its own commit with a CHANGELOG entry
   (Guiding decision 3).
9. **Tests**:
   - One conformer reproduces `compute_kie` to 1 × 10⁻¹²; the same file
     given twice gives the same KIE.
   - A conformer 10 kcal/mol higher changes nothing.
   - Synthetic two-conformer ensembles match the formula evaluated by hand.
   - Rotamers made by permuting labels reproduce the equivalent-position
     result.
   - A mismatched element order raises an error, and a duplicate triggers
     the warning.
   - `--reference`, temperature scans and EQE ensembles work.
   - **Series.** A kinetic scheme integrated for both isotopologues
     reproduces the series formula. So does the rate-matrix check above,
     with its limits C_f → 0 and C_f → ∞.
   - **Channels.** s → 0 and s → ∞ each return one channel's KIE, and
     eq S111 is reproduced.
10. **Example and validation**:
    - **Worked example** (`examples/conformers/`), cheap enough to
      regenerate in the repository. GFN2-xTB through the ASE backend,
      with reactant conformers of allyl vinyl ether and the chair and boat
      Claisen transition structures.
    - **Physical test: the Shi epoxidation.** The SI of Singleton and Wang
      2005 gives 18 transition structures, with geometries and energies, for
      trans-β-methylstyrene with the fructose-derived dioxirane. The paper
      finds only TS 10, 12 and 13 consistent with the measured KIEs, so an
      ensemble KIE over all of them tests the weighting.
    - **What that needs:** a B3LYP frequency job for each of the other 17
      transition structures. TS 10's took 10 minutes with Gaussian on a
      16-core node. The inputs are in `benchmarks/shi_epoxidation/ensemble/`,
      renumbered to TS 10's atom order. That directory also documents the
      repair of a misprint in the SI's coordinates for TS AB.
    - **Update 2026-09-27: the jobs are done.** Their logs are in that
      directory, and `analyze.py` there prototypes the ensemble formula
      (`tests/test_shi_ensemble.py` checks it).
      - All 18 structures reproduce the SI's energies and zero-point
        energies.
      - Kinisot matches the authors' per-structure predictions (SI Table 1)
        exactly for 102 of 108 values, the others within 0.0008.
      - TS 10 carries 85–90% of the rate under all four weightings, so the
        ensemble KIEs stay within 0.0004 of TS 10's.

      The Shi case therefore checks that the ensemble reduces to the
      dominant structure, but it cannot tell the weighting schemes apart.
      A reaction whose transition structures share the rate more evenly is
      still wanted for that.
    - **Channels: the nitroarene deoxygenation** (Kang and Radosevich,
      Tetrahedron 2025; `benchmarks/nitroarene_phosphetane`). ¹⁸O KIEs of
      1.033 for one oxygen and 1.066 for both, now computed from ORCA jobs.
      - In the monotopic TS1B the attacked and spectator oxygens give
        1.0474 and 1.0140. The exact average for the singly labelled
        substrate is their harmonic mean, 1.0304, which the runner now
        reports (its geometric mean gave 1.0305).
      - The paper's 1.0468 is the attacked oxygen alone, which shows what
        leaving out the average costs.
    - **Channels: the DyKAT arylation** (van Dijk et al. 2021;
      `benchmarks/dykat_allyl_arylation`). This is a real 23:77 split
      (s = 3.3), which Shi cannot provide. The measured KIEs (SI Table 4),
      the per-enantiomer KIEs (SI Tables 24 and 25, computed with Kinisot)
      and Gaussian 16 jobs at the paper's geometries are in the case, and
      `channels.json` reproduces the tables and Figure 3d within 6 × 10⁻⁴
      (mean absolute deviation from experiment 0.0030).
      - The channel formula with s = 3.3 reproduces the paper's combined
        KIEs (Figure 3d) to within the rounding of its inputs, for example
        1.0266 at C3 against 1.027.
      - The paper computes each channel from the Rh-bound substrate
        complex. Computed from free 3 instead, the KIEs include the
        complexation equilibrium isotope effect (1.008–1.012 at the alkene
        carbons) and miss the central carbon C2 by 0.011, so the data
        support the paper's choice.
    - **Series: the Wittig reaction** (Chen, Nieves-Quinones, Waas and
      Singleton, J. Am. Chem. Soc. 2014, 136, 13122;
      `benchmarks/wittig_anisaldehyde`). Table 1 already tests the formula
      without structures. From the single-structure KIEs (4‡: 1.043 and
      1.022; 6‡: 1.015 and 0.994), the series formula gives:
      - with the Figure 2 free energies (25.9 and 26.0 kcal/mol at
        340.15 K, so C_f = 0.862): 1.0280 and 1.0070, against the paper's
        1.028 and 1.008;
      - with their trajectory ratio (C_f = 128/76): 1.0326 and 1.0116,
        against 1.033 and 1.012.

      The measured 1.032–1.033 and 1.011 match the trajectory weighting.
      The statistical weights miss because most trajectories pass the
      betaine without equilibrating, so C_f must also be accepted as an
      input. Computing the single-structure KIEs needs the structures
      from the paper's SI.
    - **Series: Baeyer–Villiger.** The two cases (addition, then
      migration) are a second candidate, once their structures are
      computed.

**Decisions (confirmed 2026-09-28, the recommended options)**

1. The default weighting is qRRHO, since conformer ensembles routinely
   contain low modes that harmonic entropies mistreat; RRHO and user free
   energies are options.
2. The CLI takes `nargs="+"` per flag (backward compatible), not a separate
   `--rct-conformers` flag.
3. The benchmarks' equivalent-position averaging moves to the exact form.
   The shifts were 1 × 10⁻⁴ at most in the reported values (Diels–Alder H3
   and nitroarene TS1B), in their own commit with a CHANGELOG entry.
4. Symmetry numbers and degeneracies are supplied by the user; no
   detection.
5. Parallel ensembles first, then series and channels reusing their
   evaluation and result objects. Both were implemented the same day at
   the user's request, so all three are in 2.6.0.
6. Series and channels on the command line come from a JSON job file.

**Out of scope for the first release**:

- **Channel-specific KIEs.** Some experiments measure KIEs in a product, and
  different transition structures lead to different products. Channels
  that converge on one product, like the DyKAT arylation, are in scope
  (above). The Singleton
  group's recovered-starting-material experiments measure total consumption,
  which the sum over all transition structures describes. A product-specific
  measurement would need an ensemble per channel.
- **Slow conformer interconversion** (non-Curtin–Hammett cases).
- **Dynamical effects beyond transition-state theory.**

**Exit criteria:**

- Single-conformer inputs give today's numbers bit for bit.
- The ensemble formula is tested against hand-computed synthetic cases.
- The series and channel formulas are tested against kinetic schemes and
  eq S111.
- There is one worked example and the theory section.
- The Shi ensemble is reported (done with the prototype; the implementation
  must reproduce it).

---

## Suggested sequencing at a glance

| Phase | Deliverable | Release | Risk |
| ----- | ----------- | ------- | ---- |
| 0 | Characterization + frequency/parity tests | — | none |
| 1 | Six confirmed-bug fixes ✅ | v2.0.3 (publish now) | small documented shift |
| 2 | Validation, exceptions, CLI/output fixes ✅ | v2.1.0 | low |
| 3 | pyproject, entry point, README/examples/docs rewrite ✅ | v2.1.0 | low |
| 4 | Module split, result dataclass, Python API, JSON/CSV ✅ | v2.2.0 | medium (goldens) |
| 5 | GoodVibes dependency, Truhlar v5, ORCA backend ✅ | v2.3.0 | scaling changes documented |
| 6 | Automated publishing, Zenodo, changelog | ongoing | none |
| 7 | Isotope table + syntax, projection, tunnelling, reference ✅ | v2.4.0 | ¹⁸O default change documented |
| 8 | ASE / MLIP backend ✅ | v2.5.0 | new code path, optional extra |
| 9 | Experimental validation suite (Singleton KIEs), scaffold ✅ | v2.5.x | needs new QC calculations and transcribed experimental values |
| 10 | Conformer-ensemble KIEs and EQEs, transition structures in series, parallel channels ✅ | v2.6.0 | none for single-conformer inputs; equivalent-position averaging shifted two reported values by 1 × 10⁻⁴ |

Phases 0–1 are done; 2–3 are a few days and unlock adoption; 4 is the
largest single chunk (one reviewed PR per module move); 5 is now small
because GoodVibes carries the parsers; 7 and 8 are scoped per feature.

## Decisions (confirmed 2026-09-25)

1. **GoodVibes is a hard dependency from Phase 5** (`goodvibes>=4.4`). The
   in-repo Gaussian parser is kept only until the parity test proves the
   adapter equal, then deleted.
2. **The bare-index oxygen default changes from ¹⁷O to ¹⁸O in Phase 7**,
   with a one-release warning whenever a bare oxygen index is used and an
   explicit `--iso 5:17O` form for the old behaviour.
3. **Projection stays off by default for Gaussian and ORCA inputs through
   v2.x** (golden-number continuity) and is on by default for the ASE
   backend from its first release; the default flips for all backends in
   v3.0.
4. **`kinisot/Kinisot.py` survives as a deprecation shim through v2.x** and
   is removed in v3.0; new code lives in `cli.py`, `thermo.py`, `hessian.py`
   and `backends/` from Phase 4.
