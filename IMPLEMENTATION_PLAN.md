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

## Phase 1 — Confirmed bug fixes ✅ (v2.0.3, unreleased)

All six items done (constants block, `NameError`, exact scaling-factor
match, table row, linearity from rotational constants, CHANGELOG). Remaining:
publish v2.0.3 to PyPI and bump the conda-forge recipe. Do this before
Phase 2 lands so users get the bug fixes without waiting for the CLI
changes.

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
   The conda-forge bot picks up the PyPI release; keep `recipe/meta.yaml`
   in sync.
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
wait for the first measured values.

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
