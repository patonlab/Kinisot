# Kinisot Codebase Audit

**Date:** 2026-07-02
**Version audited:** 2.0.2 (master @ 8d46135)
**Scope:** All Python source (~410 lines), packaging, CI, and tests.

## Summary

Kinisot is a small, scientifically focused package (~410 lines across 4 modules)
that computes kinetic/equilibrium isotope effects from Gaussian frequency
calculations. The core numerics are sound — the single regression test passes
and the Bigeleisen-Mayer implementation agrees with hand calculation per the
README. However, the codebase carries several **confirmed bugs** (one crash
path, one incorrect physical constant), fragile parsing with no error handling,
**dead CI** (Travis), legacy packaging (`setup.py` only, no `pyproject.toml`),
and near-zero test coverage (1 test, happy path only).

Severity legend: 🔴 bug / crash · 🟠 correctness risk · 🟡 robustness / maintainability · 🔵 modernization

---

## 1. Confirmed bugs

### 1.1 🔴 `NameError` crash when the TS has no imaginary frequency
[Kinisot.py:191](kinisot/Kinisot.py#L191) — `compute_isotope_effect()` calls
`log.Fatal(...)`, but `log` is not defined in that scope (it is a local in
`main()`). When a user supplies a TS file without an imaginary frequency —
exactly the case this branch exists to report — the program crashes with
`NameError: name 'log' is not defined` instead of the intended message.
**Verified empirically.**

### 1.2 🔴 Wrong atomic mass unit constant
[Kinisot.py:37](kinisot/Kinisot.py#L37) — `ATOMIC_MASS_UNIT = 1.660468e-27`.
The CODATA value is `1.66053907e-27`; this looks like a transcription typo
(relative error 4.3×10⁻⁵). It propagates into the Hessian unit conversion and
shifts every computed frequency by ~2×10⁻⁵ relative. KIEs are ratios so much
of the error cancels, but the reference values baked into the regression test
were generated with the wrong constant.

### 1.3 🔴 Duplicated, conflicting physical constants block
[Kinisot.py:24-37](kinisot/Kinisot.py#L24-L37) — the constants are defined
twice; the second (lower-precision, and in the AMU case, wrong) block silently
overrides the first (correct, higher-precision) block. The first block should
be kept and the second deleted.

### 1.4 🟠 `is_linear()` misclassifies molecules
[Hess_to_Freq.py:80-91](kinisot/Hess_to_Freq.py#L80-L91) — linearity is
inferred from the string `"prolate symmetric top"` in the Gaussian output.
Non-linear prolate symmetric tops (e.g. CH₃Cl, NH₃) also print this string, so
they would be treated as linear and lose only 5 instead of 6 trans/rot modes —
one contaminating near-zero "frequency" enters the partition function products.
A robust check is a zero (or missing) rotational constant, or counting
near-zero Hessian eigenvalues.

### 1.5 🟠 Substring matching can select the wrong frequency scaling factor
[Kinisot.py:83-89](kinisot/Kinisot.py#L83-L89) — the level of theory is matched
against the Truhlar table with `str.find()` (substring). `CAM-B3LYP/ma-TZVP`
contains `B3LYP/ma-TZVP`, so a CAM-B3LYP job can pick up a plain-B3LYP factor;
the loop also keeps overwriting on multiple hits, so the *last* match in table
order wins rather than the most specific. Matching should be exact (after
normalizing Gaussian's `R`/`U` prefix) with a deliberate fallback.

### 1.6 🟠 Duplicated/incorrect row in scaling-factor table
[vib_scale_factors.py:22](kinisot/vib_scale_factors.py#L22) —
`"M06/maug-cc-pVTZ"` appears twice with different factors (0.982 and 0.971);
the second is almost certainly meant to be `M06-2X/maug-cc-pVTZ` (0.971 matches
the M06-2X family). With last-match-wins (see 1.5), M06/maug-cc-pVTZ jobs
currently get the wrong (M06-2X) factor.

---

## 2. Robustness

### 2.1 🟡 Unhandled parse failures raise `NameError`/`IndexError`
- [Hess_to_Freq.py:27](kinisot/Hess_to_Freq.py#L27) — `d_o_f` is only assigned
  if `NAtoms=` is found; otherwise the code later crashes with `NameError`.
- [Hess_to_Freq.py:78](kinisot/Hess_to_Freq.py#L78) — `level_of_theory()`
  returns `level+"/"+bs` where `bs` is unbound if no `\Freq\` archive line is
  present (e.g. a single-point job or truncated file).
- [Hess_to_Freq.py:24](kinisot/Hess_to_Freq.py#L24) — mass parsing assumes
  `and mass` lines always have the mass at token index 8.
- No validation that the file is a completed Gaussian *frequency* job, that
  `--iso` atom numbers are within range, or that a second imaginary mode is
  absent (a TS with two imaginary frequencies feeds a negative number into
  `np.log`, producing `nan` without warning).

### 2.2 🟡 `sys.exit()` and bare `except` inside library code
- [Hess_to_Freq.py:32](kinisot/Hess_to_Freq.py#L32) — `read_hess()` calls
  `sys.exit()` on parse failure, killing any host process (including pytest or
  a notebook). Library code should raise exceptions; only the CLI should exit.
- [Kinisot.py:13-18](kinisot/Kinisot.py#L13-L18) — a bare `except:` around the
  package imports masks genuine import errors (e.g. a missing dependency inside
  `Hess_to_Freq`) and, combined with `import *`, makes the namespace opaque.

### 2.3 🟡 Resource handling
- Files are opened without context managers and never closed
  ([Hess_to_Freq.py:16](kinisot/Hess_to_Freq.py#L16), 69, 85).
- `Logger` unconditionally writes `Kinisot_output.dat` to the CWD, silently
  overwriting previous results; `main()` never calls `Finalize()`, and
  `Logger.Fatal()` exits without closing the log
  ([Kinisot.py:43-70](kinisot/Kinisot.py#L43-L70)).

### 2.4 🟡 Numerical edge cases
- [Kinisot.py:121](kinisot/Kinisot.py#L121) — `np.log(math.exp(0.5*x))` can
  overflow `math.exp` at very low temperatures; it is algebraically just
  `0.5*x`. Same pattern cost/precision applies in the other factor functions.
- Frequencies are computed but the printed/rounded copies at
  [Kinisot.py:165-166](kinisot/Kinisot.py#L165-L166) are dead code.
- [Hess_to_Freq.py:35](kinisot/Hess_to_Freq.py#L35) — `hess_mat` is populated
  (upper triangle only) but never used; dead code.

### 2.5 🟡 CLI gaps
- `-s` uses `default=False` as a sentinel for a float option, so `-s 0` would
  be treated as "auto-detect" ([Kinisot.py:221](kinisot/Kinisot.py#L221)).
- If both `--ts` and `--prd` are supplied, `--prd` is silently ignored
  ([Kinisot.py:239-240](kinisot/Kinisot.py#L239-L240)).
- The "no substitution in this file" convention is the undocumented magic
  label `'0'`.
- `parse_known_args()` silently swallows typo'd flags.
- No `--version`, no console entry point (users must type `python -m kinisot`).
- Unreachable `sys.exit()` after `log.Fatal()`
  ([Kinisot.py:232-233](kinisot/Kinisot.py#L232-L233)).

---

## 3. Testing

🟡 **Coverage is one happy-path test** ([tests/test_kinisot.py](tests/test_kinisot.py)):
a single Diels-Alder KIE with a fixed scaling factor. Untested: the EQE
(`--prd`) path, tunneling correction branches, `read_hess` isotope
substitution, `level_of_theory`/scaling auto-detection, `is_linear`, the CLI
(`main()`), all error paths, linear molecules, and multiple isotopomer labels.
- Float assertions use `==` on `round(x, 6)` instead of `pytest.approx`,
  making the test brittle to platform-level last-digit differences.
- Test data lives inside the installed package
  (`kinisot/examples/`, ~MBs of Gaussian outputs shipped to every user)
  rather than in `tests/data/`.

## 4. CI and packaging

- 🔵 **Travis CI is dead.** [.travis.yml](.travis.yml) targets travis-ci.com
  (credit-based, effectively defunct for OSS) with Python 3.5-3.8 on trusty/
  xenial images. No CI currently runs. Migrate to GitHub Actions.
- 🔵 **Legacy packaging.** `setup.py`-only metadata; no `pyproject.toml`
  (PEP 517/621). `python_requires='>=3.0'` is inaccurate (README/recipe say
  ≥3.5; f-strings-free code, but NumPy 2.x requires ≥3.9 anyway).
  `setup.cfg` declares `universal=1` (py2 wheel) despite Python 3-only code.
- 🔵 **Version declared in three places** (`setup.py`, `recipe/meta.yaml`,
  `Kinisot.py:__version__`) with no single source of truth.
- 🔵 Python-2 remnants: `from __future__ import print_function`
  ([Kinisot.py:2](kinisot/Kinisot.py#L2)), `from __future__ import
  absolute_import` ([__main__.py](kinisot/__main__.py)).
- 🟡 `kinisot/__init__.py` is empty — no `__version__`, no re-export of the
  public API (`compute_isotope_effect`), forcing `from kinisot import Kinisot`.
- 🟡 Unused imports: `pathlib`, `glob` ([Kinisot.py:7-9](kinisot/Kinisot.py#L7-L9));
  `math` ([Hess_to_Freq.py:6](kinisot/Hess_to_Freq.py#L6)); `numpy` is used in
  `vib_scale_factors.py` only to hold what is naturally a list of dicts/tuples
  (with awkward `bytes.decode()` at every lookup site).

## 5. Style and structure (non-blocking)

- Mixed 3-space/4-space indentation within the same file.
- `calc_rpfr` is a class used as a function (only `__init__` does work);
  module name `Kinisot.py` shadows the package name in a confusing way
  (`from kinisot import Kinisot`).
- Function names like `calc_product_factor` say they return a "product" but
  return a log-sum (intentional, but only explained in a stray comment).
- README documents `--iso 1,2,3` usage but not the multi-reactant label
  convention or the `--prd`/EQE mode; hard-coded isotope pairs
  (²H, ¹³C, ¹⁷O only — no ¹⁸O, ¹⁵N, ³H) are a stated limitation worth listing
  in docs and as a feature gap.

## 6. What is in good shape

- The core Bigeleisen-Mayer math is clean, logarithmic throughout (good for
  numerical stability), and validated against hand calculation.
- The regression test, while lonely, is a genuine end-to-end numerical check.
- Bundled example files make it easy to build a real test matrix.
- MIT licensed, conda-forge recipe in-repo, small dependency surface (numpy).

---

See [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) for the step-by-step
remediation plan.
