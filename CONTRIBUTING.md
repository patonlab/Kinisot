# Contributing to Kinisot

## Set up

```
git clone https://github.com/patonlab/Kinisot
cd Kinisot
pip install -e ".[test,lint,ase]"
```

## Check your change

```
python -m pytest                 # 100+ tests, ~2 s; coverage must stay >= 80 %
ruff check . && ruff format --check .
examples/run_examples.sh         # regenerates examples/*/expected_output.dat
```

CI runs the same commands on Linux, macOS and Windows for Python 3.9–3.13.

## Ground rules

- **Golden numbers move only deliberately.** `tests/test_kinisot.py` pins
  the KIEs of the bundled examples to 1e-6 and `tests/test_frequencies.py`
  pins every frequency to what Gaussian prints. A change that shifts them
  needs its own commit, updated goldens, and a CHANGELOG entry stating the
  size of the shift and why.
- **Library code raises, the CLI exits.** Use the classes in
  `kinisot/exceptions.py`; `sys.exit` is only allowed in `main()`.
- **Every feature comes with a test, a CHANGELOG line and a README or
  docs sentence** in the same pull request.
- **Test data** lives in `tests/data/` and is not shipped in the wheel.
  Keep new fixtures small: trim quantum-chemistry outputs to the sections
  Kinisot reads (see `docs/file_formats.md`), or build synthetic files with
  the helpers in `tests/conftest.py`.

## Regenerating data

`python scripts/make_isotope_data.py` rebuilds `kinisot/isotope_data.py` from
the `periodictable` package (AME 2020 masses); commit the result.
`python scripts/make_xtb_claisen.py` (needs `tblite` and `sella`) rebuilds the
GFN2-xTB Claisen structures and Hessians in `tests/data/xtb/`; if they change,
update the numbers in `examples/mlip_claisen/README.md` and `tests/test_ase.py`.

## Adding a backend (ORCA, ASE, ...)

See `kinisot/backends/orca.py` for a complete example and
[IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md), Phase 8, for ASE: a
backend produces the Hessian (Hartree/Bohr²), the per-atom masses and
atomic numbers, the level of theory and the linearity of the molecule; the
physics in `Kinisot.py` is program independent. Add a fixture pair
(reactant/TS) and a golden KIE for every new backend.

## Reporting a problem

Open an [issue](https://github.com/patonlab/Kinisot/issues) with the
command you ran, the full message Kinisot printed, and, if possible, the
output files (or the sections listed in `docs/file_formats.md`).

## Releasing

Releases are published to PyPI by `.github/workflows/publish.yml` when a
version tag is pushed. PyPI trusted publishing (owner `patonlab`,
repository `Kinisot`, workflow `publish.yml`, environment `pypi`) and the
Zenodo–GitHub connection are configured, so no token or manual upload is
involved. For each release:

```
# bump __version__ in kinisot/__init__.py, move the CHANGELOG "Unreleased"
# section to the new version, update CITATION.cff, commit, then
git tag v2.1.0
git push origin v2.1.0
```

The workflow checks that the tag matches `kinisot.__version__`, runs the
tests, builds the sdist and wheel, uploads them, and then creates a GitHub
Release whose notes are the matching CHANGELOG section. Once the release is
on PyPI the conda-forge bot opens the feedstock update automatically
(`recipe/meta.yaml` in this repository is the template to keep in sync).

**Zenodo DOI per release**: the repository is connected to Zenodo, so
every GitHub Release is archived with the metadata in `.zenodo.json` and
gets its own DOI under the existing concept DOI. Update `CITATION.cff`
(`version`, `date-released`) with each release. Dependabot keeps the
Actions versions in the workflows current (monthly pull requests).
