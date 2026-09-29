# Changelog

Notable changes to Kinisot. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- Real ORCA 6.1.0 test data (`tests/data/orca/`, from the
  `feat/goodvibes-integration` branch): the TT and GG conformers of
  n-pentane at r2SCAN-3c, and the reactant and transition structure of a
  hydrogen-atom transfer at broken-symmetry M06-2X-D3/6-31+G** with SMD
  (1974.9i cm⁻¹). `tests/test_orca.py` pins the conformer EQE and the HAT
  KIEs (primary and secondary ²H, ³H, Wigner, Skodje–Truhlar, and Bell
  above 438 K), checks that the Bell correction is refused below that
  crossover temperature, and checks the scaling factor detected from the
  ORCA input line. Kinisot reproduces the values first computed for these
  files on that branch within 1.5 × 10⁻⁵ (relative).

## [2.6.0] - 2026-09-29

### Added

- Conformer ensembles (IMPLEMENTATION_PLAN.md, Phase 10). A species given
  as several files, `kinisot --rct gs_1.out gs_2.out --ts ts_*.out`, or in
  Python as a list or `Conformers(files, free_energies=None,
  degeneracy=None)`, is an ensemble of conformers with the same atom
  numbering, in fast equilibrium.
  - KIE = Σ x_i ρ_i / Σ y_j ρ‡_j, with the weights of the light isotopologue,
    which is exact; for a transition structure, the weight is its share of
    the rate, including its tunnelling factor. EQEs work the same way.
    docs/theory.md has a new section 7 with the derivation.
  - `--weights` (`weights=`): `qrrho` (default; quasi-harmonic free
    energies from Kinisot's frequencies with Grimme's entropy interpolation,
    through GoodVibes' functions), `rrho`, `user`, `lowest` or `equal`.
    `--energies TABLE` gives your own free energies and degeneracies.
  - The result, `EnsembleIsotopeEffect`, lists each conformer's free
    energy, share and own KIE, and gives the KIE of the lowest conformers,
    the effective number of transition-structure conformers and the range
    when each free energy moves by ±0.5 kcal/mol (`--weight-uncertainty`).
    The text output prints the conformer table and `KIE (ensemble) @ T`
    lines; JSON and CSV carry the same.
  - Conformers of a species must have the same atoms in the same order
    (an error names the first difference). Different bonding, duplicate
    structures and different levels of theory are warnings.
  - One conformer per species gives the ordinary `IsotopeEffect`, unchanged.
  - `kinisot.equivalent_positions()` gives the exact isotope effect of
    positions made equivalent by fast motion from one result per placement
    of the label.
  - Over the 18 Shi epoxidation transition structures, `compute_kie`
    reproduces the prototype ensemble (`benchmarks/shi_epoxidation/ensemble`)
    to 2 × 10⁻⁶ at every site. Benchmark cases may list conformers per
    species.
- Transition structures in series and parallel channels (Phase 10;
  docs/theory.md, sections 7a and 7b).
  - `compute_kie(rct, ts=Series([ts_1, ts_2, ...]), iso=...)` combines
    steps none of which alone commits the substrate:
    KIE = Σ w_n KIE_n with w_n ∝ exp(+G_n/RT)/κ_L,n. The weights come from a
    commitment factor (two steps), the steps' free energies, or computed
    free energies. The result gives each step's KIE and share, C_f, a
    sensitivity range, `kie_at(C_f)` and `commitment_for(measured)`.
  - `channels([dict(rct=..., ts=..., iso=...), ...])` combines parallel
    routes with their own reactants, transition structures or labels:
    1/KIE = Σ y_c/KIE_c. The shares are given (a measured selectivity),
    built from barriers, or computed; `amounts` covers reactants that do not
    interconvert, and `selectivity_for(measured)` inverts two channels.
  - A step or a channel may be a conformer ensemble, and a channel may be a
    series.
  - `series_kie()` and `channel_kie()` apply the formulas to KIEs from
    anywhere. They reproduce the Wittig predictions of Chen et al. (JACS
    2014; 1.028 from the free energies, 1.033 from the trajectories) and
    Figure 3d of van Dijk et al. (Nat. Catal. 2021; s = 3.3).
  - Checks: the series formula matches the slowest eigenvalue of a
    three-step rate matrix to 7 × 10⁻⁸, and channels computed from one
    reactant reproduce the conformer ensemble of their transition
    structures to 10⁻¹³.
- Worked example `examples/conformers`: eight GFN2-xTB conformers of allyl
  vinyl ether and the chair and boat Claisen transition structures, made by
  `scripts/make_conformer_example.py` (conformer search by dihedral
  rotation, duplicates and mirror images recognized, degeneracies from
  chirality, each saddle point checked to connect reactant and product).
  ¹³C and ¹⁷O KIEs hardly depend on the conformer; the C4-d₂ KIE runs from
  0.926 to 0.984, and the ensemble (0.968) differs from the lowest pair by
  0.010. `run_examples.sh` and `tests/test_examples.py` include it.
- `kinisot --job job.json`: a JSON job file with the structures, labels and
  settings, for series, channels, conformer ensembles with free energies,
  and several isotopologues in one run (docs/job_files.md). The results
  file prints the steps or channels with their shares, and `--json` and
  `--csv` carry an `isotopologue` key. Each isotopologue needs its own name.
- The benchmark runner computes a case from a job file (`"job"` in
  case.json, the format of `kinisot --job`), so transition structures in
  series and parallel channels enter the report. Each `kies` entry names its
  isotopologue in the job file.
- `--calc uma[:model]`: Meta's UMA potentials through fairchem-core
  (`pip install fairchem-core`), with the molecule (`omol`) head in double
  precision and without `torch.compile`. `model` is a fairchem pretrained
  name (default `uma-s-1p1`) or a checkpoint file already downloaded, so
  `scripts/make_claisen_structures.py --model-file` works for UMA too. The
  weights are gated on Hugging Face: accept the licence and set `HF_TOKEN`;
  a failed download says so. The charge and spin of the molecule head
  default to a neutral singlet.
- UMA Claisen structures and Hessians (`tests/data/uma/`, made by
  `scripts/make_claisen_structures.py --calc uma`): UMA is the first
  machine-learned potential tested that finds the concerted Claisen
  transition structure (C1–C6 2.19 Å, C4–O3 1.85 Å, 611i cm⁻¹). Relative to
  C5, its ¹³C and ¹⁷O KIEs match experiment with a mean absolute deviation
  of 0.0011, against 0.0009 for B3LYP and 0.0089 for GFN2-xTB
  (`benchmarks/claisen_uma`, `examples/mlip_claisen`). The overall mean
  absolute deviation of the benchmarks moves to 0.0047 (47 positions).
- `--calc sevennet[:model[:task]]`: SevenNet's multi-task models need the
  task, which selects the training data they reproduce, as in
  `sevennet:7net-omni:omol25_low` (OMol25's ωB97M-V molecules). A missing
  or unknown task is an input error that lists the valid ones.
- Six more potentials through the Claisen check (`tests/data/mlip_claisen/`,
  `examples/mlip_claisen`):
  - MACE-OMOL-0, SevenNet-Omni (`omol25_low`), AIMNet2, AIMNet2-rxn,
    ORB-v3 conservative OMol and OrbMol-v2 all find the concerted
    transition structure.
  - Relative to C5, MACE-OMOL-0 and SevenNet-Omni match the measured KIEs
    as well as UMA and B3LYP (mean absolute deviations 0.0009 and 0.0010).
    AIMNet2 and the two ORB models miss by about 0.003, and AIMNet2-rxn by
    0.0044.
  - The ORB models' energies change by up to 24 meV when the molecule is
    rotated, so their Hessians depend on its orientation. Their
    transition structures are also 4–5% too soft, which gives C5 a KIE of
    1.011–1.013 and a normal ²H₂ KIE of 1.12–1.16 at C1.
  - The example now also recommends checking the Hessian: a symmetric EQE
    and the energy of the rotated molecule.
  - These cases are not benchmark cases, so the overall mean absolute
    deviation is unchanged.
- The Wittig case (`benchmarks/wittig_anisaldehyde`) is computed from
  Gaussian 16 frequency jobs at the SI geometries, as two transition
  structures in series (`series.json`). Both transition structures
  reproduce SI Table S4 within 5 × 10⁻⁴ at every position, and the series
  reproduces the paper's weighted predictions (Table 1). Weighted by the
  paper's trajectories (C_f = 128/76), the three ¹³C KIEs match experiment
  (mean absolute deviation 0.0004). Weighted by the free energies of the two
  transition structures, as the paper shows they should not be, the case
  `wittig_anisaldehyde_statistical` misses by 0.005 at both carbons that
  change bonding; it is left out of the overall mean. The overall mean
  absolute deviation moves to 0.0051 (42 positions).
- The DYKAT case (`benchmarks/dykat_allyl_arylation`) is computed from
  Gaussian 16 frequency jobs at the paper's geometries, as two channels in
  `channels.json`. The transition structures and free 3 reproduce SI Table
  20 to the printed digits. The per-enantiomer KIEs reproduce SI Tables 24
  and 25, and the combined KIEs Figure 3d, within 6 × 10⁻⁴
  (`tests/test_dykat.py`), with unscaled frequencies as in the paper. The
  five combined ¹³C KIEs lie within 1.1 standard errors of experiment (mean
  absolute deviation 0.0030). The overall mean absolute deviation of the
  benchmarks moves from 0.0058 (34 positions) to 0.0055 (39).
  `dykat_allyl_arylation_free3` computes the same channels from free 3,
  which adds the equilibrium isotope effect of binding to Rh (1.008–1.012
  at the alkene carbons). It misses the central carbon by 0.011 and is left
  out of the overall mean.
- The results file names every labelled atom in each file, with its
  element and isotope: `Labelled atoms: claisen_gs C1 -> 13C; claisen_ts
  C1 -> 13C`, and a `Reference atoms` line with `--reference`. Ensembles,
  each step of a series and each channel get the same lines. An atom number
  that is off by one but lands on another atom of the same element passes
  every other check; this line makes it visible.
- Five benchmark cases from the SI of the PyQuiverHS paper (Grazioli, Ly,
  Sabetnejad, Mattapalli, Nguyen and O'Leary, ChemRxiv 2026), included with
  the authors' agreement. Each has its Gaussian files, the PyQuiverHS input
  and PyQuiverHS's output for 10–1000 K:
  - `dihydrophenanthrene`, `biaryl_diketone` and `metaparacyclophane`:
    conformational KIEs;
  - `sn2_chloride_methyl_bromide`: a gas-phase α-secondary KIE;
  - `tetramethylcyclohexane_eie`: a CD₃ axial/equatorial EQE.

  Six of the seven measured values are reproduced within 0.016. The SN2 KIE
  is 0.08 too high, a limit of harmonic transition-state theory at
  HF/6-31+G(d) that PyQuiverHS shares. With the nitroarene case below, the
  benchmarks now cover 34 measured positions, with a mean absolute deviation
  of 0.0058.
- `tests/test_pyquiverhs.py` checks Kinisot against those outputs for all
  11 isotopologues at every temperature, uncorrected and with Wigner and
  Bell tunnelling. Every term agrees within 0.01 K / T (3 × 10⁻⁵ at 300 K).
- The benchmark runner accepts `imag_cutoff` for transition structures with
  a small reaction-coordinate frequency (45i cm⁻¹ in the biaryl case).
- Benchmark cases `nitroarene_phosphetane` and `nitroarene_phosphetane_ts1b`
  compare the ¹⁸O KIEs of Kang and Radosevich (Tetrahedron 2025, 186,
  134892), 1.033 ± 0.003 for one oxygen and 1.066 ± 0.003 for both, with
  ORCA 6.1.0 frequency jobs at the SI geometries of nitrobenzene and the two
  candidate transition structures.
  - Kinisot reproduces the paper's PyQuiver predictions to within 0.001.
  - For the rejected monotopic TS1B, the paper's singly labelled value
    (1.0468) is the attacked oxygen alone. Averaged over the two equivalent
    oxygens it is 1.0304, within the measured value. The ¹⁸O KIEs therefore
    do not distinguish TS1B from TS2; its energy does.
  - `tests/test_nitroarene.py` checks these numbers.
- The benchmark runner accepts `alternative_to` for a case computing a
  mechanism the paper rejects. It is reported but left out of the overall
  mean absolute deviation.
- Benchmark case `wittig_anisaldehyde` holds the ¹³C KIEs of Chen,
  Nieves-Quinones, Waas and Singleton (J. Am. Chem. Soc. 2014, 136, 13122)
  for a Wittig reaction with two transition structures in series. The
  Phase 10 series formula reproduces the paper's two weighted predictions
  from its single-structure KIEs.
- Benchmark case `dykat_allyl_arylation` holds the ¹³C KIEs of van Dijk et
  al. (Nat. Catal. 2021, 4, 284) for the Rh-catalysed arylation of racemic
  3-chlorocyclohexene. Both enantiomers react, through different
  transition structures, and converge on one product. The case includes
  Gaussian frequency inputs at the paper's geometries for the five
  structures involved. The Phase 10
  channel formula reproduces the paper's combined KIEs from its
  per-enantiomer values.
- IMPLEMENTATION_PLAN.md, Phase 10: transition structures in series
  (commitment factors) and parallel channels with their own labels or
  reactants, following Dale, Leach and Lloyd-Jones (J. Am. Chem. Soc. 2021,
  143, 21079).
- docs/theory.md, section 6: free energies and enthalpy–entropy partitions
  must keep the translational terms, which do not cancel when an unlabelled
  partner (the chloride of an SN2 reaction, which Bigeleisen–Mayer lets you
  leave out) adds its mass to the transition structure. PyQuiverHS's
  enthalpy–entropy KIE for Cl⁻ + CH₃Br lacks them and is 1.3% low (0.877
  against 0.888 at 300 K). The same point is in the README and
  docs/comparison.md.
- docs/comparison.md lists QUIVER, THERMISTP and ISOEFF. The README cites
  Rzepa's 2015 comparison of the two routes for the Baeyer–Villiger
  reaction: 1.023 against 1.0226 for ¹³C.
- The Baeyer–Villiger cases cite Rzepa's ωB97XD/Def2-TZVPP models (data
  DOIs 10.14469/ch/1913xx), whose files are no longer available, so the
  structures have to be computed afresh. They also record Singleton's rule
  for which reactant to compute from.
- README, "Why the Bigeleisen–Mayer equation rather than free energies", and
  docs/theory.md section 6: the two routes differ exactly by the
  Teller–Redlich product-rule violation of the computed frequencies. A soft
  mode drops out of Bigeleisen–Mayer but not out of a free-energy
  difference. `scripts/compare_free_energy_route.py` quantifies this on the
  Claisen Hessians: printed-precision rounding, a 0.1 cm⁻¹ soft-mode error,
  and quasi-harmonic free energies. `tests/test_free_energy_route.py` pins
  the quoted numbers.
- Baeyer–Villiger benchmarks (`benchmarks/baeyer_villiger`,
  `benchmarks/baeyer_villiger_migration`) with measured intermolecular and
  intramolecular ¹³C and ²H KIEs. They come from Singleton and Szymanski,
  JACS 1999, 121, 9455 (Figure 1), and from the Crow, Hirschi, Clinton and
  Hirschi preprint, ChemRxiv 2026 (SI Tables S3a, S3b and S8b). Structures
  still need to be computed.
- The Claisen cases hold the measured ¹³C and ¹⁷O KIEs of Meyer, DelMonte
  and Singleton, JACS 1999, 121, 10865 (Table 4, 120 °C, relative to C5).
  The mean absolute deviation is 0.0009 over five positions for the B3LYP
  structures and 0.0089 for GFN2-xTB, so experiment favours the B3LYP
  transition structure.
- The first measured benchmark: the Diels–Alder case now holds the nine
  ¹³C and ²H KIEs of Singleton and Thomas, JACS 1995, 117, 9357 (Figure 1b),
  relative to the methyl group as measured. The B3LYP structures in the
  repository reproduce them with a mean absolute deviation of 0.003. Only
  the inside hydrogen on C1 is off by more than the experimental error;
  Beno, Houk and Singleton (JACS 1996, 118, 9984) likewise found the inside
  hydrogens the only misses.
- A Shi epoxidation case (Singleton and Wang, JACS 2005, 127, 6679,
  Figure 1: ¹³C KIEs at 0 °C). It holds Gaussian 16 B3LYP/6-31G(d)
  frequency jobs at the SI geometries of trans-β-methylstyrene and
  transition structure 10, with their inputs. Those jobs reproduce the SI's
  energies and zero-point energies.
  - Kinisot matches all six of the paper's QUIVER predictions to the three
    decimals given.
  - The mean absolute deviation from experiment is 0.0012. The largest miss
    is the methyl carbon (0.998 predicted, 1.001 and 1.002 measured), which
    the paper's prediction shares.
- `benchmarks/shi_epoxidation/ensemble`: frequency jobs for the other 17
  transition structures in the Shi SI. They are renumbered to TS 10's atom
  order, and TS AB includes the repair of a misprint in the SI's
  coordinates. `analyze.py` there prototypes the ensemble KIE of
  IMPLEMENTATION_PLAN.md, Phase 10.
  - Kinisot reproduces the authors' per-structure predictions for all 18
    structures: 102 of 108 values exactly at three decimals, the rest
    within 0.0008.
  - TS 10 carries 85–90% of the rate, so the ensemble KIEs differ from
    TS 10's by at most 0.0004.
- The benchmark runner lists cases without structures. It also accepts
  replicate measurements, a reference position per KIE, an averaged
  reference (a rotating methyl group), several sources and a notes
  paragraph.

- `--calc xtb[:method]`: GFN2-xTB (or GFN1-xTB) through `tblite`'s ASE
  calculator, with the SCF tightened to `accuracy=0.01`. With tblite's
  default the force noise shifts finite-difference isotope effects by up to
  a few 10⁻⁴ (an EQE between water's equivalent hydrogens came out 1.0001 or
  0.9996 instead of 1); documented in `docs/file_formats.md`.
- GFN2-xTB Claisen structures and Hessians (`tests/data/xtb/`, made by
  `scripts/make_claisen_structures.py --calc xtb`: reactant minimized with
  BFGS, transition structure refined with Sella) and a B3LYP versus GFN2-xTB
  comparison in `examples/mlip_claisen/`; a matching `benchmarks/claisen_xtb`
  case.
- `scripts/make_claisen_structures.py` takes any `--calc` and validates what
  it finds: the reactant must be a minimum, and the saddle point must have
  one imaginary mode, both the C1–C6 and C4–O3 partial bonds, and those two
  stretches dominating the imaginary mode. Otherwise it writes nothing and
  exits with status 1.
- MACE results in `examples/mlip_claisen/`. None of MACE-OFF23 (small,
  medium, large) or MACE-MP-0 (medium) has the concerted Claisen transition
  structure. Their saddle points describe C–O cleavage or a ring closure
  with C–O intact, so no MACE KIEs are reported. The rejected MACE-MP-0
  structures (`tests/data/mace_mp0_rejected/`) are the regression test for
  the checks.
- `tests/test_pyquiver.py`: cross-validation against PyQuiver on the same
  Hessians (Claisen, all positions; Diels–Alder C15/C19): uncorrected, Bell
  and Wigner KIEs agree to 2.3 × 10⁻⁶. `pyquiver-kie` joins the `test` extra.

### Fixed

- The README offered `conda install -c conda-forge kinisot`, but Kinisot is
  not on conda-forge. Install with pip. `recipe/meta.yaml` is refreshed to
  2.5.0 as the recipe to submit to conda-forge's staged-recipes, and
  CONTRIBUTING.md no longer says a conda-forge bot follows each release.
- A Gaussian log with more than one block of printed frequencies, such as
  an `opt=(calcall,ts) freq` job, gave a false "frequencies differ from the
  program's" warning. The self-check took the last 3N−5 printed values,
  one of them from the earlier block. It now takes 3N−6, or 3N−5 for a
  linear molecule. The Hessian and the KIEs were not affected.
- docs/theory.md, section 3, gave the exchange equilibrium behind an EQE in
  the wrong direction. EQE = (s/s')f_R / (s/s')f_P is the constant of
  R(light) + P(heavy) ⇌ R(heavy) + P(light): above 1, the heavy isotope
  accumulates in R. The eqe_cyclohexane example's reading was reversed to
  match. Its 1.038 means CD₃ prefers the axial methyl group, as measured,
  not the equatorial one. The numbers were right.
- The Diels–Alder example and benchmark called TS atom 19 (diene atom 10) a
  diene terminus, and concluded from its KIE of 1.001 that the transition
  structure is markedly asynchronous. Atom 19 is the methyl carbon. The
  termini are C1 (TS 15, 1.022) and C4 (TS 13, 1.018), a moderately
  asynchronous transition structure, as measured.
- `--calc orb` returned orb-models' network instead of an ASE calculator and
  could not compute anything. It now wraps the network in `ORBCalculator`
  (orb-models 0.5 and 0.6+ APIs), in float64, and accepts any
  `ORB_PRETRAINED_MODELS` name as `orb:<model>`.
- A Hessian cached next to a geometry (`--calc`) was reused after
  `--delta` changed. The cache now records the finite-difference step and
  the number of displacements and recomputes when either differs; caches
  written by 2.5.0 are recomputed once.
- `--calc orb` with orb-models' molecular models (OrbMol,
  `orb-v3-*-omol`) failed: they need the total charge and spin multiplicity
  in `atoms.info`. Kinisot now sets a neutral singlet unless the atoms carry
  their own. Charges and spins read back from extended XYZ files are numpy
  integers, which orb-models rejects, so they are converted (for UMA too).

### Changed

- `--calc aimnet2[:model]` uses the `aimnet` package (AIMNet2's current
  distribution) and falls back to the older `aimnet2calc`. `model` is a
  name from aimnet's registry, such as `aimnet2-rxn`, or a model file. A
  model that cannot be loaded is an input error. The fallback applies only
  when aimnet is not installed: an installed aimnet that fails to import
  (a missing submodule, or a missing or incompatible dependency) is
  reported with its own error rather than as a missing `aimnet2calc`.
- The calculator registry records which keyword receives the `:model` part
  of a `--calc` specification (`method` for tblite, `model` elsewhere).
- The benchmark runner averages equivalent positions (`iso_average`,
  `reference_average`) exactly: the isotope ratios are averaged over the
  placements of the label on each side, which for a difference on one side
  is the harmonic mean of the separate KIEs. It took the geometric mean
  before. Two reported values move: the Diels–Alder H3 KIE by 1 × 10⁻⁴
  (0.9904 to 0.9905) and the singly labelled TS1B nitroarene KIE by
  1 × 10⁻⁴ (1.0305 to 1.0304; 1.0321 to 1.0319 with Bell tunnelling). The
  mean absolute deviation over the 34 measured positions stays 0.0058.

## [2.5.0] - 2026-09-25

The first release since 2.0.2. It contains everything developed as the
milestones 2.0.3, 2.1.0, 2.2.0, 2.3.0 and 2.4.0 below, none of which was
published separately; read those sections for the details. In short:

- **Correctness**: the atomic-mass-unit typo and duplicated constants
  removed, CODATA 2018 constants, AME 2020 isotope masses, exact
  scaling-factor lookup, linearity from the geometry, and every frequency
  checked against what the program printed.
- **Robustness**: validated input (atom numbers, elements, already
  substituted atoms, structures that are not minima or transition
  structures, mismatched substitutions), errors instead of crashes, a
  results file that is appended to rather than overwritten.
- **Inputs**: Gaussian, ORCA (`.out` + `.hess`) and ASE (`VibrationsData`
  JSON, or any geometry with `--calc` and an ASE calculator, machine-learned
  potentials included), mixable in one run.
- **Physics options**: explicit isotope syntax with a full table (`5:13C`,
  `3:17O`, `7:D`), a bare oxygen index now meaning ¹⁸O, Eckart projection of
  external modes (default for ASE Hessians), Bell, Wigner and
  Skodje–Truhlar tunnelling, a reference isotopologue, harmonic or
  fundamental Truhlar factors, temperature scans.
- **Interfaces**: the `kinisot` console script, `--json`/`--csv` output,
  and a Python API (`compute_kie` returning an `IsotopeEffect`); the 2.0
  functions remain as deprecated shims until 3.0.
- **Project**: `pyproject.toml`, GoodVibes ≥ 4.4 as a dependency, worked
  examples with verified outputs, documentation, a benchmark scaffold, CI
  on three platforms, automated PyPI and GitHub releases.

**Numerical changes relative to 2.0.2** (all documented in the milestone
sections): the constants fix (< 3 × 10⁻⁶ in KIEs), CODATA 2018
(7 × 10⁻⁹), AME 2020 isotope masses (≤ 1.2 × 10⁻⁶), the Truhlar v5
scaling table (a few factors), and ¹⁸O instead of ¹⁷O for a bare oxygen
index (`3:17O` restores the old value).

### Fixed in the release preparation

- The wheel now includes the `kinisot.backends` subpackage (package
  discovery was limited to the top-level package, so an installed copy could
  not be imported; the in-tree test runs did not notice). Verified by
  installing the wheel and running the suite from outside the repository on
  Python 3.9, 3.10, 3.12 and 3.13.

### Phase 8 details (ASE and machine-learned potentials)


- **ASE backend** (`pip install kinisot[ase]`): Hessians from
  `VibrationsData` JSON files, or computed from a geometry with any ASE
  calculator (`--calc mace_mp:medium`, `emt`, `orb`, `sevennet`, `aimnet2`,
  `module:callable`; `--delta` for the finite-difference step, analytic
  `get_hessian` used when available) and cached next to the geometry as
  `<name>.hessian.json`. `kinisot.hessian_from_calculator`,
  `hessian_for_geometry`, `save_hessian_json`, `build_calculator` in the API;
  `compute_kie(..., calculator=)`.
- Projection of external modes is now **on by default for ASE inputs** and
  off for Gaussian/ORCA (`project=None` means "decide by backend";
  `--project`/`--no-project` override). With projection, imaginary modes
  below the cutoff are treated as real vibrations of the same magnitude,
  with a warning.
- `examples/mlip_claisen/` and ASE JSON fixtures (`tests/data/ase/`); the
  Claisen KIE from the JSON form matches the Gaussian path to 10⁻¹⁰.
- `benchmarks/`: a runner and case-file format for comparing computed with
  experimental isotope effects, with the Claisen and Diels–Alder cases
  (experimental values to be entered from the papers).
- The release workflow now also creates a GitHub Release with the changelog
  section as notes; Dependabot watches the Actions versions; `.zenodo.json`
  carries the archive metadata.

## [2.4.0] - development milestone (included in 2.5.0)

Phase 7 of the [implementation plan](IMPLEMENTATION_PLAN.md): isotopes,
external modes, tunnelling.

### Added

- **Isotope table** (`kinisot/isotope_data.py`, generated from the
  `periodictable` package: AME 2020 masses for every naturally occurring
  isotope of 83 elements plus ³H, ¹¹C, ¹⁴C, ¹³N, ¹⁵O, ¹⁸F, ³²P, ³³P, ³⁵S,
  ³⁶Cl, ¹²⁵I, ¹³¹I). Both isotopologues are built from it, so Gaussian, ORCA
  and any other program give the same numbers for the same Hessian.
- **Explicit isotope syntax**: `--iso 5:13C`, `3:17O`, `7:D`, `7:T`, or an
  explicit mass `5:13.5`; a bare number means the default heavy label
  (H, C, N, O, S, Cl, Br, Si). Labels that do not match the atom's element
  are refused.
- **Eckart projection** of translations and rotations (`--project`,
  `project=True`); external modes are removed by value and the reaction
  coordinate is the most negative remaining mode. Needs the geometry, which
  the Gaussian (archive) and ORCA (`$atoms`) readers now provide.
- **Skodje–Truhlar tunnelling** (`--tunneling skodje`, barrier from
  `--barrier KCAL` or from the electronic energies now read from the files:
  Gaussian archive `HF=`, ORCA `FINAL SINGLE POINT ENERGY`).
- **Reference isotopologue** (`--reference ATOMS`, `reference=`): the KIE
  is also reported relative to a second substitution.
- **Temperature scans**: `-t 273,298,323` or `-t 250:350:10` give one
  result line per temperature (and one JSON entry / CSV row each).

### Changed

- **Bare oxygen index means ¹⁸O** (was ¹⁷O); a note is printed once per
  file in this release. `3:17O` gives the old behaviour.
- **Isotope masses** moved from the five-decimal values of Kinisot 1.x/2.x
  (¹H 1.00783, ²H 2.0141, ¹³C 13.00335, ¹⁷O 16.9991) to AME 2020. Largest
  change in a golden KIE: 1.2 × 10⁻⁶ relative (the H7,H8 Claisen case);
  a few bundled results move in the sixth decimal. Golden tests updated.
- Projection versus the lowest-six rule on the bundled Gaussian examples:
  the corrected KIEs differ by 5 × 10⁻⁸ (Claisen C4), 4 × 10⁻⁷ (Claisen
  H7,H8), 2 × 10⁻⁸ (Diels–Alder C15) and 3 × 10⁻⁸ (CD₃ EQE). Projection
  stays off by default for Gaussian/ORCA input through 2.x, as decided.
- The `substitute()` API returns AME masses and records the isotope in each
  `Substitution`; `parse_label()` returns `(index, symbol, mass number,
  mass)` tuples.

## [2.3.0] - development milestone (included in 2.5.0)

Phase 5 of the [implementation plan](IMPLEMENTATION_PLAN.md): GoodVibes
integration and ORCA support.

### Added

- **ORCA input**: give `name.out` or `name.hess` (ORCA writes the Hessian
  to the latter); the level of theory is read from the `!` line. The
  program is detected from the file, so Gaussian and ORCA files can be
  mixed. ORCA's standard atomic weights are replaced by pure-isotope masses
  for the elements Kinisot can substitute, so the same Hessian gives the
  same numbers from either program.
- **Frequency self-check**: for every unsubstituted species Kinisot compares
  the frequencies it obtains from the Hessian with the ones the program
  printed and warns when they differ by more than 1 cm⁻¹.
- `--scale-type zpe|harm|fund` (and `scale_type=` in the API) to apply the
  harmonic or fundamental Truhlar factor instead of the ZPE one.
- Geometry-based linearity test (`kinisot.hessian.linear_from_geometry`),
  used for ORCA and available to every backend; `HessianInput` gains
  `program_frequencies`.
- `examples/orca_claisen/` and ORCA-layout test fixtures (`tests/data/orca/`).

### Changed

- **GoodVibes ≥ 4.4 is a dependency.** Scaling factors now come from the
  Truhlar database version 5 shipped with GoodVibes, with its cross-program
  canonicalization (Gaussian `PBE1PBE` = ORCA `PBE0`, `6-31G*` = `6-31G(d)`,
  ...). Differences from the version 3b2 table Kinisot carried before:
  M06-2X/6-31+G(d,p) 0.967 → 0.968; `M06-2X/maug-cc-pVTZ` and
  `PW6B95/6-31+G(d,p)` are no longer listed (factor 1.0 with a message);
  four spellings (`MN15-L`, `MN12-L`, `MN12-SX`, `M06-L(DKH2)`) are aliased
  by Kinisot until GoodVibes carries them. The bundled examples are
  unaffected (B3LYP/6-31G(d) → 0.977 in both).
- `kinisot/vib_scale_factors.py` is gone; `kinisot.scaling.find_scaling_factor`
  takes an optional scale type.

## [2.2.0] - development milestone (included in 2.5.0)

Phase 4 of the [implementation plan](IMPLEMENTATION_PLAN.md): internal
refactor and Python API. Numerically identical to 2.1.0 apart from the
constants update listed under Changed.

### Added

- `kinisot.compute_kie()` returning a frozen `IsotopeEffect` (final factors,
  per-side `SideResult`s with the light and heavy isotopologues, their kept,
  discarded and imaginary frequencies, masses and substitutions, the scaling
  choice and any warnings) with `to_dict()`, `to_json()` and
  `summary_row()`. Inputs may be file paths or `HessianInput` objects.
- `--json FILE` (full result of a run) and `--csv FILE` (one row per run).
- `--tunneling bell|wigner|none` on the command line and `tunneling=` in the
  API. The Bell correction is refused below the crossover temperature,
  where it diverges, instead of returning a meaningless number.
- `HessianInput`, the program-independent interchange type every backend
  produces (Hessian, masses, atomic numbers, level of theory, linearity,
  positions), in `kinisot/hessian.py`; `kinisot/backends/gaussian.py`
  holds the Gaussian parser (now also reads the archive geometry).
- `examples/api_example.py` and an executed `examples/examples.ipynb`.

### Changed

- Module layout: physics in `kinisot/thermo.py`, scaling in
  `kinisot/scaling.py`, isotopes in `kinisot/isotopes.py`, the CLI in
  `kinisot/cli.py`, the API in `kinisot/api.py`. `kinisot.Kinisot` and
  `kinisot.Hess_to_Freq` remain as deprecation shims (removed in 3.0):
  `compute_isotope_effect()`, `calc_rpfr`, `get_frequency_scaling()` and
  the `calc_*_factor` functions warn with `DeprecationWarning` and forward
  to the new code.
- Physical constants updated from CODATA 2010 to CODATA 2018 (h, k, E_h,
  a_0, u). Largest change in any golden quantity: 7 × 10⁻⁹ relative; no
  printed digit of the bundled examples changes.
- The scaling-factor table is a plain dictionary of named tuples
  (`kinisot.vib_scale_factors.SCALING_FACTORS`) instead of a NumPy
  structured array; factors are now exact decimals rather than float32.

## [2.1.0] - development milestone (included in 2.5.0)

Phases 2 (robustness) and 3 (packaging and documentation) of the
[implementation plan](IMPLEMENTATION_PLAN.md). Every result line of the
bundled examples is unchanged (to the 6 printed decimals, modulo the 2.0.3
constants fix).

### Added

- `kinisot` console script (`python -m kinisot` still works) and a
  `pyproject.toml` (PEP 621) build; `setup.py`/`setup.cfg` are gone. The
  version is defined once, in `kinisot/__init__.py`, which now exports the
  public API (`compute_isotope_effect`, `parse_gaussian`, the exception
  classes, ...).
- Worked examples in `examples/` (Claisen rearrangement, Diels–Alder with
  one or two reactant files, conformational EQE) with the commands, the
  expected output and literature background; `examples/run_examples.sh`
  regenerates the outputs and `tests/test_examples.py` checks them.
- `docs/theory.md` (the equations as implemented), `docs/file_formats.md`,
  `docs/faq.md`, `docs/comparison.md`; a rewritten README with a quick
  start, a guide to every output column and the Python API.
- `CITATION.cff`, `CONTRIBUTING.md`, issue templates, `ruff` lint and
  format checks in CI, and a tag-driven PyPI publishing workflow
  (`.github/workflows/publish.yml`, trusted publishing).

- Input validation with clear messages: `--iso` atom numbers out of range,
  `0` combined with atom numbers, duplicated atoms, atoms of an element
  Kinisot cannot substitute, atoms that already carry a heavy isotope in the
  Gaussian input, a reactant or product with an imaginary frequency, a
  transition-structure side with more than one file carrying an imaginary
  frequency, different substituted elements on the two sides of the reaction
  (e.g. `13C` in the reactant but `2H` in the TS), and labels that request
  no substitution at all. Previously most of these were silently ignored.
- A warning when a transition structure has more than one imaginary
  frequency beyond the cutoff; a negative frequency can no longer reach the
  partition functions unnoticed (`nan` results).
- The results table now lists, for every species, the modes kept in the
  partition function and the modes discarded as external (translation and
  rotation), so a misassigned low-frequency mode is visible.
- `--output FILE` (`-o`), `--overwrite`, `--quiet` (`-q`), `--version`,
  `--imag-cutoff` (the old spelling `--cutoff` still works), and long forms
  `--temperature` and `--scale`.
- `kinisot.exceptions` (`KinisotError`, `KinisotParseError`,
  `KinisotInputError`, `KinisotWarning`); library code raises these instead
  of calling `sys.exit()`. Both error classes also subclass `ValueError`.
- `Hess_to_Freq.parse_gaussian()` returns a `FrequencyData` record (Hessian,
  masses, atomic numbers, level of theory, linearity) from one read of the
  file; `substitute()` reports which atoms were changed.
- Tests: frequencies of every bundled output are checked against the values
  Gaussian prints (kept modes to 0.05 cm-1, discarded modes against
  Gaussian's low frequencies); error paths; CLI behaviour including a
  `python -m kinisot` run. Coverage is reported in CI and must stay at or
  above 80 %.

### Changed

- The Gaussian outputs used by the tests and examples moved from
  `kinisot/examples/` to `tests/data/`; the wheel now weighs 23 kB instead
  of shipping 2.8 MB of test data.
- **Results file**: new results are appended to `Kinisot_output.dat`
  (or `--output`) instead of silently overwriting it, matching how the
  bundled example scripts collect a series of substitutions; use
  `--overwrite` to start afresh. The `Species: ... isotopologue: ...` lines
  are now written to the results file, not only to the terminal.
- **Per-species rows** of the results table now show each species' own
  Bigeleisen-Mayer factors, i.e. light/heavy ratios of the ZPE and
  excitation terms and heavy/light ratio of the frequency product (TRPF).
  Before 2.1.0 the TRPF column of the two rows was swapped (the reactant
  row showed the transition-structure quantity). The final `KIE @` line is
  unchanged and still equals the ratio of the two rows in every column.
- File names in the table are shown without directory and extension; a
  multi-file side is shown as `a + b`. The result line of an equilibrium
  isotope effect is labelled `EQE @` instead of `KIE @`.
- The level-of-theory check covers all files (not only the first two) and
  its warning is printed to the terminal.
- Unknown command-line flags are rejected (`parse_args`), `--ts` together
  with `--prd` is rejected, and a non-positive temperature or scaling
  factor is rejected. `-s` no longer treats `0` as "auto-detect".
- Windows Gaussian archives (`|` separators) are parsed; archive entries
  wrapped across lines are joined before parsing.
- The Bigeleisen-Mayer terms are evaluated with NumPy array operations.

### Removed

- The `try/except` import shim: Kinisot is run as `python -m kinisot`
  (or, from Phase 3, the `kinisot` console script); `python Kinisot.py`
  is no longer supported.
- `Logger.Fatal()`; the logger is now a context manager.

## [2.0.3] - development milestone (included in 2.5.0)

### Fixed

- **Physical constants**: removed a duplicated constants block that silently
  overrode the correct values with lower-precision ones, including a
  transcription error in the atomic mass unit (1.660468e-27 kg instead of
  CODATA 1.660539e-27 kg). Harmonic frequencies shift by ~2×10⁻⁵ relative;
  computed KIE/EQE values shift by <3×10⁻⁶ because the error largely cancels
  in the Bigeleisen-Mayer ratios. Results printed to 6 decimal places are
  unchanged for all bundled examples except in the last digit.
- **Crash on TS without an imaginary frequency**: passing a ground-state file
  as `--ts` crashed with `NameError`; it now reports a clear error message.
- **Vibrational scaling factor auto-detection**: the level of theory is now
  matched exactly against the Truhlar database (after normalizing case,
  hyphens, and Gaussian's R/U/RO spin prefix). Previously substring matching
  could apply a wrong factor (e.g. plain B3LYP factors to CAM-B3LYP jobs) and
  the last matching table row silently won.
- **Scaling-factor table**: the entry mislabeled `M06/maug-cc-pVTZ`
  (zpe_fac 0.971) is now correctly `M06-2X/maug-cc-pVTZ`, verified against
  the Truhlar database v3b2.
- **Linear molecule detection**: linearity is now determined from the
  rotational constants instead of the "prolate symmetric top" string, which
  non-linear prolate tops (e.g. CH3Cl) also print — those molecules would
  have had one contaminating near-zero mode included in the partition
  function products.

### Added

- Characterization test suite (23 tests) covering Claisen KIEs,
  single/multi-reactant Diels-Alder KIEs, the EQE (`--prd`) path,
  scaling-factor lookup, and linearity detection.
- GitHub Actions CI (Python 3.9-3.13 on Linux/macOS/Windows), replacing the
  defunct Travis setup.

## [2.0.2] - 2021

Last release before this changelog was introduced.
