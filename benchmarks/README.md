# Benchmarks: computed versus experimental isotope effects

Each directory holds one reaction: a `case.json` describing the input
files, the substitutions, temperature and scaling, and the experimental
values with their source; `run.py` computes every case through the Kinisot
API and writes `REPORT.md` with computed-versus-measured tables. The goal
(implementation plan, Phase 9) is a set of 8–10 reactions, mostly the
natural-abundance NMR measurements of the Singleton group, so that every
release states how well it reproduces experiment.

```
python benchmarks/run.py            # writes benchmarks/REPORT.md
python benchmarks/run.py --json     # also benchmarks/report.json
```

## `case.json`

```json
{
  "name": "Claisen rearrangement of allyl vinyl ether",
  "reference": {"citation": "...", "doi": "10.1021/ja992372h"},
  "temperature": 393.0,
  "scale": 0.961,
  "level_of_theory": "B3LYP/6-31G(d)",
  "tunneling": "bell",
  "project": false,
  "reactants": ["../../tests/data/gaussian/claisen_gs.out"],
  "transition_structure": ["../../tests/data/gaussian/claisen_ts.out"],
  "product": null,
  "reference_isotopologue": null,
  "kies": [
    {"position": "C4", "iso": "4", "experimental": null, "uncertainty": null, "note": "..."}
  ]
}
```

- Paths are relative to the case directory. `iso` is one label (used for
  every file) or a list with one label per file, in Kinisot's `--iso`
  syntax. `reference_isotopologue` is a label whose KIE divides all others
  (natural-abundance NMR experiments report KIEs relative to a position
  assumed to have none); set it when the experiment did that.
- `experimental` is the measured KIE (or EQE) at `temperature`, with the
  same reference convention, and `uncertainty` its stated error. Leave
  them `null` until they have been read from the paper: **never transcribe
  experimental numbers from memory**. The report marks such rows as
  "no experimental value". When a paper reports independent measurements
  separately, give both as lists (`[1.046, 1.051]`, `[0.005, 0.004]`); the
  deviation is taken from their mean.
- Optional: `reference` in a `kies` entry overrides `reference_isotopologue`
  for that entry (13C and 2H KIEs are often measured against different
  positions). `reference_average` gives several labels instead, for a
  reference group whose positions interconvert (the three hydrogens of a
  rotating methyl group); the KIE is divided by their average. `iso_average`
  does the same for the measured position. Both average the isotope ratios
  over the placements on each side (`kinisot.equivalent_positions`), which is
  exact: the placements are conformers of equal weight.
  `reference` at the top may be a list of sources; `notes` is printed under
  the source line.
- A species may be a list of files instead of one file: its conformers,
  weighted as `compute_kie` does (`weights`: `qrrho` by default, or `rrho`,
  `lowest`, `equal`).
- Optional: `imag_cutoff` (cm⁻¹, default 50) for a transition structure whose
  reaction-coordinate frequency is small, as in some conformational
  processes.
- Optional: `alternative_to` (a case name) for a case that computes a
  mechanism the paper rejects, for comparison with the one it assigns. Its
  rows are reported but left out of the overall mean absolute deviation.
- A case whose structures do not exist yet sets `reactants` and
  `transition_structure` to `null` and each `iso` to `null`; the report
  lists its measurements and marks it "Not computed yet".

## Status

| Case | Structures | Experimental values |
| --- | --- | --- |
| claisen | in repo (B3LYP/6-31G(d)) | entered: Meyer, DelMonte, Singleton, JACS 1999, 121, 10865, Table 4 (5 positions; mean absolute deviation 0.0009) |
| claisen_xtb | in repo (GFN2-xTB, `scripts/make_claisen_structures.py`) | same as claisen (mean absolute deviation 0.0089) |
| diels_alder | in repo (B3LYP/6-31G(d)) | entered: Singleton, Thomas, JACS 1995, 117, 9357, Figure 1b (9 positions; mean absolute deviation 0.003) |
| baeyer_villiger | needed (addition of m-CPBA to cyclohexanone) | entered: Singleton, Szymanski, JACS 1999, 121, 9455, Figure 1a; Crow, Hirschi, Clinton, Hirschi, ChemRxiv 2026 (preprint), SI Tables S3a/S3b |
| baeyer_villiger_migration | needed (migration step of the Criegee intermediate) | entered: same papers, 1999 Figure 1c and 2026 SI Table S8b |
| shi_epoxidation | in repo (Gaussian 16 B3LYP/6-31G(d) at the SI geometries of the alkene and transition structure 10) | entered: Singleton, Wang, JACS 2005, 127, 6679, Figure 1 (6 positions; mean absolute deviation 0.0012; the paper's six predictions reproduced) |
| dihydrophenanthrene | in repo (PyQuiverHS SI; B3LYP/6-31G(d,p), scaled 0.97) | as quoted by Grazioli et al. 2026, Table 3 (Mislow et al., JACS 1963/1964): d6, d4, d10 at 315 K (mean absolute deviation 0.0040) |
| biaryl_diketone | in repo (PyQuiverHS SI; B3LYP/6-31G(d,p), scaled 0.97; `imag_cutoff` 30) | as quoted in their Table 4 (Mislow et al., TL 1962; JACS 1964): d8 at 368 K (deviation 0.015) |
| metaparacyclophane | in repo (PyQuiverHS SI; B3LYP/6-31G(d,p), scaled 0.97) | as quoted in their section 3.2 (Sherrod, Boekelheide, JACS 1972): 0.833 ± 0.04 at 308 K (deviation 0.0004) |
| sn2_chloride_methyl_bromide | in repo (PyQuiverHS SI; HF/6-31+G(d)) | as quoted in their Table 7 (Viggiano, Truhlar et al., JACS 1991): 0.81 ± 0.03 at 300 K (deviation 0.08: harmonic TST misses this gas-phase KIE) |
| tetramethylcyclohexane_eie | in repo (PyQuiverHS SI; B3LYP/6-311G(d)) | as quoted in their Table 8 (Anet et al., JACS 1980): EQE 1.042 ± 0.001 at 290.15 K (deviation 0.0003) |
| nitroarene_phosphetane | in repo (ORCA 6.1.0 M06-2X/6-31G(d,p) at the SI geometries; external modes projected) | entered: Kang, Radosevich, Tetrahedron 2025, 186, 134892, Tables 1 and 2: 18O KIEs 1.033 ± 0.003 (one oxygen) and 1.066 ± 0.003 (both); transition structure TS2 (deviations +0.0016 and +0.0049 with Bell tunnelling) |
| nitroarene_phosphetane_ts1b | in repo (the same jobs) | same measurements against the monotopic TS1B the paper rejects (`alternative_to`, left out of the overall mean): −0.0009 and −0.0007 |
| wittig_anisaldehyde | needed (Gaussian inputs at the SI geometries in the case directory: anisaldehyde, ylide, two transition structures in series; M06-2X/6-31+G**/PCM) and Phase 10's series treatment | entered: Chen, Nieves-Quinones, Waas, Singleton, JACS 2014, 136, 13122, Figure 1 (three positions); Table 1's predictions in the notes |
| dykat_allyl_arylation | needed (Gaussian inputs at the paper's geometries in the case directory: free 3, both reactant complexes, both anti-oxidative-addition transition structures) and Phase 10's channel treatment | entered: van Dijk et al., Nat. Catal. 2021, 4, 284, SI Table 4 (five allyl positions); the paper's combined predictions in the notes |

The five cases from the PyQuiverHS SI also keep PyQuiverHS's input and its
outputs for the same files (`pyquiverhs/`), from 10 to 1000 K;
`tests/test_pyquiverhs.py` checks Kinisot against them. Their measured
values are quoted from that paper, which cites the primary literature
listed in each `case.json`.

Further cases (an ene reaction, a hydride transfer with a large primary
²H KIE, another EQE) need new frequency calculations at a documented
level of theory; put the trimmed outputs (the sections listed in
`docs/file_formats.md`) under `benchmarks/<case>/`.
