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
  "no experimental value".

## Status

| Case | Structures | Experimental values |
| --- | --- | --- |
| claisen | in repo (B3LYP/6-31G(d)) | to be entered from Meyer, DelMonte, Singleton, JACS 1999, 121, 10865 |
| claisen_xtb | in repo (GFN2-xTB, `scripts/make_claisen_structures.py`) | same as claisen |
| diels_alder | in repo (B3LYP/6-31G(d)) | to be entered from Singleton, Thomas, JACS 1995, 117, 9357 |

Further cases (an SN2 reaction, an epoxidation, an ene reaction, a hydride
transfer, another EQE) need new frequency calculations at a documented
level of theory; put the trimmed outputs (the sections listed in
`docs/file_formats.md`) under `benchmarks/<case>/`.
