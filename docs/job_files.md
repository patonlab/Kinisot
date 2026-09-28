# Job files

Transition structures in series, parallel channels and runs over several
isotopologues need more structure than command-line flags carry. A JSON
job file holds the structures, the labels and the settings:

```
kinisot --job job.json [-o FILE] [--overwrite] [-q] [--json FILE] [--csv FILE]
```

With `--job`, the command line gives only the outputs. File names in the
job are relative to the job file. The equations are in
[theory.md, section 7](theory.md#7-conformer-ensembles).

## Settings

Every key is optional and means what the flag of the same name means.

| Key | Flag | Default |
| --- | --- | --- |
| `temperature` | `-t`: a number, a list, or a string such as `"250:350:10"` | 298.15 |
| `scale`, `scale_type` | `-s`, `--scale-type` | Truhlar ZPE factor |
| `imag_cutoff` | `--imag-cutoff` | 50 |
| `tunneling`, `barrier` | `--tunneling`, `--barrier` | `"bell"` |
| `project` | `--project` (`true` or `false`) | on for ASE inputs only |
| `calc`, `delta` | `--calc`, `--delta` | |
| `weights`, `weight_uncertainty` | `--weights`, `--weight-uncertainty` | `"qrrho"`, 0.5 |
| `description` | free text, ignored | |

## Species

A species is one of:

- a file name: `"gs.out"`;
- a list of conformers with the same atom numbering: `["gs_a.out", "gs_b.out"]`;
- an object with free energies (relative, within the species) and
  degeneracies: `{"files": ["a.out", "b.out"], "free_energies": [0.0, 0.6],
  "degeneracy": [1, 2], "energy_unit": "kcal/mol"}`.

`reactants`, `transition_structure` and `product` are lists of species.

## Structures

Exactly one of `transition_structure` (a KIE), `product` (an EQE),
`series` or `channels`. The atom numbers in the examples below are
placeholders.

**Transition structures in series.** The steps in the order the reaction
meets them, and how to weight them: the steps' free energies (of their
lowest conformers, any common zero), a commitment factor (two steps:
C_f = k₂/k₋₁), or neither, when Kinisot computes free energies (every step
must then have the same atoms). `iso` gives one label per reactant and then
one per step.

```json
{
  "temperature": 340.15,
  "reactants": ["anisaldehyde_1.log", "ylide_2.log"],
  "series": {"steps": ["ts_4.log", "ts_6.log"], "free_energies": [25.9, 26.0]},
  "iso": ["8", "0", "20", "20"]
}
```

`"commitment": 1.684` in place of `"free_energies"` weights the steps by
128 trajectories that go on against 76 that return.

**Parallel channels.** Each channel has its own `reactants` and either a
`transition_structure` or a `series`, and optionally a `name` and an
`amount` (the relative amount of its reactant, when the reactants do not
interconvert; 1 by default). The shares of the rate come from `shares`
(for example a measured selectivity), from `barriers` (kcal/mol, or
`energy_unit`), or, with neither, from free energies Kinisot computes
(every channel must then have the same molecules). `iso` has one entry per
channel: that channel's labels.

```json
{
  "temperature": 313.15,
  "channels": [
    {"name": "(S)-3", "reactants": ["allyl_chloride_3_5d.log"], "transition_structure": ["s3_anti_oa_ts.log"]},
    {"name": "(R)-3", "reactants": ["allyl_chloride_3_5d.log"], "transition_structure": ["r3_anti_oa_ts.log"]}
  ],
  "shares": [3.3, 1],
  "iso": [["1", "40"], ["3", "40"]]
}
```

Each channel's labels follow its own files' numbering, so the same product
carbon can be a different atom number in each channel.

## Isotopologues

`iso` and `reference` describe one isotopologue. For several, give a list
instead:

```json
"isotopologues": [
  {"name": "carbonyl C", "iso": ["8", "0", "20", "20"], "reference": ["3", "0", "15", "15"]},
  {"name": "ylide CH", "iso": ["0", "5", "11", "11"]}
]
```

## Output

The results file gets one table per isotopologue: the steps or channels
with their shares and KIEs, then one line per temperature with the
combined KIE, the tunnelling correction, the range when each free energy
or barrier moves by `weight_uncertainty`, and the commitment factor C_f
(two steps) or the selectivity s (two channels). `--json` writes every
result with an `isotopologue` key, `--csv` one row per isotopologue and
temperature.

From Python the same results come from
`compute_kie(rct, ts=Series([...]), iso=...)` and
`channels([dict(rct=..., ts=..., iso=...), ...], shares=...)`. The results
also give the C_f or s that reproduces a measured value
(`result.commitment_for(value)`, `result.selectivity_for(value)`), and
`series_kie()` and `channel_kie()` combine KIEs you already have.
