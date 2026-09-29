# FAQ

## Setting up a calculation

**How do I find the atom number of a position?** Kinisot numbers atoms from
1, in the order of the input geometry. The rule is the same for every
program, and it is the numbering GaussView shows. Take the numbers from
your input file or from GaussView. ORCA's printed tables count from 0, so
the atom ORCA lists as 11 is `--iso 12`.

**The atom numbering differs between my reactant and my transition
structure.** Give one `--iso` per file: the `--rct` files first, in order,
then the `--ts` or `--prd` file. Kinisot checks that both sides are labelled
with the same isotopes, but it cannot tell two carbons apart, so check
each number against your structures.

**What does `--iso 0` mean?** No labelled atom in that file. It is needed
when a side of the reaction is given as several files and only one of them
carries the label. `0` cannot be combined with atom numbers.

**Which isotope does a bare atom number mean?** The usual heavy label:
²H, ¹³C, ¹⁵N, ¹⁸O, ³⁴S, ³⁷Cl, ⁸¹Br, ²⁹Si. Anything else is asked for
explicitly (`--iso 3:17O`, `7:D`, `7:T`, `5:14C`, or a mass `5:13.5`).
Kinisot 2.3 and earlier used ¹⁷O for a bare oxygen number (the isotope
measured by ¹⁷O NMR in the Claisen study). 2.4 and later use ¹⁸O and print
a note once when they see a bare oxygen number.

**Does the temperature have to match the frequency job?** No. Kinisot
evaluates everything at `-t`; the frequency calculation itself does not
depend on temperature. Use the temperature of the experiment.

**Can I use implicit solvent?** Yes. Kinisot takes the frequencies from
the calculation as they are, so a frequency job run with an implicit
solvent model (PCM, SMD, CPCM) gives KIEs in that solvent. Kinisot adds no
solvent correction of its own.

## Checking the structures

**My transition structure has two imaginary frequencies.** Kinisot warns
and takes the larger one as the reaction coordinate. The other is discarded
along with the overall translations and rotations, so one of those (a mode
near zero) is counted as a vibration in its place. With `--project` (the
default for xTB and machine-learned potentials) nothing is discarded that
way, so the second mode stays among the vibrations and Kinisot stops with
an error. Either way, re-optimize the transition structure. If the second
mode is a genuine low-frequency torsion smaller than the cutoff (50i cm⁻¹
by default), nothing is wrong.

**My reactant has a small imaginary frequency.** Kinisot stops: reactants
and products must be minima. If the mode is a numerical artefact (a few
cm⁻¹ on a floppy molecule), raise `--imag-cutoff` above its magnitude.
Otherwise re-optimize.

**Why do my frequencies differ slightly from the ones Gaussian prints?**
Gaussian removes the overall translations and rotations before it computes
frequencies. Kinisot by default computes all of them and drops the six (five
for a linear molecule) lowest. On a converged geometry the vibrational
frequencies agree to better than 0.05 cm⁻¹ (this is tested), and the
discarded modes are the values on Gaussian's "Low frequencies" line.

**Should I use `--project`?** It removes translations and rotations exactly,
as Gaussian does. For Gaussian and ORCA outputs on converged geometries it
makes no practical difference (below 10⁻⁵ in the KIE). It is on by default
for xTB and machine-learned potentials, whose finite-difference Hessians
need it. Use it for any file whose "discarded" modes are not all within a
few tens of cm⁻¹ of zero.

## Choices that change the numbers

**Which number do I report?** `corr-KIE`, the KIE with the tunnelling
correction. The `KIE` column next to it is the value without tunnelling.

**Which tunnelling correction should I use?** Bell's infinite parabola is
the default. It is what Meyer, DelMonte and Singleton found sufficient for
heavy-atom KIEs in the Claisen rearrangement. `--tunneling wigner` is its
first-order approximation, and `--tunneling none` gives the value without
tunnelling. Bell's correction is refused below the crossover temperature,
where neither one-dimensional model is trustworthy. For primary hydrogen
KIEs with a large tunnelling contribution, any one-dimensional correction is
only approximate.

**Which scaling factor is applied?** The ZPE factor from the Truhlar
database for the level of theory Kinisot detects in the files, or 1.0 if it
is not listed. Give `-s` to override. Only the ZPE, EXC and tunnelling terms
depend on it; the V-ratio and TRPF are ratios of frequencies and cancel it.

## Comparing with experiment

**My measurements are relative to an internal standard.** Give the
standard's position with `--reference`. Kinisot then also prints the KIE
divided by the standard's KIE (the `relative to` line), which is the number
to compare with natural-abundance measurements.

**One NMR signal comes from two positions.** Examples are the two ortho
carbons of a phenyl ring, the three hydrogens of a methyl group, and the
two oxygens of a nitro group. When the positions are equivalent in the
experiment but distinct in the transition structure, compute each
placement of the label and average them.
- **Exactly:** in Python, pass one result per placement to
  `kinisot.equivalent_positions()`.
- **From the command line:** list the placements as channels with equal
  shares in a job file ([job_files.md](job_files.md)).

When the positions are equivalent in the reactant, both give the harmonic
mean of the separate KIEs, not their arithmetic mean.

**My computed KIEs disagree with experiment.** Look at the transition
structure first: heavy-atom KIEs follow its geometry closely, so a
disagreement usually means a different transition structure or a
different mechanism. Then consider conformers, a step that is not
rate-limiting alone, or competing pathways. The README's
[Comparing with experiment](../README.md#comparing-with-experiment) section
has the accuracy to expect.

## More than one structure

**Can I use several conformers?** Yes, from 2.6. Give them all after one
flag: `kinisot --rct gs_*.out --ts ts_*.out --iso 5`. They must have the
same atoms in the same order. Kinisot weights them by the free energies of
the light isotopologue. It uses quasi-harmonic free energies by default
(`--weights`), or your own values (`--energies`). The result is not the
Boltzmann average of the pairwise KIEs
([theory, section 7](theory.md#7-conformer-ensembles)).

The output lists:
- each conformer's share and its own KIE;
- the KIE of the lowest conformers alone;
- how much the result moves when each free energy moves by 0.5 kcal/mol.

[examples/conformers](../examples/conformers/README.md) works through a
case.

**What if an intermediate can return to the reactant?** Then no single
transition structure commits the substrate, and the KIE is a weighted mean
of the KIEs of the steps. Set this up as a series in a job file
([job_files.md](job_files.md)), with the steps' free energies or a
commitment factor ([theory, section 7a](theory.md#7a-transition-structures-in-series)).
Parallel pathways with their own labels, such as two enantiomers reacting
through diastereomeric transition structures, are channels
([section 7b](theory.md#7b-parallel-channels)).

## Other programs and methods

**Why not just run Gaussian with `freq=(readisotopes)` for each
isotopologue?** You can, and the KIEs agree to the last printed digit once
the same scaling is used. Kinisot avoids re-running the frequency job for
every label, which matters when you scan all positions of a molecule at
several temperatures.

**Can I use a machine-learned potential or xTB?** Yes. Optimize the
reactant and locate the transition structure with the same method, then run
`kinisot --rct rct.xyz --ts ts.xyz --iso 4 --calc mace_mp:medium` (or
`--calc xtb` for GFN2-xTB). There is no scaling factor for a potential, and
its accuracy for the curvature at the transition structure decides the
KIE, so compare against a DFT reference where you can.

Check that the saddle point is the reaction you mean. None of four MACE
foundation models (MACE-OFF23 small, medium and large, and MACE-MP-0) has
the concerted Claisen transition structure: their saddle points describe
C–O cleavage or ring closure instead. Seven others have it: UMA
(`--calc uma`), MACE-OMOL-0, SevenNet-Omni, AIMNet2, AIMNet2-rxn and two
ORB models. With UMA, MACE-OMOL-0 and SevenNet-Omni the KIEs match
experiment as well as B3LYP's; with the others they miss by 0.003–0.004
on average.

Check the Hessian as well. The ORB models' energies change when the
molecule is rotated, so their Hessians depend on its orientation, and their
transition structures are too soft. An EQE between two equivalent atoms
(exactly 1 by symmetry) and the energy of the rotated molecule are quick
tests.

[examples/mlip_claisen](../examples/mlip_claisen/README.md) compares
these potentials with GFN2-xTB and B3LYP for the Claisen rearrangement and gives the details, and
`scripts/make_claisen_structures.py` shows the whole workflow, including
the transition-structure search.

**How do I get the numbers into a script?** `kinisot ... --json run.json`
or `--csv runs.csv` (one row per run), or call `kinisot.compute_kie()` from
Python and use the result it returns (`r.kie_tunnel`, `r.to_dict()`; see
`examples/api_example.py`).
