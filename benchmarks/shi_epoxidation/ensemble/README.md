# Shi epoxidation: the other 17 transition structures

Gaussian inputs for the 17 transition structures of Singleton and Wang
(J. Am. Chem. Soc. 2005, 127, 6679) other than TS 10, which is
`../ts10.gjf`. They are the data for the conformer-ensemble validation
planned in IMPLEMENTATION_PLAN.md, Phase 10. Each input is a B3LYP/6-31G(d)
frequency job (`int=finegrid`, no optimization) at the SI geometry.

- **Names** follow the SI: A (TS 10 in the paper), B (11), EA (12), DA (13),
  EB (14), D (15), BD (16), DD (17), then AA, G, CA, H, E, C, AD, CD, AB and
  CB.
- **Titles** carry each structure's SI electronic energy and zero-point
  energy. A job at the correct geometry reproduces both, as `../ts10.log`
  does for TS A (E to 1 × 10⁻⁸ hartree, ZPE exactly).

## Atom numbering

The SI lists each structure in its own atom order. Here every file uses
TS A's order (`../ts10.gjf`), so an isotope label means the same atom in
every structure, as Kinisot's ensemble treatment will require.

The map was built by matching bonding graphs with TS A. Three constraints
fix the assignments that bonding alone leaves open:

- **Stereocentres.** Every catalyst stereocentre keeps the sign of its
  chirality volume, so the diastereotopic methyl groups and the two dioxirane
  oxygens keep their identities.
- **Phenyl carbons.** The ortho and meta carbons are assigned syn or anti to
  C-β.
- **Hydrogens.** Each hydrogen follows its parent atom, ordered by dihedral
  angle.

For all 17 structures, the map preserves every bond, all six stereocentres
and both syn/anti assignments.

In this numbering, O26 is the oxygen transferred in every structure. The
alkene is C7 (C-α) and C8 (C-β), the methyl carbon is C9, ipso is C2, the
ortho carbons are C1 and C3, the meta carbons C4 and C6, and para is C5.

`si_atom_order.json` gives, for each structure, the SI atom number of
every atom here, together with the SI energy and zero-point energy.

## A misprint in the SI (TS AB)

The SI's coordinates for TS AB contain a misprint. Four minus signs are
missing and one integer digit is wrong; the printed PDF shows the same
values, so the text extraction is not at fault.

The SI prints the AB coordinates with a constant offset added, so the last
four digits of each value show its sign: …9571/…0429 in x, …9857/…0143 in
y and …9905/…0095 in z. Four values carry a negative-pattern tail without a
minus sign:

- atom 15 (C), y and z;
- atoms 31 and 32 (H), z.

With the signs restored, the y of atom 15 is still off by exactly 1.0.
Printed as `0.1844600143`, it is corrected to `-1.1844600143`. That is the
value its three hydrogens place 1.53 Å from its neighbouring carbon; the
other two coordinates agree with the hydrogens to 2 × 10⁻⁴ Å. These atom
numbers are the SI's.

`ts_AB.gjf` contains the corrected geometry. If the repair is exact, its
job reproduces the SI energy of −1343.48585457 hartree.
