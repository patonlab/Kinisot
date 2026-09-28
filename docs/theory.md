# Theory: what Kinisot computes

This page states the equations exactly as implemented in `kinisot/Kinisot.py`
and `kinisot/Hess_to_Freq.py`, with Kinisot's sign and naming conventions.
Subscript L is the light (unsubstituted) isotopologue and H the heavy one;
‡ marks the transition structure.

## 1. Harmonic frequencies from the Hessian

Kinisot does not use the frequencies printed by the quantum chemistry
program. It reads the Cartesian force-constant matrix **H** (Hartree/Bohr²)
and the per-atom masses, and for each isotopologue builds the mass-weighted
Hessian

    H'_ij = H_ij / sqrt(m_i m_j)

with the substituted atoms carrying their heavy-isotope masses. The
eigenvalues λ of **H'** give harmonic wavenumbers

    ν = sign(λ) sqrt(|λ| · E_h / (a_0² u)) / (2π c)

in cm⁻¹, negative for imaginary modes, where E_h, a_0, u and c are the
Hartree energy, Bohr radius, atomic mass unit and speed of light (CODATA
2018 values in `kinisot/thermo.py`). The scaling factor multiplies every ν,
including the imaginary one.

**External modes.** By default translations and rotations are not
projected out: the eigenvalues are sorted and the six lowest (five for a
linear molecule) are discarded; if the lowest mode is imaginary beyond the
cutoff (default 50 cm⁻¹) it is the reaction coordinate and the six (five)
next lowest are discarded instead. On tightly converged geometries the
discarded modes are within a few cm⁻¹ of zero (Gaussian's "Low
frequencies" line, see `tests/test_frequencies.py`); Kinisot prints them so
that a genuine low-frequency vibration that has slipped below a rotational
residual can be spotted.

With `--project` (`project=True`) the six (five) external directions are
built in mass-weighted Cartesians from the geometry (translations
√m_i e_k, rotations √m_i e_k × (r_i − r_com)), orthonormalized, and
projected out: H'' = (1 − QQᵀ) H' (1 − QQᵀ). The external eigenvalues are
then zero to numerical precision and are removed by magnitude; the
reaction coordinate is the most negative remaining eigenvalue. On the
bundled Gaussian examples the two treatments agree to 4 × 10⁻⁷ in the KIE
and the projected frequencies match Gaussian's printed ones to 0.01 cm⁻¹;
projection matters for Hessians with sizeable translational/rotational
residuals (finite differences, machine-learned potentials).

**Masses.** Both isotopologues are built from the isotope table in
`kinisot/isotope_data.py` (AME 2020 masses via the `periodictable`
package): the light one from the most abundant isotope of every element
(¹H 1.0078250, ¹²C 12, ¹⁶O 15.9949146, the convention Gaussian uses) and
the heavy one with the requested atoms replaced (²H 2.0141018, ¹³C
13.0033548, ¹⁸O 17.9991596, ...). The masses reported by the program are
used only to check that no isotope was substituted in its input.

## 2. Bigeleisen–Mayer reduced partition function ratios

With u = h c ν / k T for each of the 3N−6 (3N−5, 3N−7 for a TS) real
vibrations, the reduced isotopic partition function ratio of a species is

    (s/s')f = ∏_i (u_H,i / u_L,i) · exp[(u_L,i − u_H,i)/2] · (1 − e^{−u_L,i}) / (1 − e^{−u_H,i})

Kinisot accumulates the three factors in logarithmic form for each species:

| Term | Code | Definition (per species, light/heavy) |
| --- | --- | --- |
| Teller–Redlich product factor | `PF` | Σ ln(h c ν / k) — TRPF_species = exp(PF_H − PF_L) = ∏ ν_H/ν_L |
| Zero-point energy | `ZPE` | Σ u/2 — ZPE_species = exp(ZPE_L − ZPE_H) = ∏ e^{(u_L − u_H)/2} |
| Excitation | `EXC` | Σ ln(1 − e^{−u}) — EXC_species = exp(EXC_L − EXC_H) |

These are the three numbers printed on each species row of the results
table; their product is (s/s')f for that species. When a side of the
reaction is given as several files (bimolecular reactions), the
logarithmic terms of the files are added, i.e. the partition functions are
multiplied.

## 3. Kinetic and equilibrium isotope effects

For a KIE (reactant R, transition structure ‡):

    KIE = (ν‡_L / ν‡_H) · (s/s')f_R / (s/s')f_‡

and the results line prints the ratios factor by factor:

| Column | Value |
| --- | --- |
| `V-ratio` | ν‡_L / ν‡_H, the ratio of the imaginary frequencies |
| `ZPE` | ZPE_R / ZPE_‡ |
| `EXC` | EXC_R / EXC_‡ |
| `TRPF` | TRPF_R / TRPF_‡ |
| `KIE` | V-ratio × ZPE × EXC × TRPF (semiclassical, no tunnelling) |
| `1D-tunn` | Bell tunnelling correction Q_t,L / Q_t,H (Section 4) |
| `corr-KIE` | KIE × 1D-tunn |

For an equilibrium isotope effect (reactant R, product P) there is no
reaction coordinate: `V-ratio` is empty, `1D-tunn` is 1, and

    EQE = (s/s')f_R / (s/s')f_P

is the equilibrium constant of R(light) + P(heavy) ⇌ R(heavy) + P(light),
that is K_L/K_H for R ⇌ P. It is printed in the `KIE` column of an
`EQE @` line. Above 1, the heavy isotope accumulates in R, the side with
the stiffer vibrations. Symmetry numbers are not
included: (s/s')f is by definition the *reduced* ratio. If a substitution
changes the symmetry number of a species (e.g. CH₃ → CH₂D), multiply by
the ratio s/s' yourself.

## 4. Tunnelling: Bell's infinite parabola

The one-dimensional Bell correction for a parabolic barrier with imaginary
frequency ν‡ is

    Q_t = (u‡/2) / sin(u‡/2),   u‡ = h c |ν‡| / k T

and the KIE correction is the ratio for the two isotopologues,

    1D-tunn = Q_t,L / Q_t,H = (ν‡_L / ν‡_H) · sin(u‡_H/2) / sin(u‡_L/2)

The formula diverges as u‡ → 2π (T → h c |ν‡| / 2π k, the crossover
temperature) and is meaningless below it: for a 1000i cm⁻¹ mode that is
229 K. Above that limit it reproduces the Wigner correction 1 + u‡²/24 to
first order; `--tunneling wigner` uses that expansion and `--tunneling none`
turns the correction off.

**Skodje–Truhlar** (`--tunneling skodje`) adds the barrier height V to the
parabolic model (Skodje, Truhlar, J. Phys. Chem. 1981, 85, 624). With
α = 2π / (h ν‡) and β = 1 / (k T):

    β ≤ α:  κ = (βπ/α) / sin(βπ/α) − β / (α − β) · exp[(β − α) V]
    β > α:  κ = β / (β − α) · [exp((β − α) V) − 1]

and the correction is κ_L / κ_H. V is the electronic barrier measured from
the higher of reactant and product (for a KIE Kinisot uses E(TS) − ΣE(R)
from the energies in the files, or `--barrier` in kcal/mol). For large V
the first branch reduces to Bell; the second branch stays finite below the
crossover temperature. Kinisot prints both the uncorrected and the
corrected KIE; report the one that matches the model used in the work you
compare with.

## 4a. Reference isotopologue

Natural-abundance NMR experiments report KIEs relative to a position
assumed to have no isotope effect. `--reference ATOMS` computes that
isotopologue too and reports KIE / KIE_ref (`kie_relative`,
`kie_tunnel_relative`), so computed and measured numbers are on the same
footing.

## 5. Scaling factors

When `-s` is not given, the level of theory and basis set are read from the
output (Gaussian archive entry or ORCA `!` line) and matched, after
GoodVibes' canonicalization of program-specific spellings and stripping
Gaussian's R/U/RO prefix, against the Truhlar group database (version 5,
shipped with GoodVibes; Alecu, Zheng, Zhao, Truhlar, J. Chem. Theory
Comput. 2010, 6, 2872). The **ZPE** scaling factor is used by default
because the zero-point term dominates isotope effects; `--scale-type harm`
or `fund` selects the harmonic or fundamental factor. The factor multiplies all frequencies of all species, so
only the temperature-dependent terms (ZPE, EXC, tunnelling) are affected;
TRPF and V-ratio are ratios of frequencies and cancel it.

## 6. Bigeleisen–Mayer versus free-energy differences

A KIE can also be computed as exp(ΔΔG‡/RT) from the free energies a
quantum chemistry program prints for each isotopologue (the "Sum of
electronic and thermal Free Energies" of Gaussian, for instance). There the
translational and rotational partition functions come from the masses and
moments of inertia. Bigeleisen and Mayer instead replace their isotope
ratios by the product of the vibrational frequency ratios (the
Teller–Redlich product rule, the TRPF and V-ratio terms of section 3).

For exact harmonic frequencies at an exact stationary point the two are
identical. Real frequencies obey the product rule only approximately, and
the two routes then differ by exactly the violation:

    ln(KIE_FE / KIE_BM) = −(δ_R − δ‡),   δ = ln ∏ ν_H/ν_L − ln[(M_H/M_L)^{3/2} (I_H/I_L)^{1/2} ∏ (m_L/m_H)^{3/2}]

with the product over all 3N−6 modes (the imaginary one included) and I the
product of the principal moments of inertia. Rotational symmetry numbers
are left out of both routes here, as Kinisot leaves them out (section 3).
A program that includes them adds ln[(σ_H/σ_L)‡ / (σ_H/σ_L)_R] to the right-hand
side; the term vanishes when no substitution changes a symmetry number, or
when it changes those of the reactant and the transition structure by the
same ratio, as in the Claisen example.

**Translational entropy.** The mass term (M_H/M_L)^{3/2} in δ is the
translational partition-function ratio. Bigeleisen–Mayer does not use it. An
unlabelled species' reduced ratio is exactly 1, so it can simply be left out
of the input: the chloride of an SN2 reaction, for instance. In the
free-energy route an unlabelled species' own free energy cancels between
isotopologues as well. The translational terms of the labelled species do
not: they cancel only when the reactant and the transition structure have
the same total mass.

A free-energy or enthalpy–entropy treatment that drops the translational
term is therefore wrong by (3/2) ln[(M_H/M_L)_R / (M_H/M_L)‡], with M the
total masses of the labelled reactant and of the transition structure. The
term is not zero whenever an unlabelled partner adds its mass to the
transition structure, as the chloride does in Cl⁻ + CH₃Br: the chloride's
own terms cancel, but its mass is part of the transition structure's.
PyQuiverHS's enthalpy–entropy terms are an example: they contain
vibrational and rotational entropy only. The error is 1.3% for
Cl⁻ + CH₃Br/CD₃Br. Grazioli et al. report an enthalpy–entropy KIE of 0.877
at 300 K where Bigeleisen–Mayer gives 0.888 (their Table 7;
benchmarks/sn2_chloride_methyl_bromide). Free energies from the same
frequencies, without the translational term, give exactly 0.877; with it
they give 0.8882, and the term above is the factor 1.0128 between them.
When the same atoms appear on both sides (unimolecular reactions,
conformational equilibria) the masses match, the translational terms cancel
and the omission is harmless.

The Bigeleisen–Mayer form is the robust one because of what happens to a
soft vibration. For a mode with u = hcν/kT ≪ 1 its TRPF, ZPE and EXC
factors tend to (ν_H/ν_L) · 1 · (ν_L/ν_H) = 1, so the mode drops out however
poorly its frequency is computed. In the free-energy route the same mode
contributes ln(ν_H/ν_L) in full. That contribution cancels only against
rotational terms computed from moments of inertia, so an error in a soft
mode's isotope shift goes straight into the KIE. Beno, Houk and Singleton
reported for the Diels–Alder reaction of isoprene with maleic anhydride
that errors in the low frequencies of a floppy transition structure cause
little error in Bigeleisen–Mayer KIEs (J. Am. Chem. Soc. 1996, 118, 9984,
[doi:10.1021/ja9615278](https://doi.org/10.1021/ja9615278)). Hirschi,
Takeya, Hang and Singleton give the same argument for constrained,
non-stationary structures (J. Am. Chem. Soc. 2009, 131, 2397,
[doi:10.1021/ja8088636](https://doi.org/10.1021/ja8088636)).

`scripts/compare_free_energy_route.py` computes both routes from the same
projected frequencies of the Claisen example, B3LYP/6-31G(d), at 393 K,
unscaled and without tunnelling. `tests/test_free_energy_route.py` checks
the numbers quoted here.

| Site | Kinisot (BM) | FE | FE, 6 decimals | BM, +0.1 cm⁻¹ | FE, +0.1 cm⁻¹ | BM, quasi-harmonic | FE, quasi-harmonic |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C1 | 1.01296 | 1.01269 | 1.01294 | 1.01296 | 1.01123 | 1.01292 | 1.00418 |
| C2 | 1.00192 | 1.00194 | 1.00161 | 1.00191 | 1.00051 | 1.00192 | 1.00184 |
| O3 (¹⁸O) | 1.03768 | 1.03737 | 1.03765 | 1.03768 | 1.03587 | 1.03758 | 1.01845 |
| C4 | 1.03106 | 1.03097 | 1.03100 | 1.03106 | 1.02950 | 1.03104 | 1.02697 |
| C5 | 1.00192 | 1.00192 | 1.00161 | 1.00191 | 1.00049 | 1.00192 | 1.00182 |
| C6 | 1.01511 | 1.01555 | 1.01620 | 1.01510 | 1.01409 | 1.01506 | 1.00726 |
| H7,H8 (²H₂) | 0.95064 | 0.94981 | 0.95064 | 0.95063 | 0.94841 | 0.95047 | 0.91856 |

The columns:

- **FE**: exact arithmetic. The files are converged (maximum force below
  10⁻⁴ hartree/bohr), yet the product rule is violated by up to 1 × 10⁻³ in
  ln (H7,H8 in the reactant). That is equivalent to an error of 0.07 cm⁻¹
  in the reactant's 70.5 cm⁻¹ torsion, and it moves the free-energy KIEs by
  up to 8 × 10⁻⁴.
- **FE, 6 decimals**: each free energy rounded to the six decimals Gaussian
  prints. At 393 K an error of 10⁻⁶ hartree in one of the four free
  energies changes a KIE by 8 × 10⁻⁴. Rounding shifts C6 by 6.5 × 10⁻⁴, and
  the worst case, four values each half a unit off, is 1.6 × 10⁻³. That is
  as large as the uncertainty of many natural-abundance ¹³C measurements.
- **+0.1 cm⁻¹**: the heavy isotopologue's 70.5 cm⁻¹ reactant torsion
  shifted by 0.1 cm⁻¹. The Bigeleisen–Mayer KIEs change by at most
  8 × 10⁻⁶; the free-energy KIEs change by 1.4 × 10⁻³ to 1.5 × 10⁻³, about
  180 times more.
- **Quasi-harmonic**: every real mode below 100 cm⁻¹ raised to 100 cm⁻¹
  (Truhlar's quasi-harmonic free energies; Grimme's entropy interpolation
  has the same effect). Bigeleisen–Mayer changes by at most 2 × 10⁻⁴. The
  free-energy ¹⁸O KIE drops from 1.037 to 1.018, because the soft mode's
  isotope shift leaves the vibrational term but not the rotational one.
  Quasi-harmonic free energies, such as those GoodVibes prints for
  thermochemistry, must not be used for isotope effects.

Tunnelling matters as much as the choice of equation. Meyer, DelMonte and
Singleton found that a one-dimensional tunnelling correction improves
heavy-atom KIE predictions. With it, the difference from experiment fell to
about the experimental uncertainty in the reactions they studied (J. Am.
Chem. Soc. 1999, 121, 10865,
[doi:10.1021/ja992372h](https://doi.org/10.1021/ja992372h)). Kinisot
applies Bell's correction by default.

## 7. Conformer ensembles

With transition-state theory and conformers in fast equilibrium
(Curtin–Hammett), isotopologue X reacts with

    k^X = (k T / h) Σ_j κ_j^X Q‡_j^X exp(−E‡_j / kT) / Σ_i Q_i^X exp(−E_i / kT)

where i runs over the reactant conformers and j over the
transition-structure conformers. Every conformer of a species has the same
atoms, so the Teller–Redlich mass factor is the same for all of them and
cancels, and

    KIE = Σ_i x_i ρ_i / Σ_j y_j ρ‡_j

- ρ_i is (s/s')f of reactant conformer i (Section 2).
- ρ‡_j = (s/s')f‡_j · (ν‡_H / ν‡_L)_j · (κ_H / κ_L)_j, so that ρ_i / ρ‡_j is
  the ordinary KIE of that pair of conformers.
- x_i is the population of reactant conformer i, ∝ g_i exp(−G_i / RT).
- y_j is transition structure j's share of the rate, ∝ g_j κ_L,j exp(−G‡_j / RT).
- g is a degeneracy, for example 2 for a conformer whose mirror image is
  not in the list.

The weights are those of the light isotopologue, and that is exact, not an
approximation:

    Σ_i Q_i^H e^{−E_i/kT} / Σ_i Q_i^L e^{−E_i/kT} = Σ_i x_i (Q_i^H / Q_i^L)

Three consequences follow:

1. **Only free energies within each ensemble matter.** The gap between the
   reactant and transition-structure ensembles cancels.
2. **The ensemble KIE is a ratio of means.** It is neither the Boltzmann
   average of the pairwise KIEs nor the KIE of the lowest pair. A minor
   transition structure counts in proportion to its share of the rate.
3. **Several species multiply.** With two reactants, the reactant side is
   the product of the two ensemble means. For an EQE, the product
   conformers take the place of the transition structures, with κ = 1 and
   no ν‡ ratio.

**Weights.** The isotope ratios stay harmonic (Section 6), but the weights
are free energies of the light isotopologue, for which a quasi-harmonic
treatment of soft modes is appropriate. `--weights` chooses:

| Weights | G |
| --- | --- |
| `qrrho` (default) | E + ZPE + thermal vibrational energy − T(S_vib + S_rot), with Grimme's interpolation of S_vib towards a free rotor below 100 cm⁻¹ (via GoodVibes) |
| `rrho` | the same with harmonic S_vib |
| `user` | free energies you give (`--energies`, or `Conformers(..., free_energies=...)`) |
| `lowest` | the lowest conformer of every species alone |
| `equal` | the degeneracies alone |

Translation, electronic entropy and the rotational energy are the same for
every conformer of a species and are left out. So is the rotational
symmetry number: when conformers of one species differ in it, give each
conformer 1/σ (relative to the others) as its degeneracy. The frequencies
are Kinisot's scaled ones, projected or not as for the isotope effect;
with `--project`, the relative free energies agree with GoodVibes' from
Gaussian's printed frequencies to within 0.01 kcal/mol.

**Diagnostics.** The result also gives the KIE of the lowest conformers
alone, the effective number of transition-structure conformers 1 / Σ y_j²,
and the range of the KIE when each free energy moves by ±0.5 kcal/mol in
turn (`--weight-uncertainty`).

**Equivalent positions.** Positions made equivalent by fast motion (the
three hydrogens of a rotating methyl group, the two oxygens of a nitro
group, the ortho carbons of a spinning phenyl ring) are an ensemble of the
placements of the label, with equal weights. The isotope effect is the mean
of ρ over the placements divided by the mean of ρ‡. When the positions
differ on one side only, as in a transition structure that attacks one of
two equivalent oxygens, this is the harmonic mean of the separate KIEs, not
their arithmetic or geometric mean. `kinisot.equivalent_positions` computes
it from one result per placement, and the benchmark runner uses it.

**Limits.** The conformers must interconvert faster than they react.
Transition structures in series (an intermediate that can return) and
parallel channels to different products need different formulas, planned
in IMPLEMENTATION_PLAN.md, Phase 10.

## 8. Constants (CODATA 2018)

| Constant | Value |
| --- | --- |
| h | 6.62607015 × 10⁻³⁴ J s |
| k | 1.380649 × 10⁻²³ J K⁻¹ |
| c | 2.99792458 × 10¹⁰ cm s⁻¹ |
| E_h | 4.3597447222 × 10⁻¹⁸ J |
| a_0 | 5.29177210903 × 10⁻¹¹ m |
| u | 1.66053906660 × 10⁻²⁷ kg |
