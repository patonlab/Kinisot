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

is the equilibrium constant of R(heavy) + P(light) ⇌ R(light) + P(heavy),
printed in the `KIE` column of an `EQE @` line. Symmetry numbers are not
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

## 6. Constants (CODATA 2018)

| Constant | Value |
| --- | --- |
| h | 6.62607015 × 10⁻³⁴ J s |
| k | 1.380649 × 10⁻²³ J K⁻¹ |
| c | 2.99792458 × 10¹⁰ cm s⁻¹ |
| E_h | 4.3597447222 × 10⁻¹⁸ J |
| a_0 | 5.29177210903 × 10⁻¹¹ m |
| u | 1.66053906660 × 10⁻²⁷ kg |
