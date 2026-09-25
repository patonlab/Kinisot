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
2010 values in `Kinisot.py`). The scaling factor multiplies every ν,
including the imaginary one.

**External modes.** Translations and rotations are not projected out.
The eigenvalues are sorted and the six lowest (five for a linear molecule)
are discarded; if the lowest mode is imaginary beyond the cutoff
(default 50 cm⁻¹) it is the reaction coordinate and the six (five) next
lowest are discarded instead. On tightly converged geometries the discarded
modes are within a few cm⁻¹ of zero (Gaussian's "Low frequencies" line, see
`tests/test_frequencies.py`); Kinisot prints them so that a genuine
low-frequency vibration that has slipped below a rotational residual can be
spotted. Eckart projection is planned (implementation plan, Phase 7).

**Masses.** The light isotopologue uses the program's masses (pure
most-abundant isotopes in Gaussian: ¹H 1.00783, ¹²C 12.00000, ¹⁶O
15.99491). Substitution replaces them with ²H 2.0141, ¹³C 13.00335 and
¹⁷O 16.9991 (`kinisot/isotopes.py`).

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
first order. Kinisot prints both the uncorrected and the corrected KIE;
report the one that matches the model used in the work you compare with.

## 5. Scaling factors

When `-s` is not given, the level of theory and basis set are read from the
Gaussian archive entry and matched (case-insensitively, ignoring hyphens and
Gaussian's R/U/RO prefix) against the Truhlar group database (version 3b2,
Alecu, Zheng, Zhao, Truhlar, J. Chem. Theory Comput. 2010, 6, 2872). The
**ZPE** scaling factor is used because the zero-point term dominates
isotope effects. The factor multiplies all frequencies of all species, so
only the temperature-dependent terms (ZPE, EXC, tunnelling) are affected;
TRPF and V-ratio are ratios of frequencies and cancel it.

## 6. Constants (CODATA 2010)

| Constant | Value |
| --- | --- |
| h | 6.62606957 × 10⁻³⁴ J s |
| k | 1.3806488 × 10⁻²³ J K⁻¹ |
| c | 2.99792458 × 10¹⁰ cm s⁻¹ |
| E_h | 4.35974434 × 10⁻¹⁸ J |
| a_0 | 5.2917721092 × 10⁻¹¹ m |
| u | 1.660538921 × 10⁻²⁷ kg |
