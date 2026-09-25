# FAQ

**Why do my frequencies differ slightly from the ones Gaussian prints?**
Gaussian projects translations and rotations out of the Hessian before
diagonalizing; Kinisot diagonalizes the raw mass-weighted Hessian and drops
the six (five) lowest modes. On a converged geometry the vibrational
frequencies agree to better than 0.05 cm⁻¹ (this is tested), and the
discarded modes are the values on Gaussian's "Low frequencies" line.

**Why not just run Gaussian with `freq=(readisotopes)` for each
isotopologue?** You can, and the KIEs agree to the last printed digit once
the same scaling is used. Kinisot avoids re-running the frequency job for
every label, which matters when scanning all positions of a molecule at
several temperatures.

**What does `--iso 0` mean?** No substituted atom in that file. It is
needed when a side of the reaction is given as several files and only one
of them carries the label. `0` cannot be combined with atom numbers.

**The atom numbering differs between my reactant and my TS.** Give one
`--iso` per file in the order of the `--rct` flags followed by `--ts` or
`--prd`. Kinisot checks that both sides substitute the same elements.

**My TS has two imaginary frequencies.** Kinisot warns and treats the
larger one as the reaction coordinate; the other is discarded with the
external modes, which leaves one rotational residual in the vibrational
product. Re-optimize the transition structure. If the second mode is a
genuine low-frequency torsion below the cutoff, nothing is wrong.

**My reactant has a small imaginary frequency.** Kinisot stops: reactants
and products must be minima. If the mode is a numerical artefact (a few
cm⁻¹ on a floppy molecule), raise `--imag-cutoff` above its magnitude;
otherwise re-optimize.

**Which tunnelling correction should I use?** Bell's infinite parabola is
the default and what the Claisen study of Meyer, DelMonte and Singleton
found sufficient for heavy-atom KIEs; `--tunneling wigner` is its
first-order expansion and `--tunneling none` gives the semiclassical value.
Bell is refused below the crossover temperature, where neither
one-dimensional model is trustworthy.

**Which scaling factor is applied?** The ZPE factor from the Truhlar
database for the detected level of theory, or 1.0 if it is not listed.
Give `-s` to override. Only ZPE, EXC and the tunnelling correction depend
on it; the V-ratio and TRPF are ratios of frequencies and cancel it.

**Which isotope does a bare atom number mean?** The usual heavy label:
²H, ¹³C, ¹⁵N, ¹⁸O, ³⁴S, ³⁷Cl, ⁸¹Br, ²⁹Si. Anything else is asked for
explicitly (`--iso 3:17O`, `7:D`, `7:T`, `5:14C`, or a mass `5:13.5`).
Kinisot ≤ 2.3 substituted ¹⁷O for oxygen (the isotope measured by ¹⁷O NMR
in the Claisen study); 2.4 uses ¹⁸O and prints a note once when it sees a
bare oxygen index.

**Should I use `--project`?** For Gaussian and ORCA analytic Hessians on
converged geometries it makes no practical difference (below 10⁻⁶ in the
KIE). It is on by default for Hessians from ASE calculators, whose
finite-difference residuals need it; use it for any file whose "discarded"
modes are not all within a few tens of cm⁻¹ of zero.

**Can I use a machine-learned potential?** Yes: optimize the reactant and
locate the transition structure with the potential, then
`kinisot --rct rct.xyz --ts ts.xyz --iso 4 --calc mace_mp:medium`. See
`examples/mlip_claisen/`. There is no scaling factor for a potential, and
its accuracy for the curvature at the transition structure decides the
KIE, so compare against a DFT reference where you can.

**Can I use several conformers?** Not in one run. Compute the KIE for
each reactant/TS conformer pair and Boltzmann-average the rate constants,
or use the lowest-energy pair.

**Does the temperature have to match the Gaussian job?** No. Kinisot
evaluates the partition functions at `-t`; the Hessian does not depend on
temperature.

**How do I get the numbers into a script?** `kinisot ... --json run.json`
or `--csv runs.csv` (one row per run), or call `kinisot.compute_kie()` from
Python and use the returned `IsotopeEffect` (`r.kie_tunnel`, `r.to_dict()`,
see `examples/api_example.py`).
